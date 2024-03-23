#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import yaml
import xarray as xr
import numpy as np

from src.inversion_scripts.point_sources import get_point_source_coordinates
from src.inversion_scripts.imi_preview import (
    estimate_averaging_kernel,
    map_sensitivities_to_sv,
)

# clustering
from sklearn.cluster import KMeans, MiniBatchKMeans


def cluster_data_kmeans(data, num_clusters, mini_batch=False):
    """
    Description:
        Given the sensitivities cluster grid cells into the provided
        number of state vector elements
    arguments:
        data       [][]dataarray : xarrray sensitivity data
        num_clusters         int : number of labels to assign data to
        mini_batch          bool : optionally use MiniBatchKmeans to
                                   speed up clustering algorithm
    Returns:         [][]ndarray : labeled data
    """
    # Get the latitude and longitude coordinates as separate arrays
    latitudes = data.coords["lat"].values
    longitudes = data.coords["lon"].values

    # Get the sensitivity values as a 1D array
    Z = data.values.flatten()
    # labels shape for later
    labels = np.zeros(Z.shape)
    valid_indices = ~np.isnan(Z)

    # Flatten the latitude and longitude arrays into a 2D grid
    # only keeping valid indices
    X, Y = np.meshgrid(longitudes, latitudes)
    X = X.flatten()[valid_indices]
    Y = Y.flatten()[valid_indices]
    Z = Z[valid_indices]

    # Stack the X, Y, and Z arrays to create a (n_samples, n_features) array
    features = np.column_stack((X, Y, Z))

    # Cluster the features using KMeans
    # Mini-Batch k-means is much faster, but with less accuracy
    if mini_batch:
        kmeans = MiniBatchKMeans(n_clusters=num_clusters, random_state=0)
    else:
        kmeans = KMeans(n_clusters=num_clusters, random_state=0)

    cluster_labels = kmeans.fit_predict(features)

    # fill labels on corresponding valid indices of label array
    # add 1 to labels so they start with 1
    labels[valid_indices] = cluster_labels + 1

    # reconstruct 2D grid
    cluster_labels = labels.reshape(data.shape)

    return cluster_labels


def get_highest_labels(labels, sensitivities, n):
    """
    Description:
        Get the n labels that have the highest avg sensitivity
        per grid cell. Returned in descending order.
    arguments:
        labels            [][]ndarray : state vector labels
        sensitivities   [][]dataarray : xarrray sensitivity data
        n                         int : number of labels to return
    Returns:      [] : list of n labels with highest sensitivities
    """
    sensitivity_dict = {}
    max_label = int(np.nanmax(labels))

    # calculate avg dofs per gridcell for each cluster
    for i in range(1, max_label + 1):
        indices = np.where(labels == i)
        cluster_sensitivities = sensitivities[indices]
        num_ind = cluster_sensitivities.size
        sensitivity_dict[i] = np.sum(cluster_sensitivities) / num_ind

    # sort by maximum sensitivity and then return
    # the corresponding n highest labels
    return sorted(sensitivity_dict, key=sensitivity_dict.get, reverse=True)[:n]


def zero_buffer_elements(clusters, num_buffer_elems):
    """
    Description:
        Return clusters with buffer elements set to 0 and
        buffers with roi set to 0
    arguments:
        clusters            [][] : xarray dataarray: the mapping
                                   from the grid cells to state vector number
        num_buffer_elems    int  : number of buffer elements in the state vector
    Returns:                tuple: zeroed buffer area and zeroed roi
    """
    buffer_threshold = int(clusters.max()) - num_buffer_elems
    result = (clusters <= buffer_threshold) * clusters
    buffers = (clusters > buffer_threshold) * clusters
    return result, buffers


def scale_buffer_elements(scale_number, buffer):
    """
    Description:
        Scales buffer elements by the difference in elements between
        the original and scaled cluster elements
    arguments:
        scale_number      int : how much to reduce buffer elements by
        buffer            [][]: ndarray buffer element statevector
    Returns:              [][]: ndarray of scaled buffer elements
    """
    # replace 0's with nan, scale buffer, then replace nan with 0
    buffer = buffer.where(buffer != 0)
    buffer = buffer - scale_number
    buffer = buffer.fillna(0)

    return buffer


def find_cluster_pairs(
    sorted_sensitivities,
    max_dofs,
    desired_elements,
    max_aggregation_level,
    cluster_pairs=None,
):
    """
    Description:
        Recursively generate best-guess at information content distributed
        clustering pairs based on the maximum dofs allowed per cluster and
        the desired number of state vector elements.
    arguments:
        sorted_sensitivities     [] : ndarray of sensitivities sorted in descending order
        max_dofs              float : maximum dofs per state vector element
        desired_elements        int : number of desired elements in state vector
        max_aggregation_level   int : maximum number of elements to aggregate per cluster
        cluster_pairs          dict : clustering pairs
    Returns:                   dict : information content informed clustering pairs
    """
    # Handle initial call
    if cluster_pairs is None:
        cluster_pairs = {}

    # the number of elements that would be needed to create a background of 4x5 degree elements
    background_elements_needed = np.ceil(
        len(sorted_sensitivities) / max_aggregation_level
    )

    # determine the number of elements that should be aggregated
    # based on the number of elements left to be distributed
    # and the maximum dofs per cluster
    # we save enough elements to create a background of 4x5
    # degree state vector elements
    elements_left = desired_elements - sum(cluster_pairs.values())
    if elements_left == 1:
        # aggregate remaining grid cells into a single element
        subset = np.arange(len(sorted_sensitivities))
    elif elements_left == background_elements_needed:
        # start aggregating remainder of elements into clusters
        # of size max_aggregation_level
        subset = np.arange(0, max_aggregation_level)
    elif elements_left < background_elements_needed:
        # allow aggregation of greater than max_aggregation_level
        # This may occur due to rounding
        excess = len(sorted_sensitivities) % max_aggregation_level
        subset = np.arange(max_aggregation_level + excess)
    else:
        # aggregate the number of cells that have a cumulative sum
        # below the dofs threshold
        # assert a minimum of 1 cell per cluster
        cumsum_sensitivities = np.cumsum(sorted_sensitivities)
        subset = np.where(cumsum_sensitivities < max_dofs)[0]
        subset = subset if len(subset) > 0 else [0]
        if len(subset) > max_aggregation_level:
            subset = np.arange(0, max_aggregation_level)

    # handle cases where too few clusters would be created by
    # adding additional native resolution clusters
    if (elements_left - 1) > (len(sorted_sensitivities) - len(subset)):
        subset = [0]

    # update dictionary with the new cluster pairing
    native_cells = len(subset)
    if native_cells in cluster_pairs.keys():
        cluster_pairs[native_cells] = cluster_pairs[native_cells] + 1
    else:
        cluster_pairs[native_cells] = 1
    # delete the sensitivities that have been assigned a pairing
    sorted_sensitivities = np.delete(sorted_sensitivities, subset)
    # recursively find the next pairing
    if len(sorted_sensitivities) == 0:
        return cluster_pairs
    else:
        return find_cluster_pairs(
            sorted_sensitivities,
            max_dofs,
            desired_elements,
            max_aggregation_level,
            cluster_pairs,
        )


def get_max_aggregation_level(config, sensitivities, desired_element_num):
    """
    Description:
        Returns the maximum aggregation level based on the number of desired
        elements and the resolution. By default, if there are enough elements
        we default to using a max aggregation level corresponding to a 8x10
        grid cell (for global).
    arguments:
        config             {dict} : imi config file
        sensitivities    [double] : list of avging kernel senstivities
        desired_element_num   int : desired number of state vector elements
    Returns:                  int : max gridcells per cluster
    """
    if config["Res"] == "2.0x2.5":
        max_aggregation_level = 16 # Setting background to 8x10 for global
    elif config["Res"] == "0.5x0.625":
        max_aggregation_level = 64
    elif config["Res"] == "0.25x0.3125":
        max_aggregation_level = 256

    background_elements_needed = np.ceil(len(sensitivities) / max_aggregation_level)
    if background_elements_needed > desired_element_num:
        print(
            "Warning: too few clusters to create a background of 4x5 degree state vector elements."
            + " More state vector elements recommended. Increasing aggregation level threshold."
        )
        # if there are too few clusters then we set the max aggregation level
        # to either total_native_elements/8 or total_native_elements
        denominator = 8 if desired_element_num > 8 else 1
        max_aggregation_level = np.ceil(len(sensitivities) / denominator)
        print(
            f"Max aggregation level set to: {max_aggregation_level} elements in a cluster"
        )
    print(
        f"Max aggregation level set to: {max_aggregation_level} elements in a cluster"
    )
    return max_aggregation_level


def generate_cluster_pairs(config, sensitivities):
    """
    Description:
        Generate information content informed clustering pairs
    arguments:
        config           {dict} : imi config file
        sensitivities        [] : averaging kernel sensitivities
    Returns:          [(tuple)] : information content informed clustering pairs
    """
    desired_element_num = config["NumberOfElements"] - config["nBufferClusters"]
    # Error handling
    if desired_element_num < 0:
        raise Exception(
            "Error in clustering algorithm: too few clusters requested."
            + f"requested {desired_element_num} clusters."
            + "Remember to take into account the number of buffer elements."
        )
    if desired_element_num > len(sensitivities):
        raise Exception(
            "Error in clustering algorithm: too many clusters requested."
            + f" {len(sensitivities)} native resolution elements and "
            + f"requested {desired_element_num} elements."
            + "Remember to take into account the number of buffer elements."
        )
    # sort sensitivities in descending order
    sensitivities = np.sort(sensitivities)[::-1]

    # maximum number of native elements per cluster
    max_aggregation_level = get_max_aggregation_level(
        config, sensitivities, desired_element_num
    )

    # temporarily set the upper bound limit to prevent recursion error
    # for large domains. python has a default recursion limit of 1000
    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(len(sensitivities))

    # determine dofs threshold for each cluster and create cluster pairings
    target_dofs_per_cluster = sum(sensitivities) / desired_element_num # alternatively hardcode
    print(f"Sum of sensitivities is {sum(sensitivities)}.")
    print(f"Target DOFS per cluster is {target_dofs_per_cluster}.") # make higher threshold to get more aggregated elements
    cluster_pairs = find_cluster_pairs(
        sensitivities,
        target_dofs_per_cluster,
        desired_element_num,
        max_aggregation_level,
    )
    sys.setrecursionlimit(limit)

    # put cluster pairs into format expected by clustering algorithm
    cluster_pairs = list(cluster_pairs.items())
    print(f"Generated cluster pairings: {cluster_pairs}")
    return sorted(cluster_pairs, key=lambda x: x[0])


def force_native_res_pixels(config, clusters, sensitivities):
    """
    Description:
        Forces native resolution for specified coordinates in config file
        by setting the corresponding sensitivity to 1.0
    arguments:
        config           {dict}   : imi config file
        clusters         [][]     : ndarray statevector clusters
        sensitivities    [double] : list of avging kernel senstivities
        cluster_pairs    [(tuple)]: cluster pairings
    Returns:             [double] : updated sensitivities
    """
    coords = get_point_source_coordinates(config)

    if len(coords) == 0:
        # No forced pixels inputted
        print(
            f"No forced native pixels specified or in {config['PointSourceDatasets']} dataset."
        )
        return sensitivities

    if config["Res"] == "0.25x0.3125":
        lat_step = 0.25
        lon_step = 0.3125
    elif config["Res"] == "0.5x0.625":
        lat_step = 0.5
        lon_step = 0.625
    elif config["Res"] == "2.0x2.5":
        lat_step = 2.0
        lon_step = 2.5
    elif config["Res"] == "4.0x5.0":
        lat_step = 4.0
        lon_step = 5.0

    for lat, lon in coords:
        lon = np.floor(lon / lon_step) * lon_step
        lat = np.floor(lat / lat_step) * lat_step

    # Remove any duplicate coordinates within the same gridcell.
    coords = sorted(set(map(tuple, coords)), reverse=True)
    coords = [list(coordinate) for coordinate in coords]

    if len(coords) > config["NumberOfElements"]:
        coords = coords[0 : config["NumberOfElements"] - 1]

    for lat, lon in coords:
        binned_lon = np.floor(lon / lon_step) * lon_step
        binned_lat = np.floor(lat / lat_step) * lat_step

        try:
            cluster_index = int(
                clusters.sel(lat=binned_lat, lon=binned_lon).values.flatten()[0]
            )
            # assign higher than 1 to ensure first assignment
            sensitivities[cluster_index - 1] = 1.1
        except:
            print(
                f"Warning: not forcing pixel at (lat, lon) = ({lat}, {lon})"
                + " because it is not in the specified region of interest."
            )
    return sensitivities


def update_sv_clusters(config, flat_sensi, orig_sv, cluster_pairs):
    """
    Description:
        Reduce a given statevector based on the averaging kernel sensitivities and inputted
        aggregation scheme
    arguments:
        config                 dict : dict constructed from yaml config file
        flat_sensi        []ndarray : 1D array of sensitivity values for each grid cell
        orig_sv       dataarray[][] : original state vector in native resolution
        cluster_pairs     [[tuple]] : list of tuples representing the desired
                                      aggregation scheme eg. [[1,15], [2,23]] would result in
                                      15 native resolution elements and 23 2-cell elements
    Returns:           [][]datarray : reduced dimension state vector
    """
    # check clustering method
    if config["ClusteringMethod"] == "kmeans":
        mini_batch = False
    elif config["ClusteringMethod"] == "mini-batch-kmeans":
        mini_batch = True
    else:
        raise (
            "Error: Invalid Clustering Method. Valid values are: 'kmeans', 'mini-batch-kmeans'."
        )

    desired_num_labels = config["NumberOfElements"] - config["nBufferClusters"]
    last_ROI_element = int(orig_sv["StateVector"].max() - config["nBufferClusters"])

    # copy original sv for updating at end
    new_sv = orig_sv.copy()

    # set buffer elements to 0, retain buffer labels for rejoining
    _, buffer_labels = zero_buffer_elements(
        orig_sv.copy()["StateVector"], config["nBufferClusters"]
    )
    # sv with no buffer elements
    orig_sv_cp = orig_sv.copy()["StateVector"]
    sv = new_sv["StateVector"].where(orig_sv_cp <= last_ROI_element)

    # match sensitivities with coordinates
    sensi = map_sensitivities_to_sv(flat_sensi, orig_sv, last_ROI_element)

    # initialize labels as 0 everywhere in the ROI
    # labels are NaN outside of ROI
    labels = xr.where(buffer_labels == 0, buffer_labels, np.nan)

    # for each agg_level, cluster the data and assign the n_labels
    # with highest total sensitivity to the new label dataset
    for agg_level, n_labels in cluster_pairs:
        # number of unassigned native resolution elements
        elements_left = np.count_nonzero(labels.values == 0)

        # number of clusters yet to be assigned
        clusters_left = desired_num_labels - int(labels.max())

        if clusters_left < n_labels:
            # if there are fewer clusters left to assign than n_labels
            # then evenly distribute the remaining clusters
            # prevents the algorithm from generating one massive cluster
            out_labels = cluster_data_kmeans(
                sensi["Sensitivities"].where(labels == 0), clusters_left, mini_batch
            )
            n_labels = clusters_left
        # clustering for agg_level 1 is just the state vector
        elif agg_level == 1:
            out_labels = sv.values
        else:
            out_labels = cluster_data_kmeans(
                sensi["Sensitivities"].where(labels == 0),
                int(np.round(elements_left / agg_level)),
                mini_batch,
            )
        n_max_labels = get_highest_labels(out_labels, sensi["Sensitivities"], n_labels)

        # assign the n_max_labels to the labels dataset
        # starting from the highest sensitivity label in the dataset
        label_start = int(labels.max()) + 1
        for max_label in n_max_labels:
            # get indices for label assignment
            if label_start != desired_num_labels:
                ind_max = np.where(out_labels == max_label)
            else:
                # if the desired number of labels is reached,
                # assign all unlabeled elements to the last label
                ind_max = np.where(labels.values == 0)
            labels.values[ind_max] = label_start
            label_start += 1

    # scale buffer elements to correct label range
    cluster_number_diff = last_ROI_element - int(labels.max())
    buffer_labels = scale_buffer_elements(cluster_number_diff, buffer_labels)

    # add buffer labels to new state vector
    new_sv["StateVector"] = (
        new_sv["StateVector"].dims,
        np.where(buffer_labels == 0, labels, buffer_labels),
    )

    # format state vector metadata
    new_sv.StateVector.attrs["units"] = "none"
    new_sv.StateVector.encoding["_FillValue"] = -9999
    new_sv.StateVector.encoding["missing_value"] = -9999

    return new_sv


if __name__ == "__main__":
    inversion_path = sys.argv[1]
    config_path = sys.argv[2]
    state_vector_path = sys.argv[3]
    preview_dir = sys.argv[4]
    tropomi_cache = sys.argv[5]
    kf_index = int(sys.argv[6]) if len(sys.argv) > 6 else None
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    output_file = open(f"{inversion_path}/imi_output.log", "a")
    sys.stdout = output_file
    sys.stderr = output_file

    original_clusters = xr.open_dataset(state_vector_path)
    print("Starting aggregation")
    sensitivity_args = [config, state_vector_path, preview_dir, tropomi_cache, False]

    # dynamically generate sensitivities with only a
    # subset of the data if kf_index is not None
    if kf_index is not None:
        print(f"Dynamically generating clusters for period: {kf_index}.")
        sensitivity_args.append(kf_index)

    sensitivities = estimate_averaging_kernel(*sensitivity_args)

    # force point sources to be high resolution by updating sensitivities
    sensitivities = force_native_res_pixels(
        config, original_clusters["StateVector"], sensitivities
    )
    cluster_pairs = generate_cluster_pairs(config, sensitivities)

    print(
        "Creating new clusters based on cluster pairings. "
        + "Run time needed will vary with state vector size."
        + "\nUsing ClusteringMethod: 'mini-batch-kmeans' will "
        + "perform faster, but may reduce cluster accuracy."
    )

    new_sv = update_sv_clusters(config, sensitivities, original_clusters, cluster_pairs)
    original_clusters.close()

    # replace original statevector file
    print(f"Saving file {state_vector_path}")
    new_sv.to_netcdf(
        state_vector_path,
        encoding={v: {"zlib": True, "complevel": 9} for v in new_sv.data_vars},
    )
