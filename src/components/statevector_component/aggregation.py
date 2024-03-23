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


def get_highest_labels_threshold(labels, sensitivities, threshold):
    """
    Description:
        Get labels with sensitivity above the threshold. 
        Returned in descending order.
    arguments:
        labels          [][]  ndarray : state vector labels
        sensitivities   [][]dataarray : xarrray sensitivity data
        threshold               float : number between 0 and 1
    Returns: ([], int, []) : tuple with highest sensitivity labels, 
                             number of labels(n), and list of sensitivity 
                             values
    """
    sensitivity_dict = {}
    max_label = int(np.nanmax(labels))
    n = 0

    total_sensis = []
    # calculate avg dofs per gridcell for each cluster
    for i in range(1, max_label + 1):
        indices = np.where(labels == i)
        cluster_sensitivities = sensitivities[indices]
        num_ind = cluster_sensitivities.size
        total_sensi = np.sum(cluster_sensitivities)
        # avg_sensi = total_sensi / num_ind
        sensitivity_dict[i] = total_sensi
        total_sensis.append(np.round(total_sensi, 2))
        if total_sensi >= threshold:
            n += 1

    n_clusters = sorted(sensitivity_dict, key=sensitivity_dict.get, reverse=True)
    inds = np.array(n_clusters) - 1
    n_sensis = np.array(total_sensis)[inds]
    # sort by maximum sensitivity and then return
    # the corresponding n highest labels
    return (n_clusters[:n], n, n_sensis[:n])


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

def get_max_cluster_size(config, sensitivities, desired_element_num):
    """
    Description:
        Returns the maximum number of elements allowed in a single cluster
        based on the number of desired elements and the resolution. 
        By default, if there are enough elements we default to using a max 
        cluster size of 64 elements.
    arguments:
        config             {dict} : imi config file
        sensitivities    [double] : list of avging kernel senstivities
        desired_element_num   int : desired number of state vector elements
    Returns:                  int : max gridcells per cluster
    """
    max_cluster_size = config["MaxClusterSize"] if "MaxClusterSize" in config.keys() else 64
    
    background_elements_needed = np.ceil(len(sensitivities) / max_cluster_size)
    if background_elements_needed > desired_element_num:
        raise (
            "Warning: too few clusters to have background state vector elements\n"
            + f"aggregating {max_cluster_size} native resolution elements.\n"
            + "More state vector elements recommended. Alternatively, raise the\n"
            + "\"MaxClusterSize\" allowed for the region of interest."
        )
    print(
        f"MaxClusterSize set to: {max_cluster_size} elements in a cluster"
    )
    return max_cluster_size


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
    
    # make sure elements are native res by asserting higher sensitivity
    # than clustering threshold
    if "ClusteringThreshold" in config.keys():
        dofs_max = float(config["ClusteringThreshold"]) + 0.1
    else:
        dofs_max = 1.1
        
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
            sensitivities[cluster_index - 1] = dofs_max
        except:
            print(
                f"Warning: not forcing pixel at (lat, lon) = ({lat}, {lon})"
                + " because it is not in the specified region of interest."
            )
    return sensitivities


def update_sv_clusters(config, flat_sensi, orig_sv):
    """
    Description:
        Reduce a given statevector based on the averaging kernel sensitivities and inputted
        aggregation scheme
    arguments:
        config                 dict : dict constructed from yaml config file
        flat_sensi        []ndarray : 1D array of sensitivity values for each grid cell
        orig_sv       dataarray[][] : original state vector in native resolution
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

    # set dofs threshold based on user preferences
    if "ClusteringThreshold" in config.keys():
        dofs_threshold = float(config["ClusteringThreshold"])
    else:
        # default is to use the avg dofs per element
        dofs_threshold = sum(sensitivities) / desired_num_labels
    
    print(f"Target DOFS per cluster (ClusteringThreshold): {dofs_threshold}")

    # max cluster size based on user preferences
    max_cluster_size = (
        config["MaxClusterSize"] if "MaxClusterSize" in config.keys() else 64
    )
    max_cluster_size = get_max_cluster_size(config, sensitivities, desired_num_labels)

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

    # get list of different aggregateion levels
    # where an aggregation level is (roughly) the number of 
    # elements in a cluster
    cluster_pairs = np.arange(1, max_cluster_size + 1)
    
    # boolean for whether to evenly fill the grid with the 
    # remaining cluster elements
    fill_grid = False
    
    print(f"Reducing to {desired_num_labels} elements")

    # total sensitivities of clusters added to the state vector
    cluster_sensis = []
    # for each agg_level, cluster the data and assign the n_labels
    # with highest total sensitivity to the new label dataset
    for agg_level in cluster_pairs:
        # fill grid if we are on the final round of clusters
        if agg_level == max_cluster_size:
            fill_grid = True
            
        # number of unassigned native resolution elements
        elements_left = np.count_nonzero(labels.values == 0)
        
        # check if no more clustering to do, exit the loop
        if elements_left == 0:
            break

        # number of clusters yet to be assigned
        clusters_left = desired_num_labels - int(labels.max())

        # approximate number labels needed to fill the grid
        # with max_cluster_size clusters
        backfill_num = int(elements_left / max_cluster_size)

        if fill_grid or clusters_left <= backfill_num:
            # if there are fewer clusters left to assign than n_labels
            # then evenly distribute the remaining clusters
            # prevents the algorithm from generating one massive cluster
            out_labels = cluster_data_kmeans(
                sensi["Sensitivities"].where(labels == 0), clusters_left, mini_batch
            )
        # clustering for agg_level 1 is just the state vector
        elif agg_level == 1:
            out_labels = sv.values
        else:
            # generate clusters that are approximately agg_level in size
            out_labels = cluster_data_kmeans(
                sensi["Sensitivities"].where(labels == 0),
                int(np.round(elements_left / agg_level)),
                mini_batch,
            )
            
        # assign all remaining clusters if filling the grid
        # by assigning dofs_threshold to artificially low value
        if fill_grid:
            dofs_threshold = -1

        # get the n_highes labels with sensitivities above the 
        # dofs threshold
        n_max_labels, n_highest, _ = get_highest_labels_threshold(
            out_labels, sensi["Sensitivities"], dofs_threshold
        )
        
        # if too many labels to assign, then we need to assign
        # only some of them and leave enough labels to fill
        # the grid with the max_cluster_size clusters
        if (clusters_left - backfill_num) < n_highest and not fill_grid:
            fill_grid = True
            n_ind = clusters_left - backfill_num
            n_max_labels = n_max_labels[:n_ind]

        print(f"assigning {len(n_max_labels)} labels with agg level: {agg_level}")
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
    print(
        "Creating new clusters based on cluster pairings. "
        + "Run time needed will vary with state vector size."
        + "\nUsing ClusteringMethod: 'mini-batch-kmeans' will "
        + "perform faster, but may reduce cluster accuracy."
    )
    # generate multi resolution state vector
    new_sv = update_sv_clusters(config, sensitivities, original_clusters)
    original_clusters.close()

    # replace original statevector file
    print(f"Saving file {state_vector_path}")
    new_sv.to_netcdf(
        state_vector_path,
        encoding={v: {"zlib": True, "complevel": 9} for v in new_sv.data_vars},
    )
