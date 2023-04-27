#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
import numpy as np
import yaml
import copy
import sys
import time
from inversion_scripts.imi_preview import estimate_averaging_kernel

# clustering
from sklearn.cluster import KMeans


def match_data_to_clusters(data, clusters, default_value=0):
    """
    Description:
        match 1D list of elements to 2D statevector
    arguments:
        data                [int]: data to infill into 2D statevector
        clusters            [][] : xarray dataarray: the mapping from
                                   the grid cells to state vector number
    Returns:                [][] : filled with new values
    """
    reference = clusters.copy()  # reference copy for determining mapped elements
    result = clusters.copy()
    # map sensitivities onto corresponding xarray DataArray
    for i in range(1, int(reference.max()) + 1):
        mask = reference == i
        result = xr.where(mask, data[i - 1], result)
    return result


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


def update_sv_clusters(
    orig_xr_clusters,
    sensitivities,
    clustering_options,
    num_buffer_elems,
    method="kmeans",
):
    """
    Description:
        Reduce a given statevector based on the averaging kernel sensitivities and inputted
        aggregation scheme
    arguments:
        orig_xr_clusters    [][]: xarray dataset: the mapping from the grid cells to state
                            vector number
        sensitivities       [float32]: the native resolution diagonal of the averaging
                            kernel estimate. Only sensitivites for the ROI
        clustering_options  [[tuple]]: list of tuples representing the desired
                                aggregation scheme eg. [[1,15], [2,46]] would result in
                                15 native resolution elements and 23 2-cell elements
        num_buffer_elems    int: number of buffer elements in statevector
    """
    if method == "kmeans":
        clusters_copy = orig_xr_clusters.copy()
        # set buffer elements to 0, retain buffer labels for rejoining
        roi_labels, buffer_labels = zero_buffer_elements(
            clusters_copy["StateVector"], num_buffer_elems
        )
        clusters_copy["StateVector"] = roi_labels

        new_sv = clusters_copy.copy()
        sv_element_list = np.arange(1, clusters_copy["StateVector"].max() + 1, 1)
        aggregation_sizes = [item[0] for item in clustering_options]
        num_aggregation_elements = [item[1] for item in clustering_options]

        new_cluster_labels = aggregate_cells(
            clusters_copy["StateVector"],
            sv_element_list,
            sensitivities,
            num_aggregation_elements,
            aggregation_sizes,
        )

        new_sv_labels = match_data_to_clusters(
            new_cluster_labels, clusters_copy["StateVector"]
        )

        # scale buffer elements to correct label range
        cluster_number_diff = int(clusters_copy["StateVector"].max()) - int(
            new_sv_labels.max()
        )
        buffer_labels = scale_buffer_elements(cluster_number_diff, buffer_labels)

        # add buffer labels to new state vector
        new_sv["StateVector"] = (
            new_sv["StateVector"].dims,
            np.where(new_sv_labels == 0, buffer_labels, new_sv_labels),
        )
    else:
        raise ("Error: Invalid Clustering Method Inputted. Valid values are: 'kmeans'.")

    # format state vector metadata
    new_sv.StateVector.attrs["units"] = "none"
    new_sv.StateVector.encoding["_FillValue"] = -9999
    new_sv.StateVector.encoding["missing_value"] = -9999

    return new_sv


def kmeans_clustering(label_idx, clusters, n_cluster_size):
    """
    Description:
        cluster remaining unlabeled pixels into cluster of approximately n_cluster_size
    arguments:
        label_idx         [] : indices of gridcells with sorted by sensitivity magnitude
        clusters        [][] : ndarray mapping from the grid cells to state vector number
        n_cluster_size   int : desired cluster size
    Returns:              [] : new labels for gridcells
    """

    # Initialize the labels and set the labels of interest to
    # non zero values
    labels = np.zeros(int(clusters.max()))  # Get zeros of nstate dimension
    labels[label_idx] = label_idx + 1  # Fix for pythonic indexing

    # Move the labels into a 2D cluster format to get lat/lon information
    # so we can cluster proximate grid cells
    labels = match_data_to_clusters(labels, clusters)

    # Move to a dataframe and only choose the labels of interest (label_idx)
    labels = labels.to_dataframe("labels").reset_index()
    labels = labels[labels["labels"] > 0]
    labels["labels"] -= 1  # Undo pythonic indexing fix

    # Now do kmeans clustering
    ## The numer of clusters fed to the algorithm is given by the number
    ## of native resolution grid cells available (len(label_idx)) divided
    ## by the desired cluster size
    n_clusters = int(len(label_idx) / n_cluster_size)

    ## These two lines do the Kmeans clustering and then apply the
    ## labels to the grid cells.
    ## I specify the random_state to avoid variations between test rounds
    ## since Kmeans has a degree of randomness
    labels_new = KMeans(n_clusters=n_clusters, random_state=0)
    labels_new = labels_new.fit(labels[["lat", "lon"]])

    # Save the information
    labels = labels.assign(new_labels=labels_new.labels_ + 1)  # Pythonic indexing
    labels[["labels", "new_labels"]] = labels[["labels", "new_labels"]].astype(int)

    return labels


def aggregate_cells(
    clusters, orig_state_vector, sensitivities, n_cells, n_cluster_size=None
):
    """
    This function generates a multi-scale Jacobian on the basis
    of the information content, as given by the diagonal elements
    of the averaging kernel. It maintains the native resolution
    for the grid cells with the highest information content while
    aggregating together grid cells with lower information content
    using K means clustering. The user must provide the number of
    desired number of grid cells and the number of grid cells per
    cluster for each aggregation level. It accepts the following
    arguments:
        clusters           [][]     : ndarray mapping from the grid cells to state vector number
        orig_state_vector  [int]    : the previous state vector
        sensitivities      [float32]: the native resolution diagonal of the averaging
                                      kernel estimate
        n_cells            [int]    : the number of native resolution grid cells to be used
                                      in the aggregation scheme (integer or list of integers);
                                      defaults to [100, 200]
        n_cluster_size     [int]    : the number of native resolution grid cells to aggregate
                                      together at each level of aggregation; defaults to [1, 2]
    Example:
        Passing n_cells=[100, 200], n_cluster_size=[1, 2] will
        generate a state vector where the 100 grid cells with highest
        information content maintain the native resolution (i.e.
        cluster size = 1) and the 200 grid cells with next highest
        information content will be aggregated into clusters with size
        2. The final state vector will have dimension 200.
    """

    new_sv = copy.deepcopy(orig_state_vector)  # (New state vector)
    sensitivities = copy.deepcopy(sensitivities)

    # Get the indices associated with the most significant
    # grid boxes
    sig_idx = sensitivities.argsort()[::-1]

    # Initialize the new state vector
    new_sv = np.zeros(len(orig_state_vector))

    # Iterate through n_cells
    n_cells = np.append(0, n_cells)
    nidx = np.cumsum(n_cells)

    for i, _ in enumerate(n_cells[1:]):
        # Get cluster size
        if n_cluster_size is None:
            raise ("Error: n_cluster_size is None")
        else:
            cluster_size = n_cluster_size[i]

        # Select the indices corresponding to the largest availble
        # DOFS elements
        sub_sig_idx = sig_idx[nidx[i] : nidx[i + 1]]

        # Get the new labels
        new_labels = kmeans_clustering(sub_sig_idx, clusters, cluster_size)

        # Adjust the labels (which are from 1 to m) to match the range
        # of values already in the state vector
        new_sv[new_labels["labels"]] = new_labels["new_labels"] + new_sv.max()

    return new_sv


def find_cluster_pairs(
    sorted_sensitivities,
    max_dofs,
    desired_elements,
    max_aggregation_level,
    cluster_pairs={},
):
    """
    Description:
        Recursively generate optimal clustering pairs based on the
        maximum dofs allowed per cluster and the desired number of
        state vector elements.
    arguments:
        sorted_sensitivities     [] : ndarray of sensitivities sorted in descending order
        max_dofs              float : maximum dofs per state vector element
        desired_elements        int : number of desired elements in state vector
        max_aggregation_level   int : maximum number of elements to aggregate per cluster
        cluster_pairs           dict: optimally distributed clustering pairs
    Returns:                    dict: optimal cluster pairings
    """
    # the number of elements that would be needed to create a background of 4x5 degree elements
    background_elements_needed = np.ceil(len(sorted_sensitivities) / max_aggregation_level)
    
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
        excess = len(sorted_sensitivities)%max_aggregation_level
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
        we default to using a max aggregation level of corresponding to a 4x5
        grid cell.
    arguments:
        config             {dict} : imi config file
        sensitivities    [double] : list of avging kernel senstivities
        desired_element_num   int : desired number of state vector elements
    Returns:                  int : max gridcells per cluster
    """
    if config["Res"] == "0.25x0.3125":
        max_aggregation_level = 256
    elif config["Res"] == "0.5x0.625":
        max_aggregation_level = 64
    
    background_elements_needed = len(sensitivities) / max_aggregation_level
    if background_elements_needed > desired_element_num:
        print(
            "Warning: too few clusters to create a background of 4x5 degree state vector elements."
            + " More state vector elements recommended. Increasing aggregation level threshold."
        )
        # if there are too few clusters then we set the max aggregation level
        # to either total_native_elements/8 or total_native_elements
        denominator = 8 if desired_element_num > 8 else 1
        max_aggregation_level = len(sensitivities) / denominator
        print(
            f"Max aggregation level set to: {max_aggregation_level} elements in a cluster"
        )
    print(f"Max aggregation level set to: {max_aggregation_level} elements in a cluster")
    return max_aggregation_level


def generate_cluster_pairs(config, sensitivities):
    """
    Description:
        Generate optimal clustering pairs
    arguments:
        config           {dict} : imi config file
        sensitivities        [] : averaging kernel sensitivities
    Returns:          [(tuple)] : optimal cluster pairings
    """
    desired_element_num = config["NumberOfElements"] - config["nBufferClusters"]
    # Error handling
    if desired_element_num < 0:
        raise Exception(
            f"Error in clustering algorithm: too few clusters requested."
            + f"requested {desired_element_num} clusters."
            + "Remember to take into account the number of buffer elements."
        )
    if desired_element_num > len(sensitivities):
        raise Exception(
            f"Error in clustering algorithm: too many clusters requested."
            + f" {len(sensitivities)} native resolution elements and "
            + f"requested {desired_element_num} elements."
            + "Remember to take into account the number of buffer elements."
        )
    # sort sensitivities in ascending order
    sensitivities = np.sort(sensitivities)[::-1]

    # maximum number of native elements per cluster
    max_aggregation_level = get_max_aggregation_level(
        config, sensitivities, desired_element_num
    )

    # determine dofs threshold for each cluster and create optimal pairings
    target_dofs_per_cluster = sum(sensitivities) / desired_element_num
    cluster_pairs = find_cluster_pairs(
        sensitivities,
        target_dofs_per_cluster,
        desired_element_num,
        max_aggregation_level,
    )

    # put cluster pairs into format expected by clustering algorithm
    cluster_pairs = list(cluster_pairs.items())
    print(f"Generated cluster pairings: {cluster_pairs}")
    new_cluster_pairs = []

    for cells_per_cluster, num_clusters in cluster_pairs:
        native_cells = num_clusters * cells_per_cluster
        new_cluster_pairs.append((cells_per_cluster, native_cells))

    # sort cluster pairs in ascending order
    return sorted(new_cluster_pairs, key=lambda x: x[0])


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
    coords = config["ForcedNativeResolutionElements"]

    if coords is None:
        # No forced pixels inputted
        return sensitivities

    if config["Res"] == "0.25x0.3125":
        lat_step = 0.25
        lon_step = 0.3125
    elif config["Res"] == "0.5x0.625":
        lat_step = 0.5
        lon_step = 0.625

    for lat, lon in coords:
        binned_lon = np.floor(lon / lon_step) * lon_step
        binned_lat = np.floor(lat / lat_step) * lat_step
        cluster_index = int(
            clusters.sel(lat=binned_lat, lon=binned_lon).values.flatten()[0]
        )
        sensitivities[cluster_index - 1] = 1.0
    return sensitivities


if __name__ == "__main__":
    inversion_path = sys.argv[1]
    config_path = sys.argv[2]
    state_vector_path = sys.argv[3]
    preview_dir = sys.argv[4]
    tropomi_cache = sys.argv[5]
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    output_file = open(f"{inversion_path}/imi_output.log", "a")
    sys.stdout = output_file
    sys.stderr = output_file

    original_clusters = xr.open_dataset(state_vector_path)
    print("Starting aggregation")
    sensitivities = estimate_averaging_kernel(
        config, state_vector_path, preview_dir, tropomi_cache
    )
    if "ForcedNativeResolutionElements" in config.keys():
        sensitivities = force_native_res_pixels(
            config, original_clusters["StateVector"], sensitivities
        )
    cluster_pairs = generate_cluster_pairs(config, sensitivities)

    print("Creating new clusters based on cluster pairings. "+"Run time needed will vary with state vector size.")
    
    new_sv = update_sv_clusters(
        original_clusters, sensitivities, cluster_pairs, config["nBufferClusters"]
    )
    original_clusters.close()

    # replace original statevector file
    print(f"Saving file {state_vector_path}")
    new_sv.to_netcdf(state_vector_path)
