#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SBATCH -N 1
# SBATCH -n 1

import xarray as xr
import numpy as np
import pandas as pd
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
    result = clusters.copy()
    c_array = result.values
    c_idx = np.where(c_array > 0)
    c_val = c_array[c_idx]
    row_idx = [r for _, r, _ in sorted(zip(c_val, c_idx[0], c_idx[1]))]
    col_idx = [c for _, _, c in sorted(zip(c_val, c_idx[0], c_idx[1]))]
    idx = (row_idx, col_idx)

    d_idx = np.where(c_array == 0)

    c_array[c_idx] = data
    c_array[d_idx] = default_value
    result.values = c_array

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


############################################
### REDUCED DIMENSION JACOBIAN FUNCTIONS ###
############################################
def kmeans_clustering(label_idx, clusters, n_cluster_size):
    #   n_cells, n_cluster_size=None):

    # get indices of native resolution clusters
    # u, ind, counts = np.unique(clusters, return_index=True, return_counts=True)
    # unique_mask = [c == 1 for c in counts]
    # label_idx = ind[unique_mask]

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

    # Print out some information
    label_stats = np.unique(labels_new.labels_, return_counts=True)
    print("Number of clusters: %d" % len(label_stats[0]))
    print("Cluster size: %d" % n_cluster_size)
    print("Maximum number of grid boxes in a cluster: %d" % max(label_stats[1]))
    print("Average number of grid boxes in a cluster: %.2f" % np.mean(label_stats[1]))
    print("...")

    # Save the information
    labels = labels.assign(new_labels=labels_new.labels_ + 1)  # Pythonic indexing
    labels[["labels", "new_labels"]] = labels[["labels", "new_labels"]].astype(int)

    return labels


def aggregate_cells(clusters, orig_state_vector, sensitivities, n_cells, n_cluster_size=None):
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

    # Print information about state vector
    ux = np.unique(new_sv)
    print("Number of state vector elements: %d" % len(ux))

    return new_sv


def generate_cluster_pairs(clusters, num_buffer_cells, cluster_pairs):
    """
    Description:
        Generate cluster pairs expected by aggregation algorithm
        and validate inputted values
    arguments:
        clusters         [][]     : ndarray of statevector clusters
        num_buffer_cells int      : num buffer elements in inversion domain
        cluster_pairs    [(tuple)]: cluster pairings
    Returns:             [(tuple)]: updated cluster pairings
    """
    native_num_clusters = int(clusters.max()) - num_buffer_cells
    new_cluster_pairs = []
    total_native_cells_requested = 0

    for cells_per_cluster, num_clusters in cluster_pairs:
        native_cells = num_clusters * cells_per_cluster
        total_native_cells_requested += native_cells
        new_cluster_pairs.append((cells_per_cluster, native_cells))

    remainder = native_num_clusters - total_native_cells_requested
    if remainder < 0:
        raise Exception(
            f"Error in cluster pairs: too many pixels requested."
            + f" {native_num_clusters} native resolution pixels and "
            + f"requested cluster pairings use {total_native_cells_requested} pixels."
        )
    elif remainder != 0:
        print(
            "Warning: Cluster pairings do not use all native pixels."
            + f"Adding additional cluster pairing: {[remainder, 1]}."
        )
        new_cluster_pairs = new_cluster_pairs.append((remainder, remainder))

    return new_cluster_pairs


def force_native_res_pixels(config, clusters, sensitivities, cluster_pairs):
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
    coords = config["ForcedNativeResolutionPixels"]

    # Error Handling
    num_native_pixels = [pair[1] for pair in cluster_pairs if pair[0] == 1]
    if len(coords) > num_native_pixels[0]:
        message = "Error: Not enough native resolution pixels for forced coordinates."
        +f"{len(coords)} forced coordinates, but only {num_native_pixels[0]} native pixels."
        raise (Exception(message))

    if config["Res"] == "0.25x0.3125":
        lat_step = 0.25
        lon_step = 0.3125
    elif config["Res"] == "0.5x0.625":
        lat_step = 0.5
        lon_step = 0.625

    for lat, lon in coords:
        binned_lon = np.floor(lon / lon_step) * lon_step
        binned_lat = np.floor(lat / lat_step) * lat_step
        cluster_index = (
            int(clusters.sel(lat=binned_lat, lon=binned_lon).values.flatten()[0])
        )
        # TODO understand why the indexing is backwards
        sensitivities[-cluster_index] = 1.0
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
    cluster_pairs = config["ClusteringPairs"]
    cluster_pairs = generate_cluster_pairs(
        original_clusters["StateVector"], config["nBufferClusters"], cluster_pairs
    )
    print("Starting aggregation")
    tic = time.perf_counter()
    sensitivities = estimate_averaging_kernel(
        config, state_vector_path, preview_dir, tropomi_cache
    )
    toc = time.perf_counter()
    agg_start = time.perf_counter()
    print(f"estimate_averaging_kernel time: {toc-tic}")

    sensitivities = force_native_res_pixels(
        config, original_clusters["StateVector"], sensitivities, cluster_pairs
    )
    new_sv = update_sv_clusters(
        original_clusters, sensitivities, cluster_pairs, config["nBufferClusters"]
    )
    original_clusters.close()
    agg_end = time.perf_counter()
    print(f"update_sv_cluster time: {agg_end-agg_start}")

    # replace original statevector file
    print(f"Saving file {state_vector_path}")
    new_sv.to_netcdf(state_vector_path)
