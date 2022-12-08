import xarray as xr
import numpy as np
import pandas as pd
import copy
import math

# clustering
from sklearn.cluster import KMeans

# Import information for plotting in a consistent fashion
import sys

sys.path.append("./python/")

##########################
### PLOTTING FUNCTIONS ###
##########################


def match_data_to_clusters(data, clusters, default_value=0):
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

def remove_buffer_elements(state_vector, lats, lons):
    """
    Description: 
        Remove buffer elements from state vector so that only region of interest elements 
        remain.
    arguments:
        state_vector    [][] xarray dataset: the mapping from the grid cells to state
                            vector number
    Returns:            xarray dataset containing state vector elements without buffer elements
    """
    # TODO: handle more than just rectangular regions of interest
    return state_vector.where((state_vector.lat <= max(lats) 
        & (state_vector.lat >= min(lats)) 
        & (state_vector.lon >= min(lons)) 
        & (state_vector.lon <= max(lons)))).dropna("lat", how="all").dropna("lon", how="all")

def update_sv_clusters(orig_xr_clusters, sensitivities, clustering_options, method="kmeans"):
    """
    Description: 
        Reduce a given statevector based on the averaging kernel sensitivities and inputted 
        aggregation scheme
    arguments:
        orig_xr_clusters    [][] xarray dataset: the mapping from the grid cells to state
                            vector number
        sensitivities       [int]: the native resolution diagonal of the averaging
                            kernel estimate
        clustering_options  [(tuple)]: list of tuples representing the desired
                                aggregation scheme eg. [(1,15), (2,46)] would result in
                                15 native resolution elements and 23 2-cell elements
    """
    if method == "kmeans":
        new_sv = orig_xr_clusters.copy()
        sv_element_list = np.arange(1, orig_xr_clusters["StateVector"].max() + 1, 1)
        aggregation_sizes = [item[0] for item in clustering_options]
        num_aggregation_elements = [item[1] for item in clustering_options]

        new_cluster_labels = aggregate_cells(
            orig_xr_clusters["StateVector"],
            sv_element_list,
            sensitivities,
            num_aggregation_elements,
            aggregation_sizes,
        )
        new_sv["StateVector"] = match_data_to_clusters(
            new_cluster_labels, orig_xr_clusters["StateVector"]
        )
    else:
        raise("Error: Invalid Clustering Method Inputted. Valid values are: \'kmeans\'.")
    return new_sv


# def plot_multiscale_grid(self, clusters, **kw):
#     # Get KW
#     title = kw.pop('title', '')
#     fig_kwargs = kw.pop('fig_kwargs', {})
#     title_kwargs = kw.pop('title_kwargs', {})
#     map_kwargs = kw.pop('map_kwargs', {})
#     kw['colors'] = kw.pop('colors', 'black')
#     kw['linewidths'] = kw.pop('linewidths', 1)

#     # Plot
#     nstate = len(np.unique(self.state_vector)[1:])
#     data = self.match_data_to_clusters(self.state_vector,
#                                        clusters, default_value=0)
#     data_zoomed = zoom(data.values, 50, order=0, mode='nearest')
#     fig, ax = fp.get_figax(maps=True, lats=data.lat, lons=data.lon,
#                            **fig_kwargs)
#     ax.contour(data_zoomed, levels=np.arange(0, nstate, 1),
#                extent=[data.lon.min(), data.lon.max(),
#                        data.lat.min(), data.lat.max()],
#                **kw)
#     ax = fp.add_title(ax, title, **title_kwargs)
#     ax = fp.format_map(ax, data.lat, data.lon, **map_kwargs)

#     return fig, ax

############################################
### REDUCED DIMENSION JACOBIAN FUNCTIONS ###
############################################
def kmeans_clustering(label_idx, clusters, n_cluster_size):
    #   n_cells, n_cluster_size=None):

    # * dofs isn't used???
    # * is label_idx the indices of clusters with native res?
    # get indices of native resolution clusters
    # u, ind, counts = np.unique(clusters, return_index=True, return_counts=True)
    # unique_mask = [c == 1 for c in counts]
    # label_idx = ind[unique_mask]

    # * clusters is 2d and labels is 1d?
    # Initialize the labels and set the labels of interest to
    # non zero values
    labels = np.zeros(int(clusters.max()))  # Get zeros of nstate dimension
    labels[label_idx] = label_idx + 1  # Fix for pythonic indexing

    # * not sure if I need to use this func because my original sv is already labelled?
    # * so labels might just be sv.values?
    # Move the labels into a 2D cluster format to get lat/lon information
    # so we can cluster proximate grid cells
    labels = match_data_to_clusters(labels, clusters)

    # * are labels of interest the ones with high dofs? Not sure where dofs comes in here
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
    # * As is, this looks like it would fit based on spatial proximity?
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


def aggregate_cells(clusters, orig_state_vector, dofs, n_cells, n_cluster_size=None):
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
        clusters            [][] ndarray: the mapping from the grid cells to state
                            vector number
        orig_state_vector   [int]: the previous state vector #* just a list of 1:num_vs_element?
        dofs                [int]: the native resolution diagonal of the averaging
                            kernel estimate #* reasonable to use random 0-1 values? 2d like the clusters?
        rf                  int: the native resolution regularization factor #* not actually a parameter?
        n_cells             [int]: the number of native resolution grid
                            cells to be used in the aggregation
                            scheme (integer or list of integers);
                            defaults to [100, 200] #* what are subsequent indices used for?
        n_cluster_size      [int]: the number of native resolution grid
                            cells to aggregate together at each level
                            of aggregation; defaults to [1, 2] #* what is this used for?
    Example:
        Passing n_cells=[100, 200], n_cluster_size=[1, 2] will
        generate a state vector where the 100 grid cells with highest
        information content maintain the native resolution (i.e.
        cluster size = 1) and the 200 grid cells with next highest
        information content will be aggregated into clusters with size
        2. The final state vector will have dimension 200.
    """

    print("... Generating multiscale grid ...")
    rf = 1.0

    new_sv = copy.deepcopy(orig_state_vector)  # (New state vector)
    dofs = copy.deepcopy(dofs)
    orig_elements, orig_counts = np.unique(orig_state_vector, return_counts=True)

    # For any grid cells previously added to the state vector,
    # set significance to 0
    # if np.any(orig_counts > 1): # If we aggregated together any grid cells
    #     print('Ignoring previously optimized grid cells.')
    #     dofs[orig_counts == 1] = 0 # Then assume that any grid cells with
    # only one element are previously optimized
    ## Note to Lucas: I think this might be wrong...?

    # Get the indices associated with the most significant
    # grid boxes
    sig_idx = dofs.argsort()[::-1]

    # Initialize the new state vector
    new_sv = np.zeros(len(orig_state_vector))

    # Iterate through n_cells
    n_cells = np.append(0, n_cells)
    nidx = np.cumsum(n_cells)
    # * I dont understand why this is in a loop (1,100), (2,200) or maybe what ncells is?
    for i, n in enumerate(n_cells[1:]):
        # Get cluster size
        if n_cluster_size is None:
            cluster_size = i + 2
        else:
            cluster_size = n_cluster_size[i]

        # Select the indices corresponding to the largest availble
        # DOFS elements
        sub_sig_idx = sig_idx[nidx[i] : nidx[i + 1]]

        # Get the new labels
        # * params don't match function?
        # new_labels = kmeans_clustering(sub_sig_idx, clusters, cluster_size, n_cells)
        new_labels = kmeans_clustering(sub_sig_idx, clusters, cluster_size)

        # Adjust the labels (which are from 1 to m) to match the range
        # of values already in the state vector
        new_sv[new_labels["labels"]] = new_labels["new_labels"] + new_sv.max()

    # Print information about state vector
    ux, ux_cnt = np.unique(new_sv, return_counts=True)
    summ, summ_cnt = np.unique(ux_cnt, return_counts=True)

    # Adjust the rf
    new_elements, counts = np.unique(new_sv, return_counts=True)
    new_rf = rf * len(new_elements) / len(orig_elements)

    print("Number of state vector elements: %d" % len(ux))
    print(summ, summ_cnt)
    print("... Complete ...\n")

    return new_sv


def calculate_k_ms(forward_model, state_vector):
    k_ms = pd.DataFrame(forward_model)
    k_ms = k_ms.groupby(state_vector, axis=1).sum()
    k = np.array(k_ms)
    return k


def calculate_prior_ms(xa_abs, sa_vec, state_vector):
    xa_abs_ms = pd.DataFrame(xa_abs).groupby(state_vector).sum()

    xa_abs = np.array(xa_abs_ms).reshape(
        -1,
    )
    xa = np.ones(len(xa_abs)).reshape(
        -1,
    )
    sa_vec = 0.25 * np.ones(len(xa_abs)).reshape(
        -1,
    )

    return xa, xa_abs, sa_vec
