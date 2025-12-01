import os
import numpy as np
import xarray as xr
from scipy.sparse import csr_matrix, coo_matrix

def get_regrid_weights_ESMF(config, inv_directory, ref_config):
    CSgridDir = os.path.join(os.path.expandvars(inv_directory), 'CS_grids')
    
    ref_CS_RES = ref_config['CS_RES']
    ref_STRETCH_FACTOR = ref_config['STRETCH_FACTOR']
    ref_TARGET_LAT = ref_config['TARGET_LAT']
    ref_TARGET_LON = ref_config['TARGET_LON']
    src_shape = (6, ref_CS_RES, ref_CS_RES)
    
    if config['UseGCHP']:
        CS_RES = config['CS_RES']
        regrid_weight_fpath = f"regrid_weights_c{ref_CS_RES}_s{ref_STRETCH_FACTOR:.1f}_{ref_TARGET_LAT:.1f}N_{ref_TARGET_LON:.1f}E_to_c{CS_RES}_conserve.nc"
        sv_shape = (6, CS_RES, CS_RES)
    else:
        sv_fpath = os.path.join(os.path.expandvars(inv_directory), 'StateVector.nc')
        with xr.open_dataset(sv_fpath).squeeze() as sv_ds:
            sv_lat = sv_ds['lat'].values
            sv_lon = sv_ds['lon'].values
        sv_shape = (len(sv_lat), len(sv_lon))
        nlat = len(sv_lat)
        nlon = len(sv_lon)
        regrid_weight_fpath = f"regrid_weights_c{ref_CS_RES}_s{ref_STRETCH_FACTOR:.1f}_{ref_TARGET_LAT:.1f}N_{ref_TARGET_LON:.1f}E_to_latlon_{nlat}x{nlon}_conserve.nc"
    dst_shape = sv_shape
    
    regrid_weight_ds = xr.open_dataset(os.path.join(CSgridDir, regrid_weight_fpath))
    
    return src_shape, dst_shape, regrid_weight_ds

def get_superobs_index_GCC(config, inv_directory, GC_index):
    RunName = config['RunName']
    gc_fpath = os.path.join(os.path.expandvars(inv_directory), 
        f"jacobian_runs/{RunName}_0000/OutputDir/GEOSChem.SpeciesConc.{config['StartDate']}_0000z.nc4")
    with xr.open_dataset(gc_fpath).squeeze() as gc_ds:
        gc_lat = gc_ds['lat'].values
        gc_lon = gc_ds['lon'].values
    gc_shape = (len(gc_lat), len(gc_lon))
    lati, loni = np.unravel_index(GC_index, gc_shape)
    
    sv_fpath = os.path.join(os.path.expandvars(inv_directory), 'StateVector.nc')
    with xr.open_dataset(sv_fpath).squeeze() as sv_ds:
        sv_lat = sv_ds['lat'].values
        sv_lon = sv_ds['lon'].values
    sv_shape = (len(sv_lat), len(sv_lon))
    
    # Check if grids are truly identical
    grids_identical = (
        gc_shape == sv_shape
        and np.allclose(gc_lat, sv_lat)
        and np.allclose(gc_lon, sv_lon)
    )

    if grids_identical:
        # grids identical: index mapping is 1:1
        superobs_index = GC_index
    else:
        # general case: sv_ds may be larger or smaller or just different
        ind_lat = np.abs(sv_lat[:, None] - gc_lat[lati]).argmin(axis=0)
        ind_lon = np.abs(sv_lon[:, None] - gc_lon[loni]).argmin(axis=0)
        superobs_index = np.ravel_multi_index((ind_lat, ind_lon), sv_shape)
    
    return superobs_index
    
def get_regrid_weights_jacobian_row(config, inv_directory, GC_index, ref_config, ref_GC_index):
    # The row of Jacobian matrix K is super observations, 
    # or satellite observatiosn sampled at simulation grid
    
    # GC_index: flattened simulation index with valid super observations
    # note that GC_index for GCC is from the simulation diagnostics, 
    # it is not necessarily the same as the state vector netCDF file 
    # as state vector at GCC could be trimmed;
    # GC_index is consistent between simulation diagnostics and state vector file for GCHP
    
    # The source grid to create regridding weights is from state vector netCDF file.
    # To avoid creating regridding too many regridding weights, 
    # here I just reuse the regridding weights created during 
    # the regridding process of state vector in IMI setup
    
    # preprocess the destination GC_index at state vector grid for GCC
    # the resulting sv_index contains the indicies of a grid cell with valid super observations
    if not config['UseGCHP']:
        superobs_index = get_superobs_index_GCC(config, inv_directory, GC_index)
    else:
        superobs_index = GC_index
    
    src_shape, dst_shape, regrid_weight_ds = get_regrid_weights_ESMF(config, inv_directory, ref_config)
    # The index in the ESMF regridding weights netCDF file is 1-based
    src_ind = regrid_weight_ds['col'].values - 1
    dst_ind = regrid_weight_ds['row'].values - 1
    regrid_weights = regrid_weight_ds['S'].values
    
    n_dst = np.prod(dst_shape)
    n_src = np.prod(src_shape)
    regrid_weights_csr = csr_matrix((regrid_weights, (dst_ind, src_ind)), 
                                    shape=(n_dst, n_src))
    
    src_mask = np.zeros(np.prod(src_shape), dtype=bool)
    src_mask[ref_GC_index] = True
    
    dst_mask = np.zeros(np.prod(dst_shape), dtype=bool)
    dst_mask[superobs_index] = True
    
    # Subset sparse matrix to grid cell with valid obs
    W = regrid_weights_csr[dst_mask, :][:, src_mask] # (n_dst_valid, n_src_ref)
    
    # renormalize weights to sum to 1 over source grid cells for each destination grid cell
    row_sums = np.asarray(W.sum(axis=1)).ravel()
    W.data /= np.repeat(row_sums, np.diff(W.indptr)) 
    
    # --- enforce same order as GC_index and ref_GC_index ---
    dst_positions = np.where(dst_mask)[0]
    src_positions = np.where(src_mask)[0]
    
    # Build mapping from position -> sequential index in mask
    dst_order = np.searchsorted(dst_positions, superobs_index)
    src_order = np.searchsorted(src_positions, ref_GC_index)

    # Reindex W_sub to match the desired order
    W_ordered = W[dst_order, :][:, src_order]

    return W_ordered

def sum_and_sort_along_statevector(val, sv, fill_value=np.nan):
    # This is useful for getting state vector area
    if np.isnan(fill_value):
        mask = ~np.isnan(sv)
    else:
        mask = sv != fill_value
    sv_valid = sv[mask]
    val_valid = val[mask]
    _, inv = np.unique(sv_valid, return_inverse=True)
    sv_val = np.bincount(inv, weights=val_valid)
    return sv_val

def median_and_sort_along_statevector(val, sv, fill_value=np.nan):
    """
    Compute the median of values grouped by state vector ID.

    Parameters
    ----------
    val : array_like
        Values (e.g., per gridcell or pixel) to aggregate.
    sv : array_like
        State vector IDs (same shape as val). May include NaN or fill_value.
    fill_value : scalar, default np.nan
        Fill value marking invalid SV entries. If NaN, NaN entries are ignored.

    Returns
    -------
    sv_unique : np.ndarray
        Sorted unique valid SV IDs.
    val_median : np.ndarray
        Median of `val` per SV ID, same order as sv_unique.
        
    This function is used for getting prior emissions per state vector ID
    when state vector clustering is applied.
    """
    val = np.asarray(val)
    sv = np.asarray(sv)

    # --- Mask invalid entries ---
    if np.isnan(fill_value):
        mask = ~np.isnan(sv)
    else:
        mask = sv != fill_value

    sv_valid = sv[mask].astype(np.int64, copy=False)
    val_valid = val[mask]

    # --- Sort by ID so equal IDs are contiguous ---
    order = np.argsort(sv_valid)
    sv_sorted = sv_valid[order]
    val_sorted = val_valid[order]

    # --- Find boundaries between IDs ---
    sv_unique, idx_start = np.unique(sv_sorted, return_index=True)

    # --- Split and compute median per ID ---
    splits = np.split(val_sorted, idx_start[1:])
    val_median = np.array([np.nanmedian(s) for s in splits], dtype=val_valid.dtype)

    return val_median

def sum_and_sort_along_statevector_matrix(W_csr, dst_ids_raw, src_ids_raw,
                                          sort_dst=True, sort_src=True):
    """
    Aggregate overlap matrix by state vector IDs, optionally along dst and/or src.

    Parameters
    ----------
    W_csr : csr_matrix, shape (n_dst, n_src)
        Overlap-area matrix for valid (masked) state vector elements.
    dst_ids_raw : (n_dst,) array
        Destination state vector IDs (1-based, >0), aligned with W_csr rows.
    src_ids_raw : (n_src,) array
        Source state vector IDs (1-based, >0), aligned with W_csr cols.
    sort_dst : bool, default True
        If True:
            - Bin/sum rows by dst_ids_raw.
            - Output row order = sorted unique dst IDs.
        If False:
            - Keep original row indices.
            - No aggregation by ID on dst axis.
    sort_src : bool, default True
        If True:
            - Bin/sum cols by src_ids_raw.
            - Output col order = sorted unique src IDs.
        If False:
            - Keep original col indices.
            - No aggregation by ID on src axis.

    Returns
    -------
    W_ids : csr_matrix
        Aggregated matrix:
          n_rows = n_unique(dst_ids_raw) if sort_dst else n_dst
          n_cols = n_unique(src_ids_raw) if sort_src else n_src
    dst_ids_out : ndarray
        If sort_dst:
            sorted unique dst IDs (len = n_rows).
        Else:
            dst_ids_raw as 1D array (len = n_rows), aligned with rows.
    src_ids_out : ndarray
        If sort_src:
            sorted unique src IDs (len = n_cols).
        Else:
            src_ids_raw as 1D array (len = n_cols), aligned with cols.
    
    This function is used to aggregate the overlap-area matrix by state vector IDs
    """
    W = W_csr.tocoo()

    # ----- destination axis (rows) -----
    if sort_dst:
        # Bin by ID; bins ordered by sorted unique IDs
        dst_ids_unique, dst_inv = np.unique(dst_ids_raw, return_inverse=True)
        r_bins = dst_inv[W.row]                # indices 0..n_dst_unique-1
        n_rows = dst_ids_unique.size
        dst_ids_out = dst_ids_unique
    else:
        # Keep as-is; assume dst_ids_raw aligned with rows
        r_bins = W.row
        n_rows = W.shape[0]
        dst_ids_out = np.asarray(dst_ids_raw)

    # ----- source axis (cols) -----
    if sort_src:
        src_ids_unique, src_inv = np.unique(src_ids_raw, return_inverse=True)
        c_bins = src_inv[W.col]                # indices 0..n_src_unique-1
        n_cols = src_ids_unique.size
        src_ids_out = src_ids_unique
    else:
        c_bins = W.col
        n_cols = W.shape[1]
        src_ids_out = np.asarray(src_ids_raw)

    # ----- build aggregated matrix -----
    W_ids = coo_matrix(
        (W.data.astype(np.float32, copy=False), (r_bins, c_bins)),
        shape=(n_rows, n_cols),
        dtype=np.float32,
    ).tocsr()

    return W_ids, dst_ids_out, src_ids_out

def aggregate_jacobian_along_dst_sv(jacobian_row_col, dst_ids_raw):
    """
    Aggregate a dense Jacobian along destination state vector IDs (column-wise),
    using a fully vectorized grouped reduction.

    Parameters
    ----------
    jacobian_row_col : np.ndarray, shape (n_dst_obs, n_dst_sv)
        Dense Jacobian; each column corresponds to a destination SV element
        at gridcell level (may include duplicates / clusters).
    dst_ids_raw : array_like, shape (n_dst_sv,)
        Destination state vector IDs (>0, 1-based ints, may repeat),
        aligned with the columns of jacobian_row_col.

    Returns
    -------
    J_agg : np.ndarray, shape (n_dst_obs, n_dst_sv_unique)
        Column-aggregated Jacobian. Each column is the sum over all original
        columns sharing that destination ID. Columns ordered by sorted unique IDs.
    dst_ids_unique : np.ndarray, shape (n_dst_sv_unique,)
        Sorted unique destination IDs corresponding to columns of J_agg.
    
    This function is used to aggregate the Jacobian matrix by destination state vector IDs
    """
    J = np.asarray(jacobian_row_col)
    ids = np.asarray(dst_ids_raw)

    if J.ndim != 2:
        raise ValueError(f"jacobian_row_col must be 2D, got shape {J.shape}")
    if J.shape[1] != ids.size:
        raise ValueError(
            f"dst_ids_raw must have same length as number of columns in jacobian_row_col "
            f"({ids.size} vs {J.shape[1]})."
        )

    # 1) Sort columns by ID so equal IDs are contiguous
    order = np.argsort(ids, kind="mergesort")  # stable; ids_sorted is sorted by ID
    ids_sorted = ids[order]
    J_sorted = J[:, order]

    # 2) Unique IDs and start indices of each group
    dst_ids_unique, idx_start = np.unique(ids_sorted, return_index=True)

    mask = np.isfinite(J_sorted)              # treat inf/-inf like NaN
    J_filled = np.where(mask, J_sorted, 0.0)

    J_sum   = np.add.reduceat(J_filled, idx_start, axis=1)
    # number of finite values for each row
    counts  = np.add.reduceat(mask.view(np.int8), idx_start, axis=1)

    # If a whole group was all-NaN, keep it NaN instead of 0
    J_sum[counts == 0] = np.nan
    J_agg = J_sum

    return J_agg, dst_ids_unique

def get_jacobian_scale(config, inv_directory, prior, sort_by_sv=False):
    # If there is state vector clustering:
    #   Get the jacobian scale on the grid cell level by default, or by state vector ID level if sort_by_sv=True
    
    # It is assumed that prior emission from the reference is already median and sorted by state vector IDs if clustering is used.
    # For the destination state vector, it is assumed that prior emission is gridded
    
    # For state vector clusters, we are using the median perturbation scale factor within the cluster,
    # which means the target emissions are not necessarily the stated target emissions in configuration file
    # To accommodate this, we get the real target_emission = perturbation scale factor * prior emission,
    # and then the scale to reconcile the jacobian column becomes:
    # scale = (perturbation scale factor) / (target emission) / (state vector area)
    #       = (perturbation scale factor) / (perturbation scale factor * prior emission) / (state vector area)
    #       = 1 / prior emission / (state vector area)
    # This works best when the reference precomputed Jacobian is calculated without any clustering, 
    # as we have to aggregate the Jacobian scale over the reference state vector.
    # If the reference precomputed Jacobian is calculated with clustering, we use median prior emissions
    # for the calculation of jacobian scale, which is an approximation.
    # It is OK to have destination state vector clustering, as we are applying the scale for each destination grid cell
    # and then aggregating the Jacobian column by destination state vector IDs afterwards.
    
    # get the area ratio between the reference grid cell and the destination grid cell
    RunDirs=f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}"
    CSgridDir=f"{RunDirs}/CS_grids"
    
    if config['UseGCHP']:
        sv_fpath = os.path.join(os.path.expandvars(inv_directory), 'StateVector.nc')
        with xr.open_dataset(sv_fpath).squeeze() as sv_ds:
            sv = sv_ds['StateVector'].values
            
        CS_RES = config['CS_RES']
        grid_path = f"{CSgridDir}/grids.c{CS_RES}.nc"
    else:
        sv_fpath = os.path.join(os.path.expandvars(inv_directory), 'StateVector.nc')
        with xr.open_dataset(sv_fpath).squeeze() as sv_ds:
            sv_lat = sv_ds['lat'].values
            sv_lon = sv_ds['lon'].values
            sv = sv_ds['StateVector'].values
            
        nlat = len(sv_lat)
        nlon = len(sv_lon)
        grid_path = f"{CSgridDir}/regular_lat_lon_{nlat}x{nlon}.nc"
    
    with xr.open_dataset(grid_path).squeeze() as grid_ds:
        area = grid_ds['area'].values.copy() # m2
        if sort_by_sv:
            # sum and sort area over state vector IDs
            area = sum_and_sort_along_statevector(area, sv, np.nan)
    
    # --- avoid division by zero ---
    # Identify where prior == 0
    zero_mask = (prior == 0) | np.isnan(prior)

    # Initialize sf_K as NaN everywhere first
    sf_K = np.full_like(prior, np.nan, dtype=np.float32)

    # Compute only for valid (nonzero) entries
    valid_mask = ~zero_mask
    sf_K[valid_mask] = 1.0 / prior[valid_mask] / area[valid_mask]

    sf_K = np.nan_to_num(sf_K, nan=1.0)
    
    return np.asarray(sf_K).ravel().astype(np.float32, copy=False)

def get_regrid_weights_jacobian_col(config, inv_directory, ref_config, ref_directory, ref_prior):
    
    """Getting overlap area between reference state vector and the destination state vector

    Returns:
        overlap_area: sparse matrix of overlap area
    """
    
    # get the regridding weights
    src_shape, dst_shape, regrid_weight_ds = get_regrid_weights_ESMF(config, inv_directory, ref_config)
    n_dst = int(np.prod(dst_shape))
    n_src = int(np.prod(src_shape))
    
    # ESMF rows=dst, cols=src (1-based to 0-based)
    src_ind = regrid_weight_ds['col'].values - 1
    dst_ind = regrid_weight_ds['row'].values - 1
    regrid_weights = regrid_weight_ds['S'].values.astype(np.float32, copy=False) # (n_dst, n_src)
    regrid_weights_csr = csr_matrix((regrid_weights, (dst_ind, src_ind)), 
                                    shape=(n_dst, n_src))
    
    dst_area = regrid_weight_ds['area_b'].values # (n_dst)
    
    # sort and sum value over state vector elements
    sv_fpath = os.path.join(os.path.expandvars(inv_directory), 'StateVector.nc')
    with xr.open_dataset(sv_fpath).squeeze() as sv_ds:
        sv = sv_ds['StateVector'].values
    flat_sv = sv.ravel(order="C")
    
    ref_sv_fpath = os.path.join(os.path.expandvars(ref_directory), 'StateVector.nc')
    with xr.open_dataset(ref_sv_fpath).squeeze() as ref_sv_ds:
        ref_sv = ref_sv_ds['StateVector'].values
    flat_ref_sv = ref_sv.ravel(order="C")
    
    dst_mask = flat_sv > 0
    src_mask = flat_ref_sv > 0
    
    W = regrid_weights_csr[dst_mask, :][:, src_mask] # (n_dst_sv × n_src_sv)
    
    # using overlap_area as the regridding weights for 
    # multiple overlapping reference state vector to destination single state vector
    overlap_area = W.multiply(dst_area[dst_mask][:, None]).tocsr()  # (n_dst_sv × n_src_sv)
    
    # sum and sort over state vector ID (Aggregate by IDs (sum over duplicate IDs when there is sv cluster))
    # Get corresponding state vector IDs for masked entries
    # IDs for valid entries
    dst_ids_raw = flat_sv[dst_mask].astype(np.int64, copy=False) # (n_dst_sv, )
    src_ids_raw = flat_ref_sv[src_mask].astype(np.int64, copy=False) # (n_src_sv, )

    # aggregate and sort overlap area by source state vector IDs as this is how reference jacobian columns are ordered
    # We defer the aggregation on destination state vector IDs to later, after jacobian ratio is applied
    # overlap_area in unit of square radians
    overlap_area_src, _, _ = sum_and_sort_along_statevector_matrix(
        overlap_area, dst_ids_raw, src_ids_raw, sort_dst=False, sort_src=True) # (n_dst_sv × n_src_sv_u)
    
    # For single destination grid cell:
    # jacobian ratio = (perturbation scale factor)_ref / (perturbation scale factor)_dst \
    #     / ( (target emission)_ref / (target emission)_dst ) \
    #     / ( (state vector area)_ref / (state vector area)_dst )
    # we will do row scale by 1 / dst_scale later outside the loop for each reference directory 
    # where this function of get_regrid_weights_jacobian_col is called.
    # 1 / dst_scale is applied during regrid_jacobian_row_col function
    src_scale = get_jacobian_scale(ref_config, ref_directory, ref_prior, sort_by_sv=True) # (n_src_sv_u,)
    
    # to make it memory-efficient
    # make a copy to keep overlap_area unchanges
    overlap_area_jacobian_ratio_src = overlap_area_src.copy() # (n_dst_sv × n_src_sv_u)
    overlap_area_jacobian_ratio_src = overlap_area_jacobian_ratio_src.tocsc()
    overlap_area_jacobian_ratio_src.data *= np.repeat(src_scale,
                                                  np.diff(overlap_area_jacobian_ratio_src.indptr))
    overlap_area_jacobian_ratio_src = overlap_area_jacobian_ratio_src.tocsr()

    
    return overlap_area_jacobian_ratio_src, overlap_area_src

def regrid_jacobian_row_col(
    jacobian_RegridRow,                  # J: (n_dst_obs, n_ref_sv_u)  dense np.ndarray
    overlap_area_jacobian_ratio_src,     # R: (n_dst_sv, n_ref_sv_u) sparse
    overlap_area_src,                    # W: (n_dst_sv, n_ref_sv_u) sparse
    inv_directory, config, prior
):
    """
    Compute weighted, jacobian-corrected averages per destination SV.

    For each destination grid cell {dsti}:
        jacobian_row_col[:, {dsti}] = sum_{refi} ( J[:, {refi}] * R[{dsti}, {refi}] ) / sum_{refi} W[{dsti}, {refi}]

    Args
    ----
    jacobian_RegridRow : np.ndarray, shape (n_dst_obs, n_ref_sv_u)
        Dense matrix of per-ref-SV Jacobian rows.
    overlap_area_jacobian_ratio_src : scipy.sparse.spmatrix, shape (n_dst_sv, n_ref_sv_u)
        Sparse matrix of the combined overlap area corrected by the jacobian ratio (row-area-scale by dst, col-area-scale by ref).
    overlap_area_src : scipy.sparse.spmatrix, shape (n_dst_sv, n_ref_sv_u)
        Sparse matrix of overlap area.

    Returns
    -------
    jacobian_row_col_sv : np.ndarray, shape (n_dst_obs, n_dst_sv_u)
        Weighted averages for each destination SV.
    """

    n_dstobs, n_ref_sv = jacobian_RegridRow.shape
    n_dst_sv_W, n_ref_sv_W = overlap_area_src.shape
    n_dst_sv_R, n_ref_sv_R = overlap_area_jacobian_ratio_src.shape

    # --- Sanity check for shape matching ---
    if (n_ref_sv_W != n_ref_sv) or (n_ref_sv_R != n_ref_sv) or (n_dst_sv_R != n_dst_sv_W):
        raise ValueError(
            f"Shape mismatch: jacobian_RegridRow{jacobian_RegridRow.shape}, "
            f"overlap_area_jacobian_ratio_src{overlap_area_jacobian_ratio_src.shape}, overlap_area_src{overlap_area_src.shape}. \
             Expected same n_ref_sv and n_dst_sv in overlap_area_jacobian_ratio_src/overlap_area_src."
        )

    # --- ensure good sparse formats & dtypes ---
    # Row ops fast in CSR; column ops on R.T fast in CSC
    R = overlap_area_jacobian_ratio_src.tocsr()
    # scale each row by 1/dst_scale,
    # sort and sum value over state vector elements
    sv_fpath = os.path.join(os.path.expandvars(inv_directory), 'StateVector.nc')
    with xr.open_dataset(sv_fpath).squeeze() as sv_ds:
        sv = sv_ds['StateVector'].values
        sv_lat = sv_ds['lat'].values
        sv_lon = sv_ds['lon'].values
        sv_lon[sv_lon>180] -= 360

    flat_sv = sv.ravel(order="C")
    dst_mask = flat_sv > 0
    # get dst_scale for all grid cells first
    dst_scale = get_jacobian_scale(config, inv_directory, prior, sort_by_sv=False)     # (n_dst,)
    # subset to destination state vector (with potential duplicates due to clustering)
    dst_scale = dst_scale[dst_mask] # (n_dst_sv, )
    
    R.data *= np.repeat(1.0 / dst_scale, np.diff(R.indptr))
    
    W = overlap_area_src.tocsr()
    # Ensure a common float32 dtype to save memory
    if jacobian_RegridRow.dtype != np.float32:
        J = jacobian_RegridRow.astype(np.float32, copy=False)
    else:
        J = jacobian_RegridRow
    
    # --- numerator: J @ R.T ---
    # (n_dst_obs, n_ref_sv_u) @ (n_ref_sv_u, n_dst_sv) = (n_dst_obs, n_dst_sv)
    # Use CSC (Compressed Sparse Column) for efficient right-transpose multiply
    # best for operations like dense_matrix (J) @ sparse_matrix (R) (where R’s columns are multiplied)
    num = J @ R.T.tocsc()

    # --- denominator: sum of overlap area per destination SV ---
    denom = np.asarray(W.sum(axis=1)).ravel()  # (n_dst_sv,)

    # Avoid divide-by-zero: neutralize zero-weight columns
    denom_safe = denom.copy()
    zero_cols = (denom_safe == 0)
    if np.any(zero_cols):
        denom_safe[zero_cols] = 1.0

    # Weighted average in shape of (n_dst_obs, n_dst_sv)
    jacobian_row_col = num / denom_safe[None, :]  # broadcast over columns

    # nan columns with zero weights (no contributors)
    if np.any(zero_cols):
        jacobian_row_col[:, zero_cols] = np.nan

    # Ensure dense ndarray output
    if not isinstance(jacobian_row_col, np.ndarray):
        jacobian_row_col = np.asarray(jacobian_row_col)
    if jacobian_row_col.dtype != np.float32:
        jacobian_row_col = jacobian_row_col.astype(np.float32, copy=False)

    # sum and sort over destination state vector IDs (n_dst_obs, n_dst_sv) -> (n_dst_obs, n_dst_sv_u)
    # Get corresponding state vector IDs for masked entries
    # IDs for valid entries
    # Note: dst IDs may have duplicates (clusters)
    dst_ids_raw = flat_sv[dst_mask].astype(np.int64, copy=False) # (n_dst_sv, )
    
    # remove the one more grid cell along each domain edge
    res = config['Res']
    sv_lonm, sv_latm = np.meshgrid(sv_lon, sv_lat)
    sv_lonm = sv_lonm.ravel(order="C")[dst_mask]
    sv_latm = sv_latm.ravel(order="C")[dst_mask]
    # initialize shaved off degrees
    degx = 0
    degy = 0
    if config['isRegional']:
        if "0.125x0.15625" in res:
            degx = 4 * 0.15625
            degy = 4 * 0.125
        elif "0.25x0.3125" in res:
            degx = 4 * 0.3125
            degy = 4 * 0.25
        elif "0.5x0.625" in res:
            degx = 4 * 0.625
            degy = 4 * 0.5
    xlim = [sv_lon.min() + degx, sv_lon.max() - degx]
    ylim = [sv_lat.min() + degy, sv_lat.max() - degy]
    # set to nan for grid cells within the [4 4 4 4] grid cells along each regional domain edge
    indices = np.where(
        (sv_lonm < xlim[0]) | (sv_lonm > xlim[1]) |
        (sv_latm < ylim[0]) | (sv_latm > ylim[1])
    )[0]
    jacobian_row_col[:, indices] = np.nan
    jacobian_row_col_sv, _ = aggregate_jacobian_along_dst_sv(jacobian_row_col, dst_ids_raw)
    
    return jacobian_row_col_sv # (n_dst_obs, n_dst_sv_u)