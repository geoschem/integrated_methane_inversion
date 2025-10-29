import os
import numpy as np
import xarray as xr
from scipy.sparse import csr_matrix, coo_matrix, spdiags
from src.utilities.regrid_state_vector_file import(
    get_gridspec_prefix,
)

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
    
    if len(gc_lat) != len(sv_lat):
        ind_lat = np.abs(sv_lat[:, None] - gc_lat[lati]).argmin(axis=0)
        ind_lon = np.abs(sv_lon[:, None] - gc_lon[loni]).argmin(axis=0)
        superobs_index = np.ravel_multi_index((ind_lat, ind_lon), sv_shape)
    else:
        # grids identical
        superobs_index = GC_index
    
    return superobs_index
    
def get_regrid_weights_jacobian_row(config, inv_directory, GC_index, ref_config, ref_GC_index):
    # The row of Jacobian matrix K is super observations, 
    # or satellite observatiosn sampled at simulation grid
    
    # GC_index: flattened simulation index with valid super observations
    # note that GC_index for GCC is from the simulation diagnostics, 
    # it is not necessarily the same as the state vector netCDF file 
    # as state vector at GCC is trimed;
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
    
    return W

def sum_and_sort_along_statevector(val, sv, fill_value=np.nan):
    if np.isnan(fill_value):
        mask = ~np.isnan(sv)
    else:
        mask = sv != fill_value
    sv_valid = sv[mask]
    val_valid = val[mask]
    _, inv = np.unique(sv_valid, return_inverse=True)
    sv_val = np.bincount(inv, weights=val_valid)
    return sv_val

def get_jacobian_scale(period_number, config, inv_directory):
    # scale = (perturbation scale factor) / (target emission) / (state vector area)
    pert_sf_path = os.path.join(
        inv_directory, "archive_perturbation_sfs", f"pert_sf_{period_number}.npz"
    )
    pert_sf_dict = np.load(pert_sf_path)
    pert_sf = pert_sf_dict["effective_pert_sf"]
    
    # Get the ratio of the targetted emissions in the target and reference inversions
    target_emis = pert_sf_dict["target_emission"]
    
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
        area = grid_ds['area'].values # m2
        sv_area = sum_and_sort_along_statevector(area, sv, np.nan)
    
    sf_K = pert_sf / target_emis / sv_area
    
    return np.asarray(sf_K)

def get_regrid_weights_jacobian_col(period_number, config, inv_directory, ref_config, ref_directory):
    
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
    regrid_weights = regrid_weight_ds['S'].values # (n_dst, n_src)
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
    dst_ids_raw = flat_sv[dst_mask].astype(np.int64)
    src_ids_raw = flat_ref_sv[src_mask].astype(np.int64)

    # Map to sorted unique IDs and get per-row/col bin indices
    dst_ids_sorted, dst_inv = np.unique(dst_ids_raw, return_inverse=True)  # len = n_dst_sv_unique
    src_ids_sorted, src_inv = np.unique(src_ids_raw, return_inverse=True)  # len = n_src_sv_unique
    n_dst_sv_u = dst_ids_sorted.size
    n_src_sv_u = src_ids_sorted.size

    # One-hot aggregation matrices (gridcell -> ID bin)
    # Shapes: (n_dst_sv, n_dst_sv_unique) and (n_src_sv, n_src_sv_unique)
    A_dst = coo_matrix(
        (np.ones_like(dst_inv, dtype=np.float64),
         (np.arange(dst_inv.size), dst_inv)),
        shape=(dst_inv.size, n_dst_sv_u),
    ).tocsr()

    A_src = coo_matrix(
        (np.ones_like(src_inv, dtype=np.float64),
         (np.arange(src_inv.size), src_inv)),
        shape=(src_inv.size, n_src_sv_u),
    ).tocsr()

    # Aggregate by IDs: A_dst.T @ overlap_area @ A_src
    overlap_area_sv = (A_dst.T @ overlap_area @ A_src).tocsr()  # (n_dst_sv_unique, n_src_sv_unique)
    
    # reconcile jacobian sensitivities across grid
    # Running sensitivity inversions with pre-constructed Jacobian may require scaling the
    # Jacobian to match updated prior estimates. This is because the Jacobian is defined
    # by relative perturbations of emissions, so that the sensitivities depend on
    # the magnitude of emissions in the prior inventory.
    
    # For single state vector element:
    # jacobian ratio = (perturbation scale factor)_ref / (perturbation scale factor)_dst \
    #     / ( (target emission)_ref / (target emission)_dst ) \
    #     / ( (state vector area)_ref / (state vector area)_dst )
    dst_scale = get_jacobian_scale(period_number, config, inv_directory)     # (n_dst_sv_unique,)
    src_scale = get_jacobian_scale(period_number, ref_config, ref_directory) # (n_src_sv_unique)
    
    # sanity check
    assert dst_scale.shape == (n_dst_sv_u,)
    assert src_scale.shape == (n_src_sv_u,)

    # multiply cols by src_scale, and divide rows by dst_scale
    # Using the diagnostic matrices should give the same results and more math-beautiful, 
    # but this broadcasting multiplication method would be faster than using the sparse diagnostic matrices
    overlap_area_jacobian_ratio = (
        overlap_area_sv
        .multiply(src_scale[None, :])             # scale columns
        .multiply((1.0 / dst_scale)[:, None])     # scale rows
        .tocsr()
    )
    
    # # Build (very) sparse diagonals
    # dst_scale_diag  = spdiags(1.0 / dst_scale, 0, n_dst_sv_u, n_dst_sv_u)   # divide rows by dst
    # src_scale_diag = spdiags(src_scale,       0, n_src_sv_u, n_src_sv_u)    # multiply cols by ref

    # # Order of ops matters for performance but result is the same
    # overlap_area_jacobian_ratio = (dst_scale_diag @ overlap_area_sv @ src_scale_diag).tocsr() # (n_dst_sv_unique × n_src_sv_unique)
    
    return overlap_area_jacobian_ratio, overlap_area_sv

def regrid_jacobian_row_col(
    jacobian_RegridRow,              # J: (n_dst_obs, n_ref_sv)  dense np.ndarray
    overlap_area_jacobian_ratio,     # R: (n_dst_sv, n_ref_sv) sparse
    overlap_area_sv,                 # W: (n_dst_sv, n_ref_sv) sparse
):
    """
    Compute weighted, jacobian-corrected averages per destination SV.

    For each destination SV {dsti}:
        jacobian_row_col[:, {dsti}] = sum_{refi} ( J[:, {refi}] * R[{dsti}, {refi}] ) / sum_{refi} W[{dsti}, {refi}]

    Args
    ----
    jacobian_RegridRow : np.ndarray, shape (n_dst_obs, n_ref_sv)
        Dense matrix of per-ref-SV Jacobian rows.
    overlap_area_jacobian_ratio : scipy.sparse.spmatrix, shape (n_dst_sv, n_ref_sv)
        Sparse matrix of the combined overlap area corrected by the jacobian ratio (row-area-scale by dst, col-area-scale by ref).
    overlap_area_sv : scipy.sparse.spmatrix, shape (n_dst_sv, n_ref_sv)
        Sparse matrix of overlap area.

    Returns
    -------
    jacobian_row_col : np.ndarray, shape (n_dst_obs, n_dst_sv)
        Weighted averages for each destination SV.
    """

    n_dstobs, n_ref_sv = jacobian_RegridRow.shape
    n_dst_sv_W, n_ref_sv_W = overlap_area_sv.shape
    n_dst_sv_R, n_ref_sv_R = overlap_area_jacobian_ratio.shape

    # --- Sanity check for shape matching ---
    if (n_ref_sv_W != n_ref_sv) or (n_ref_sv_R != n_ref_sv) or (n_dst_sv_R != n_dst_sv_W):
        raise ValueError(
            f"Shape mismatch: jacobian_RegridRow{jacobian_RegridRow.shape}, "
            f"overlap_area_jacobian_ratio{overlap_area_jacobian_ratio.shape}, overlap_area_sv{overlap_area_sv.shape}. \
             Expected same n_ref_sv and n_dst_sv in overlap_area_jacobian_ratio/overlap_area_sv."
        )

    # --- ensure good sparse formats & dtypes ---
    # Row ops fast in CSR; column ops on R.T fast in CSC
    R = overlap_area_jacobian_ratio.tocsr()
    W = overlap_area_sv.tocsr()
    # Ensure a common float32 dtype to save memory
    if jacobian_RegridRow.dtype != np.float32:
        J = jacobian_RegridRow.astype(np.float32, copy=False)
    else:
        J = jacobian_RegridRow
    
    # --- numerator: J @ R.T ---
    # (n_dst_obs, n_ref_sv) @ (n_ref_sv, n_dst_sv) = (n_dst_obs, n_dst_sv)
    # Use CSC (Compressed Sparse Column) for efficient right-transpose multiply
    # best for operations like dense_matrix (J) @ sparse_matrix (R) (where R’s columns are multiplied)
    num = J @ R.T.tocsc()

    # --- denominator: sum of weights per destination SV ---
    denom = np.asarray(W.sum(axis=1)).ravel()  # (n_dst_sv,)

    # Avoid divide-by-zero: neutralize zero-weight columns
    denom_safe = denom.copy()
    zero_cols = (denom_safe == 0)
    if np.any(zero_cols):
        denom_safe[zero_cols] = 1.0

    # Weighted average
    jacobian_row_col = num / denom_safe[None, :]  # broadcast over columns

    # nan columns with zero weights (no contributors)
    if np.any(zero_cols):
        jacobian_row_col[:, zero_cols] = np.nan

    # Ensure dense ndarray output
    if not isinstance(jacobian_row_col, np.ndarray):
        jacobian_row_col = np.asarray(jacobian_row_col)
    if jacobian_row_col.dtype != np.float32:
        jacobian_row_col = jacobian_row_col.astype(np.float32, copy=False)

    return jacobian_row_col