import numpy as np
import xarray as xr
import os
import subprocess
import pygeohash as pgh
import yaml
from scipy.sparse import csc_matrix
from src.utilities.make_state_vector_file import(
    make_state_vector_file,
)

def _format_stretch_factor(sf: float) -> str:
    """
    Format stretch factor like the Bash version:
    1.5 -> '1d50' (i.e., two decimals, then '.' -> 'd').
    """
    s = f"{sf:.2f}"
    return s.replace(".", "d")


def get_gridspec_prefix(
    cs_res,
    stretch_grid,
    stretch_factor: float = 1.0,
    target_lat: float = -90.0,
    target_lon: float = 170.0,
    geohash_precision=12
) -> str:
    """
    Build the GridSpec prefix.

    Examples
    --------
    >>> get_gridspec_prefix(48, False)
    'c48'
    >>> get_gridspec_prefix(36, True, 1.5, 32.0, -103.0)
    'c36_s1d50_t<geohash>'
    """
    if not isinstance(cs_res, int) or cs_res <= 0:
        raise ValueError("cs_res must be a positive integer.")

    if not stretch_grid:
        return f"c{cs_res}"

    sf_formatted = _format_stretch_factor(float(stretch_factor))
    target_lon = ((float(target_lon) + 180) % 360) - 180

    geoh = pgh.encode(float(target_lat), target_lon, precision=geohash_precision)
    return f"c{cs_res}_s{sf_formatted}_t{geoh}"

def subset_ensemble_sv(grid_sv_ds, all_ref_sv_ds, save_subset_sv_path):
    region_lat = grid_sv_ds['lat'].values
    region_lon = grid_sv_ds['lon'].values
    
    if config["isRegional"]:
        if config["Res"] == "4.0x5.0":
            deg_lat, deg_lon = 4.0, 5.0
        elif config["Res"] == "2.0x2.5":
            deg_lat, deg_lon = 2.0, 2.5
        elif config["Res"] == "0.5x0.625":
            deg_lat, deg_lon = 0.5, 0.625
        elif config["Res"] == "0.25x0.3125":
            deg_lat, deg_lon = 0.25, 0.3125
        elif config["Res"] == "0.125x0.15625":
            deg_lat, deg_lon = 0.125, 0.15625
        buffer_deg_lat = deg_lat * 3
        buffer_deg_lon = deg_lon * 3

    extents = [region_lon.min() - buffer_deg_lon, region_lon.max() + buffer_deg_lon,
               region_lat.min() - buffer_deg_lat, region_lat.max() + buffer_deg_lat]
    
    # get target face indices with stretched face (index 5) within the extents
    ref_lats = all_ref_sv_ds['lats'].values[:,5,...]
    ref_lons = np.array(all_ref_sv_ds['lons'][:,5,...])
    ref_lons[ref_lons>180] -= 360
    ind = np.where((ref_lons>=extents[0]) & (ref_lons<=extents[1]) & 
                   (ref_lats>=extents[2]) & (ref_lats<=extents[3]))
    target_faces = np.unique(ind[0])
    subset_ref_sv_ds = all_ref_sv_ds.isel(target_face=list(target_faces))
    subset_ref_sv_ds = subset_ref_sv_ds.assign_coords(
        target_face=all_ref_sv_ds['target_face'].isel(target_face=list(target_faces))
    )
    if save_subset_sv_path is not None:
        print("Saving subset ensemble state vector file {}".format(save_subset_sv_path))
        subset_ref_sv_ds.to_netcdf(
            save_subset_sv_path,
            encoding={
                v: {"zlib": True, "complevel": 1} for v in subset_ref_sv_ds.data_vars
            },
        )
    return subset_ref_sv_ds

def regrid_state_vector_file(config, grid_sv_ds):
    RunDirs=f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}"
    CSgridDir=f"{RunDirs}/CS_grids"
    
    ref_sv_fpath = config['StateVectorFile']
    # do not interpret fillvalue as nan, so that regridded state vector 
    # would have value when it contains nan values but not all as nan
    with xr.open_dataset(ref_sv_fpath, mask_and_scale=False) as tmp:
        all_ref_sv_ds = tmp.load()
    if config['isRegional']:
        ref_sv_fname = os.path.basename(ref_sv_fpath).replace('combined', 'subset')
        save_subset_sv_path = f"{CSgridDir}/{ref_sv_fname}"
        ref_sv_ds = subset_ensemble_sv(grid_sv_ds, all_ref_sv_ds, save_subset_sv_path)
    else:
        ref_sv_ds = all_ref_sv_ds
    CS_RES = int(ref_sv_ds.attrs['CS_RES'])
    STRETCH_FACTOR = float(ref_sv_ds.attrs['STRETCH_FACTOR'])
    STRETCH_GRID = (STRETCH_FACTOR - 1.0) > 1e-2
    TARGET_LATs = np.array(ref_sv_ds['TARGET_LAT'])
    TARGET_LONs = np.array(ref_sv_ds['TARGET_LON'])
    TARGET_LONs[TARGET_LONs>180] -= 360
    src_shape = (6, CS_RES, CS_RES)
    n_src = int(np.prod(src_shape))
    
    if config['UseGCHP']:
        lats = grid_sv_ds['lats'].values
        lons = grid_sv_ds['lons'].values
    else:
        lat = grid_sv_ds['lat'].values
        lon = grid_sv_ds['lon'].values
        lons, lats = np.meshgrid(lon, lat)
    dst_shape = lats.shape
    n_dst = int(np.prod(dst_shape))
    
    if config['UseGCHP']:
        dst_CS_RES = config['CS_RES']
        dst_prefix = get_gridspec_prefix(dst_CS_RES, False)
        dst_grid_path = f"{dst_prefix}_gridspec.nc"
        if not os.path.exists(os.path.join(CSgridDir, dst_grid_path)):
            subprocess.run([
                'gridspec-create', 'gcs',
                f'{dst_CS_RES}'
                ],
                check=True, cwd=CSgridDir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,)
    else:
        dst_lat = grid_sv_ds['lat'].values
        dst_lon = grid_sv_ds['lon'].values
        nlat = len(dst_lat)
        nlon = len(dst_lon)
        delta_lon = np.median(np.diff(dst_lon))
        dst_grid_path = f"regular_lat_lon_{nlat}x{nlon}.nc"
        if not os.path.exists(os.path.join(CSgridDir, dst_grid_path)):
            subprocess.run([
                'gridspec-create', 'latlon',
                '-b', f'{dst_lon.min()}', f'{dst_lat.min()}',
                f'{dst_lon.max()+delta_lon}', f'{dst_lat.max()}',
                '-pc', '-hp', '-dc',
                f'{nlat}', f'{nlon}',
                ],
                check=True, cwd=CSgridDir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,)
    
    # collect per-face regrid matrices after src masking and sum over the target face
    # (n_dst, 1) dense matrix accumulator
    regrid_weights_row_sum = np.zeros((n_dst,), dtype=np.float32)
    for i in range(len(TARGET_LATs)):
        TARGET_LAT = TARGET_LATs[i]
        TARGET_LON = TARGET_LONs[i]
        prefix = get_gridspec_prefix(CS_RES, STRETCH_GRID, STRETCH_FACTOR, TARGET_LAT, TARGET_LON)
        if not os.path.exists(os.path.join(CSgridDir, f"{prefix}_gridspec.nc")):
            subprocess.run([
                'gridspec-create', 'sgcs',
                '-s', f"{STRETCH_FACTOR}",
                '-t', f"{TARGET_LAT}", f"{TARGET_LON}", f"{CS_RES}"
                ],
                check=True, cwd=CSgridDir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,)

        src_grid_path = f"{prefix}_gridspec.nc"
        if config['UseGCHP']:
            regrid_weight_fpath = f"regrid_weights_c{CS_RES}_s{STRETCH_FACTOR:.1f}_{TARGET_LAT:.1f}N_{TARGET_LON:.1f}E_to_c{dst_CS_RES}_conserve.nc"
        else:
            regrid_weight_fpath = f"regrid_weights_c{CS_RES}_s{STRETCH_FACTOR:.1f}_{TARGET_LAT:.1f}N_{TARGET_LON:.1f}E_to_latlon_{nlat}x{nlon}_conserve.nc"
        if not os.path.exists(os.path.join(CSgridDir, regrid_weight_fpath)):
            subprocess.run([
                "ESMF_RegridWeightGen",
                "-s", src_grid_path,
                "-d", dst_grid_path,
                "-m", "conserve",
                "--ignore_unmapped",
                "-w", regrid_weight_fpath
                ], 
                check=True, cwd=CSgridDir, 
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL, 
            )

        with xr.open_dataset(os.path.join(CSgridDir, regrid_weight_fpath)) as regrid_weight_ds:
            # ESMF rows=dst, cols=src (1-based to 0-based)
            src_ind = regrid_weight_ds['col'].values - 1
            dst_ind = regrid_weight_ds['row'].values - 1
            regrid_weights = regrid_weight_ds['S'].values 
        regrid_weights_csc = csc_matrix((regrid_weights, (dst_ind, src_ind)), 
                                        shape=(n_dst, n_src)) # (n_dst, n_src) sparse matrix
        
        ref_sv_ds_face = ref_sv_ds.isel(target_face=i).squeeze()
        ref_sv_face = ref_sv_ds_face['StateVector'].values
        flat_ref_sv_face = ref_sv_face.ravel(order="C")
        
        src_mask = flat_ref_sv_face > 0
        # sum over each target face where state vector ID > 0
        W = np.asarray(regrid_weights_csc[:, src_mask].sum(axis=1)).squeeze().astype('float32') # (n_dst, ) dense matrix
        # Sum regrid_weights along the target face
        # to get total contribution from ensemble target faces per destination
        regrid_weights_row_sum += W

    # Find destination indices with any non-zero total weight and 
    # valid state vector IDs when creating normal state vector at the destination grid
    grid_sv = grid_sv_ds['StateVector'].values
    regrid_sv_mask = ((regrid_weights_row_sum > 0).reshape(dst_shape)) & (grid_sv > 0) & (~np.isnan(grid_sv))

    if config['isRegional']:
        if config['UseGCHP']:
            raise ValueError("Regional masking currently implemented only for regular lat-lon output.")
        lat_min = config["LatMin"]
        lat_max = config["LatMax"]
        lon_min = config["LonMin"]
        lon_max = config["LonMax"]
        lonm, latm = np.meshgrid(dst_lon, dst_lat)
        valid_mask_indomain = (regrid_sv_mask) & (lonm>=lon_min) & (lonm<=lon_max) & (latm>=lat_min) & (latm<=lat_max)
        valid_mask_buffer = regrid_sv_mask & ~valid_mask_indomain
        
        grid_sv_new = np.empty(dst_shape)
        grid_sv_new[:] = np.nan

        if valid_mask_indomain.any():
            uniq, inv = np.unique(grid_sv[valid_mask_indomain], return_inverse=True)  # uniq is sorted ascending
            N = uniq.size
            grid_sv_new[valid_mask_indomain] = inv + 1                           # 1..N in same order as uniq
        else:
            N = 0

        # Shift buffer labels so they come right after in-domain, preserving their order
        # (uses the max original in-domain label to compute a constant offset)
        if valid_mask_buffer.any():
            if N > 0:
                max_in_domain_orig = grid_sv[valid_mask_indomain].max()
                offset = max_in_domain_orig - N
            else:
                # No in-domain cells: keep buffer labels as-is
                offset = 0
            grid_sv_new[valid_mask_buffer] = grid_sv[valid_mask_buffer] - offset
        
    else:
        valid_mask = regrid_sv_mask

        # zeros as sentinel
        grid_sv_new = np.empty(dst_shape)
        grid_sv_new[:] = np.nan

        if valid_mask.any():
            # uniq sorted ascending; inverse gives compact IDs in that order
            _, inv = np.unique(grid_sv[valid_mask], return_inverse=True)
            grid_sv_new[valid_mask] = inv + 1  # 1..N

    
    refyear = 2000
    fillvalue = -9999
    # Make dataset
    if config['UseGCHP']:
        gridfpath=f"{CSgridDir}/grids.c{config['CS_RES']}.nc"
        gridds = xr.open_dataset(gridfpath)
    
        # add time dimension
        da_statevector = xr.DataArray(grid_sv_new[None,...], dims=['time', 'nf', 'Ydim', 'Xdim'],
                                    coords=dict(time=(['time'], [0.]), lats=(['nf', 'Ydim', 'Xdim'], gridds['lats'].values),
                                                lons=(['nf', 'Ydim', 'Xdim'], gridds['lons'].values)),
                                    attrs=dict(units='1', missing_value=fillvalue, _FillValue=fillvalue))
        ds_statevector = xr.Dataset({'StateVector': da_statevector})

        # Add attribute metadata
        ds_statevector.lats.attrs["units"] = "degrees_north"
        ds_statevector.lats.attrs["long_name"] = "Latitude"
        ds_statevector.lons.attrs["units"] = "degrees_east"
        ds_statevector.lons.attrs["long_name"] = "Longitude"
        ds_statevector['time'].attrs = dict(units='days since {}-01-01 00:00:00'.format(refyear),
                                            delta_t='0000-01-00 00:00:00', axis='T', standard_name='Time',
                                            long_name='Time', calendar='standard')
        ds_statevector['corner_lats'] = gridds['corner_lats']
        ds_statevector['corner_lons'] = gridds['corner_lons']
    else:
        da_statevector = xr.DataArray(grid_sv_new[None,...], dims=['time', 'lat', 'lon'],
                                    coords=dict(time=(['time'], [0.]), lat=(['lat'], grid_sv_ds['lat'].values),
                                                lon=(['lon'], grid_sv_ds['lon'].values)),
                                    attrs=dict(units='1', missing_value=fillvalue, _FillValue=fillvalue))
        ds_statevector = xr.Dataset({'StateVector': da_statevector})

        # Add attribute metadata
        ds_statevector['time'].attrs = dict(units='days since {}-01-01 00:00:00'.format(refyear),
                                            delta_t='0000-01-00 00:00:00', axis='T', standard_name='Time',
                                            long_name='Time', calendar='standard')
        ds_statevector.lat.attrs["units"] = "degrees_north"
        ds_statevector.lat.attrs["long_name"] = "Latitude"
        ds_statevector.lon.attrs["units"] = "degrees_east"
        ds_statevector.lon.attrs["long_name"] = "Longitude"

    return ds_statevector

if __name__ == "__main__":
    import sys

    config_path = sys.argv[1]
    land_cover_pth = sys.argv[2]
    hemco_diag_pth = sys.argv[3]
    save_pth = sys.argv[4]

    grid_sv_ds = make_state_vector_file(
        config_path,
        land_cover_pth,
        hemco_diag_pth,
    ).squeeze()
    
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    
    ds_statevector = regrid_state_vector_file(config, grid_sv_ds)
    
    if save_pth is not None:
        if os.path.exists(save_pth):
            os.remove(save_pth)
        print("Saving file {}".format(save_pth))
        ds_statevector.to_netcdf(
            save_pth,
            encoding={
                v: {"zlib": True, "complevel": 1} for v in ds_statevector.data_vars
            },
        )