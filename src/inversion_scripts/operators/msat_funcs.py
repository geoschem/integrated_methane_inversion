from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from scipy.constants import gas_constant as Rgas
from scipy.stats import binned_statistic_2d
from scipy.constants import N_A
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import datetime
#from operator_utilities import *
import os
#import interpolation
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
import xarray as xr
from src.inversion_scripts.utils import get_strdate


# def make_custom_bc_files():



def adjust_column_average_concentration_time_avg(met_file, conc_file, target_column_avg, output_file):
    """
    Adjust the column-average concentration of a chemical species to a uniform value,
    using time-averaged data for calculations but outputting data with the original time dimension.

    Parameters:
    met_file (str): Path to the meteorological file (NetCDF).
    conc_file (str): Path to the concentration file (NetCDF).
    target_column_avg (float): Target column-average concentration (units must match original data).
    output_file (str): Path to save the modified concentration file.

    Returns:
    None: Saves the modified concentration file to output_file.
    """
    # Open the files
    met_ds = xr.open_dataset(met_file)
    conc_ds = xr.open_dataset(conc_file)
    
    # Extract relevant variables
    air_density = met_ds['Met_AIRDEN']  # Air density
    air_volume = met_ds['Met_AIRVOL']  # Air volume
    species_conc = conc_ds['SpeciesBC_CH4']  # Species concentration
    
    # Time-average the air density, air volume, and species concentration
    air_density_avg = air_density.mean(dim='time')
    air_volume_avg = air_volume.mean(dim='time')
    species_conc_avg = species_conc.mean(dim='time')
    
    # Calculate the time-averaged column-average concentration
    column_avg = (species_conc_avg * air_density_avg * air_volume_avg).sum(dim='lev') / (air_density_avg * air_volume_avg).sum(dim='lev')

    
    # Compute scaling factor to achieve target column average
    scaling_factor = target_column_avg / column_avg
    
    # Apply scaling to the time-averaged concentration
    adjusted_conc_avg = species_conc_avg * scaling_factor
    
    # Broadcast the adjusted concentration to all time steps in the original dataset
    adjusted_conc = adjusted_conc_avg.expand_dims(dim={"time": conc_ds.dims['time']}).copy()
    
    # Save adjusted concentrations into a new file
    conc_ds['SpeciesBC_CH4'][:] = adjusted_conc
    conc_ds.to_netcdf(output_file)
    print(f"Modified concentration file saved to {output_file}")

def process_files_in_date_range(start_date, end_date, target_column_avg, conc_path, met_path, output_dir):
    """
    Adjust column-average concentration for all files within a given date range.

    Parameters:
    start_date (str): Start date in the format 'YYYY-MM-DD'.
    end_date (str): End date in the format 'YYYY-MM-DD'.
    target_column_avg (float): Target column-average concentration.
    conc_template (str): Template path for concentration files with $YYYY, $MM, $DD.
    met_template (str): Template path for meteorological files with $YYYY, $MM, $DD.
    output_dir (str): Directory to save the adjusted concentration files.

    Returns:
    None
    """
    # Parse start and end dates
    start_date = datetime.datetime.strptime(start_date, "%Y%m%d")
    end_date = datetime.datetime.strptime(end_date, "%Y%m%d")
    start_date-= datetime.timedelta(days=31)
    
    current_date = start_date

    met_template = "GEOSChem.StateMet.$YYYY$MM$DD_0000z.nc4"
    conc_template = "GEOSChem.BoundaryConditions.$YYYY$MM$DD_0000z.nc4"
    while current_date <= end_date:
        # Format date components
        yyyy = current_date.strftime("%Y")
        mm = current_date.strftime("%m")
        dd = current_date.strftime("%d")
        
        # Construct file paths
        conc_file = conc_template.replace("$YYYY", yyyy).replace("$MM", mm).replace("$DD", dd)
        met_file = met_template.replace("$YYYY", yyyy).replace("$MM", mm).replace("$DD", dd)
        output_file = os.path.join(output_dir, f"GEOSChem.BoundaryConditions.{yyyy}{mm}{dd}_0000z.nc4")

        conc_file = f"{conc_path}/{conc_file}"
        met_file = f"{met_path}/{met_file}"

        print(met_file, os.path.exists(met_file))
        print(conc_file, os.path.exists(conc_file))
        # Check if input files exist
        if os.path.exists(conc_file) and os.path.exists(met_file):
            print(f"Processing date: {current_date.strftime('%Y-%m-%d')}")
            adjust_column_average_concentration_time_avg(met_file, conc_file, target_column_avg, output_file)
        else:
            print(f"Skipping date {current_date.strftime('%Y-%m-%d')}: Input files not found.")
        
        # Move to the next date
        current_date += datetime.timedelta(days=1)


def get_sv_ll_bnds(state_vector):
    sv_ll = {}
    
    sv = xr.open_dataset(state_vector)
    lat_bnds = np.zeros(len(sv['lat']+1))
    lon_bnds = np.zeros(len(sv['lon']+1))
    
        
    lat_inc = (sv['lat'][1]-sv['lat'][0])/2
    lon_inc = (sv['lon'][1]-sv['lon'][0])/2


    print(lon_inc)
    for j in range(len(lat_bnds)-1):
        lat_bnds[j] = sv['lat'][j]-lat_inc
    for j in range(len(lon_bnds)-1):
        lon_bnds[j] = sv['lon'][j]-lon_inc

    lat_bnds[-1] = sv['lat'][-1]+lat_inc
    lon_bnds[-1] = sv['lon'][-1]+lon_inc

    sv_ll['lon_bnds'] = lon_bnds
    sv_ll['lat_bnds'] = lat_bnds
    

    return sv_ll

def get_sv_ll(state_vector):
    sv_ll = {}
    
    sv = xr.open_dataset(state_vector)


    sv_ll['lon'] = sv['lon']
    sv_ll['lat'] = sv['lat']
    

    return sv_ll

def _legacy_average_methanesat_observations(
    data_path,
    state_vector,
    species="CH4",
    time_threshold=None,
    file_id='MSAT',
):
    
    directory_path = data_path

    
    gc_ll_bnds = get_sv_ll_bnds(state_vector)

    lat_edges_list = gc_ll_bnds['lat_bnds']
    lon_edges_list = gc_ll_bnds['lon_bnds']
    



    xch4_out = np.zeros((len(lat_edges_list)-1, len(lon_edges_list)-1))
    ak_ch4_out = np.zeros((len(lat_edges_list)-1, len(lon_edges_list)-1,19))
    p_surf_out = np.zeros((len(lat_edges_list)-1, len(lon_edges_list)-1))
    p_trop_out = np.zeros((len(lat_edges_list)-1, len(lon_edges_list)-1))
    prior_out = np.zeros((len(lat_edges_list)-1, len(lon_edges_list)-1,19))
    c_out = np.zeros((len(lat_edges_list)-1, len(lon_edges_list)-1))
    time_out = np.zeros((len(lat_edges_list)-1, len(lon_edges_list)-1))

    file_num = 0

    class RestartLoop(Exception):
        pass


    file_path = data_path

    print(file_num)
    file_num  = file_num+1

    #small snippet for testing if you want to limit the number of files used
    # if file_num>50:
    #     break

    with nc.Dataset(file_path, 'r') as ncfile:

        latitudes = ncfile['lat'][:]
        longitudes = ncfile['lon'][:]

        longitudes, latitudes = np.meshgrid(longitudes, latitudes)




        # Accessing variables from 'co2proxy_fit_diagnostics' group
        # Data from apriori

        #THESE NEED TO BE FIGURED OUT

        example_ak = [0.54948103, 0.5463017 , 0.5425764 , 0.5630824 ,0.586732  , 0.61969376, 0.6751325 , 0.7441332 ,0.7991072 , 0.8353182 , 0.87577033, 0.92045873,0.9373637 , 0.97244173, 0.99397856, 1.0006262 , 1.0219611 , 1.0428349 , 1.0573764]
        example_prior = [ 224.05298,  685.95667, 1449.2792 , 1765.1627 ,1894.5641 , 1929.3362 , 1993.0267 , 2023.5267 ,2023.5922 , 2023.6508 , 2024.5969 , 2024.6519 ,2024.6962 , 2024.9379 , 2024.9979 , 2025.0558 ,2038.752  , 2040.5107 , 2040.6727 ]

        sp = ncfile['apriori_data']['surface_pressure']
        mtp = ncfile['co2proxy_fit_diagnostics']['ch4_averaging_kernel_mid_tropo_p']
        diff = sp[:]-mtp[:]
        tp = mtp[:]-diff
        tp = tp*0
        tp[tp==0] = 100

        print(np.nansum(tp))
        
        ptrop, c_ak = group_data_in_gridbox(tp, latitudes, longitudes, lat_edges_list, lon_edges_list)

        coverage = get_coverage(tp, latitudes, longitudes, lat_edges_list, lon_edges_list)
        
        p_trop_out = p_trop_out+ptrop


        prior_out[:, :, :] = example_prior


        ak_ch4_out[:,:,:] = example_ak

        # Access latitude and longitude data and time data
        group_name_geolocation = 'geolocation'

        time_data = ncfile['time'][:]

        xch4_co2proxy_data = ncfile['xch4'][:]
        print(xch4_co2proxy_data)
        xch4, c_ak = group_data_in_gridbox(xch4_co2proxy_data, latitudes, longitudes, lat_edges_list, lon_edges_list)
        xch4_out = xch4_out + xch4

        retrieved_surface_pressure_data = ncfile['apriori_data']['surface_pressure'][:]

        psurf, c_ak = group_data_in_gridbox(retrieved_surface_pressure_data, latitudes, longitudes, lat_edges_list, lon_edges_list)

        p_surf_out = p_surf_out+psurf

        time, c_ak = group_data_in_gridbox(time_data, latitudes, longitudes, lat_edges_list, lon_edges_list)
        time_out = time_out+time
        c_out = c_out+c_ak

    print('hello')       
    c_out[c_out==0] = np.nan
    superobs = {}
    superobs['count'] = c_out
    superobs['psurf'] = p_surf_out
    superobs['ptrop'] = p_trop_out
    superobs['ptrop'][:] = 100
    superobs['xch4'] = xch4_out
    superobs['time'] = time_out
    superobs['datetime'] = np.array([[datetime.datetime.utcfromtimestamp(float(ts)) if not np.isnan(float(ts)) else None for ts in inner] for inner in superobs['time']])


    c = np.repeat(c_out[:, :, np.newaxis], 19, axis=2)
    
    superobs['ak'] = np.flip(ak_ch4_out,2)
    superobs['prior'] = np.flip(prior_out, 2)

    superobs['pedge'] = gosat_pressure_grid(superobs['psurf'], superobs['ptrop'])
    

    out_p = compute_profile_constgrav(superobs['pedge'], g0=9.80665)
    
    gc_ll = get_sv_ll(state_vector)
    output = {}
    
    output["LATITUDE_BNDS"] = lat_edges_list
    output["LONGITUDE_BNDS"] = lon_edges_list
    
    print(gc_ll['lat'])
    
    output["LATITUDE"] = gc_ll['lat']
    output["LONGITUDE"] = gc_ll['lon']
    output["PRESSURE_EDGES"] = superobs['pedge']
    output["TIME"] = superobs['datetime']
    output["DRY_AIR_VCD"] = out_p['dry_air_vcd']
    output['PRESSURE_WEIGHT'] = out_p['dry_air_vcd']/np.tile(np.expand_dims(np.sum(out_p['dry_air_vcd'],2), axis=2), (1, 1, 19))
    output["AVERAGING_KERNEL"] = superobs['ak']
    output["PRIOR_PROFILE"] = superobs['prior']
    output["SATELLITE_COLUMN"] = superobs['xch4']
    output['OBS_COUNT'] = c_out
    print("coverage", coverage[0])
    output['OBS_ERROR'] = superobs['xch4']*0+10.0**2#((1-coverage[0])*4.7)**2+5**2
    
    #dummy
    output['COVERAGE'] = c_out
    output['COVERAGE'] = coverage
    
    valid = (
        np.isfinite(output["OBS_COUNT"])
        & (output["OBS_COUNT"] > 0)
        & np.isfinite(output["SATELLITE_COLUMN"])
    )
    j_gc, i_gc = np.where(valid)
    n_edges = output["PRESSURE_EDGES"].shape[-1]
    n_layers = output["AVERAGING_KERNEL"].shape[-1]

    dtype = [
        ("iGC", "i4"), ("jGC", "i4"),
        ("lat_sat", "f4"), ("lon_sat", "f4"),
        (species, "f4"), ("time", "U13"),
        ("p_sat", "f4", (n_edges,)),
        ("surface_pressure", "f4"),
        ("nir_albedo", "f4"), ("swir_albedo", "f4"),
        ("dry_air_subcolumns", "f4", (n_layers,)),
        ("apriori", "f4", (n_layers,)),
        ("avkern", "f4", (n_layers,)),
        ("layer", "f4", (n_layers,)),
        ("observation_count", "f4"),
        ("lat", "f4"), ("lon", "f4"),
    ]
    observations = np.zeros(len(i_gc), dtype=dtype)

    lat = np.asarray(output["LATITUDE"])
    lon = np.asarray(output["LONGITUDE"])
    observations["iGC"] = i_gc
    observations["jGC"] = j_gc
    observations["lat"] = lat[j_gc]
    observations["lon"] = lon[i_gc]
    observations["lat_sat"] = observations["lat"]
    observations["lon_sat"] = observations["lon"]
    observations[species] = output["SATELLITE_COLUMN"][valid]
    observations["p_sat"] = output["PRESSURE_EDGES"][valid]
    observations["surface_pressure"] = observations["p_sat"][:, 0]
    observations["dry_air_subcolumns"] = output["DRY_AIR_VCD"][valid]
    observations["avkern"] = output["AVERAGING_KERNEL"][valid]
    observations["layer"] = np.arange(n_layers, dtype=np.float32)
    observations["observation_count"] = output["OBS_COUNT"][valid]

    # The common IMI in-memory representation stores the prior as a partial
    # column. The canonical writer converts it back to mol/mol.
    observations["apriori"] = (
        output["PRIOR_PROFILE"][valid]
        * 1e-9
        * observations["dry_air_subcolumns"]
    )
    observations["nir_albedo"] = np.nan
    observations["swir_albedo"] = np.nan

    for index, value in enumerate(output["TIME"][valid]):
        timestamp = pd.Timestamp(value)
        observations["time"][index] = (
            get_strdate(timestamp, time_threshold)
            if time_threshold is not None
            else timestamp.round("60min").strftime("%Y%m%d_%H")
        )

    return observations


def average_methanesat_observations(
    data_path,
    state_vector,
    species="CH4",
    time_threshold=None,
    gc_startdate=None,
    gc_enddate=None,
    file_id="MSAT",
    row_chunk_size=256,
):
    """Average a rectilinear MethaneSAT L3 file onto the state-vector grid.

    The input can be much larger than memory. Only coordinate rows intersecting
    the state-vector domain are read, and the two-dimensional retrieval fields
    are accumulated in row chunks. A single joint validity mask is applied to
    XCH4, time, and surface pressure so all averaged fields and counts describe
    the same contributing pixels.
    """
    del file_id  # Retained for API compatibility.

    with xr.open_dataset(state_vector) as state_vector_ds:
        gc_lats = np.asarray(state_vector_ds["lat"].values, dtype=np.float64)
        gc_lons = np.asarray(state_vector_ds["lon"].values, dtype=np.float64)

    if gc_lats.ndim != 1 or gc_lons.ndim != 1:
        raise ValueError("MethaneSAT averaging requires 1-D state-vector lat/lon")
    if len(gc_lats) < 2 or len(gc_lons) < 2:
        raise ValueError("State-vector lat/lon must each contain at least two cells")

    def coordinate_edges(centers):
        differences = np.diff(centers)
        if not np.all(differences > 0):
            raise ValueError("State-vector coordinates must be strictly increasing")
        edges = np.empty(len(centers) + 1, dtype=np.float64)
        edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
        edges[0] = centers[0] - differences[0] / 2
        edges[-1] = centers[-1] + differences[-1] / 2
        return edges

    lat_edges = coordinate_edges(gc_lats)
    lon_edges = coordinate_edges(gc_lons)
    n_lat = len(gc_lats)
    n_lon = len(gc_lons)
    n_cells = n_lat * n_lon

    def epoch_seconds(value):
        if value is None:
            return None
        return pd.Timestamp(value).timestamp()

    start_seconds = epoch_seconds(gc_startdate)
    end_timestamp = pd.Timestamp(gc_enddate) if gc_enddate is not None else None
    end_is_date = end_timestamp is not None and end_timestamp == end_timestamp.normalize()
    end_seconds = (
        (end_timestamp + pd.Timedelta(days=1)).timestamp()
        if end_is_date else epoch_seconds(gc_enddate)
    )

    sums_xch4 = np.zeros(n_cells, dtype=np.float64)
    sums_time = np.zeros(n_cells, dtype=np.float64)
    sums_surface_pressure = np.zeros(n_cells, dtype=np.float64)
    sums_latitude = np.zeros(n_cells, dtype=np.float64)
    sums_longitude = np.zeros(n_cells, dtype=np.float64)
    counts = np.zeros(n_cells, dtype=np.int64)

    def values_with_nan(variable, key):
        values = variable[key]
        if np.ma.isMaskedArray(values):
            values = values.filled(np.nan)
        return np.asarray(values, dtype=np.float64)

    with nc.Dataset(data_path, "r") as ncfile:
        required_root = {"lat", "lon", "time", "xch4"}
        missing_root = required_root.difference(ncfile.variables)
        if missing_root:
            raise ValueError(f"MSAT file is missing variables: {sorted(missing_root)}")
        if "apriori_data" not in ncfile.groups:
            raise ValueError("MSAT file is missing the apriori_data group")
        if "surface_pressure" not in ncfile.groups["apriori_data"].variables:
            raise ValueError("MSAT file is missing apriori_data/surface_pressure")

        source_lats = np.asarray(ncfile["lat"][:], dtype=np.float64)
        source_lons = np.asarray(ncfile["lon"][:], dtype=np.float64)
        lat_selection = np.flatnonzero(
            (source_lats >= lat_edges[0]) & (source_lats <= lat_edges[-1])
        )
        lon_selection = np.flatnonzero(
            (source_lons >= lon_edges[0]) & (source_lons <= lon_edges[-1])
        )
        if len(lat_selection) == 0 or len(lon_selection) == 0:
            return _empty_methanesat_observations(species)

        # Rectilinear L3 coordinates form contiguous domain subsets.
        lat_start, lat_stop = lat_selection[0], lat_selection[-1] + 1
        lon_start, lon_stop = lon_selection[0], lon_selection[-1] + 1
        selected_lons = source_lons[lon_start:lon_stop]
        lon_bins = np.searchsorted(lon_edges, selected_lons, side="right") - 1
        valid_lon = (lon_bins >= 0) & (lon_bins < n_lon)

        for row_start in range(lat_start, lat_stop, row_chunk_size):
            row_stop = min(row_start + row_chunk_size, lat_stop)
            selected_lats = source_lats[row_start:row_stop]
            lat_bins = np.searchsorted(lat_edges, selected_lats, side="right") - 1

            key = (slice(row_start, row_stop), slice(lon_start, lon_stop))
            xch4 = values_with_nan(ncfile["xch4"], key)
            times = values_with_nan(ncfile["time"], key)
            surface_pressure = values_with_nan(
                ncfile.groups["apriori_data"]["surface_pressure"], key
            )

            valid = (
                np.isfinite(xch4) & (xch4 > 0)
                & np.isfinite(times) & (times > 0)
                & np.isfinite(surface_pressure) & (surface_pressure > 0)
            )
            if start_seconds is not None:
                valid &= times >= start_seconds
            if end_seconds is not None:
                valid &= times < end_seconds if end_is_date else times <= end_seconds
            valid &= valid_lon[np.newaxis, :]
            valid &= (lat_bins[:, np.newaxis] >= 0)
            valid &= (lat_bins[:, np.newaxis] < n_lat)
            if not np.any(valid):
                continue

            cell_indices = (
                lat_bins[:, np.newaxis] * n_lon + lon_bins[np.newaxis, :]
            )
            flat_cells = np.broadcast_to(cell_indices, valid.shape)[valid]
            counts += np.bincount(flat_cells, minlength=n_cells)
            sums_xch4 += np.bincount(
                flat_cells, weights=xch4[valid], minlength=n_cells
            )
            sums_time += np.bincount(
                flat_cells, weights=times[valid], minlength=n_cells
            )
            sums_surface_pressure += np.bincount(
                flat_cells, weights=surface_pressure[valid], minlength=n_cells
            )
            chunk_latitudes = np.broadcast_to(
                selected_lats[:, np.newaxis], valid.shape
            )
            chunk_longitudes = np.broadcast_to(
                selected_lons[np.newaxis, :], valid.shape
            )
            sums_latitude += np.bincount(
                flat_cells, weights=chunk_latitudes[valid], minlength=n_cells
            )
            sums_longitude += np.bincount(
                flat_cells, weights=chunk_longitudes[valid], minlength=n_cells
            )

    populated = counts > 0
    if not np.any(populated):
        return _empty_methanesat_observations(species)

    mean_xch4 = sums_xch4[populated] / counts[populated]
    mean_time = sums_time[populated] / counts[populated]
    mean_surface_pressure = sums_surface_pressure[populated] / counts[populated]
    mean_latitude = sums_latitude[populated] / counts[populated]
    mean_longitude = sums_longitude[populated] / counts[populated]
    j_gc, i_gc = np.unravel_index(np.flatnonzero(populated), (n_lat, n_lon))

    # Fixed profiles retained from the original MethaneSAT implementation.
    fixed_ak = np.array([
        0.54948103, 0.5463017, 0.5425764, 0.5630824, 0.586732,
        0.61969376, 0.6751325, 0.7441332, 0.7991072, 0.8353182,
        0.87577033, 0.92045873, 0.9373637, 0.97244173, 0.99397856,
        1.0006262, 1.0219611, 1.0428349, 1.0573764,
    ], dtype=np.float64)[::-1]
    fixed_prior = np.array([
        224.05298, 685.95667, 1449.2792, 1765.1627, 1894.5641,
        1929.3362, 1993.0267, 2023.5267, 2023.5922, 2023.6508,
        2024.5969, 2024.6519, 2024.6962, 2024.9379, 2024.9979,
        2025.0558, 2038.752, 2040.5107, 2040.6727,
    ], dtype=np.float64)[::-1]

    tropopause_pressure = np.full_like(mean_surface_pressure, 100.0)
    pressure_edges = gosat_pressure_grid(
        mean_surface_pressure[:, np.newaxis],
        tropopause_pressure[:, np.newaxis],
    )[:, 0, :]
    dry_air = compute_profile_constgrav(
        pressure_edges[:, np.newaxis, :]
    )["dry_air_vcd"][:, 0, :]

    observations = _empty_methanesat_observations(species, len(i_gc))
    observations["iGC"] = i_gc
    observations["jGC"] = j_gc
    observations["lat"] = gc_lats[j_gc]
    observations["lon"] = gc_lons[i_gc]
    observations["lat_sat"] = mean_latitude
    observations["lon_sat"] = mean_longitude
    observations[species] = mean_xch4
    observations["p_sat"] = pressure_edges
    observations["surface_pressure"] = mean_surface_pressure
    observations["dry_air_subcolumns"] = dry_air
    observations["apriori"] = fixed_prior * 1e-9 * dry_air
    observations["avkern"] = fixed_ak
    observations["layer"] = np.arange(19, dtype=np.float32)
    observations["observation_count"] = counts[populated]
    observations["nir_albedo"] = np.nan
    observations["swir_albedo"] = np.nan

    for index, seconds in enumerate(mean_time):
        timestamp = pd.to_datetime(seconds, unit="s", utc=True).tz_localize(None)
        observations["time"][index] = (
            get_strdate(timestamp, time_threshold)
            if time_threshold is not None
            else timestamp.round("60min").strftime("%Y%m%d_%H")
        )
    return observations


def _empty_methanesat_observations(species, size=0):
    """Allocate the standard IMI representation for 19-layer MSAT data."""
    dtype = [
        ("iGC", "i4"), ("jGC", "i4"),
        ("lat_sat", "f4"), ("lon_sat", "f4"),
        (species, "f4"), ("time", "U13"),
        ("p_sat", "f4", (20,)), ("surface_pressure", "f4"),
        ("nir_albedo", "f4"), ("swir_albedo", "f4"),
        ("dry_air_subcolumns", "f4", (19,)),
        ("apriori", "f4", (19,)), ("avkern", "f4", (19,)),
        ("layer", "f4", (19,)), ("observation_count", "f4"),
        ("lat", "f4"), ("lon", "f4"),
    ]
    return np.zeros(size, dtype=dtype)


def gosat_pressure_grid(psurf,ptrop):

    ''' UoL GOSAT Pressure grid 

      ARGS:
        psurf[x,t]: Surface pressure [hPa]
        ptrop[x,t]: Tropopause pressure [hPa]

      RETURNS:
        pedge[x,t,z]: UoL Proxy GOSAT Pressure grid

    '''

    # Allocate output pressure edges
    lmx = 19 ; lmx_t = 13
    xmx = psurf.shape[0]
    tmx = psurf.shape[1]
    pedge = np.zeros((xmx, tmx, lmx+1))

    # Tropospheric sigma levels
    tsig = np.linspace(0,1,lmx_t+1)[::-1]
    for l in range(lmx_t+1):
        pedge[:,:,l] = ptrop + tsig[l]*(psurf-ptrop)

    # Constant stratospheric levels
    p_strat = [8.0000000e1, 5.0000000e1, 1.0000000e1, 1.0000000, 1.0000000e-1]
    for l,pl in zip(range(lmx_t+2,lmx+1),p_strat):
        pedge[:,:,l] = pl

    # Intermediate strat level
    pedge[:,:,lmx_t+1] = 0.5*(pedge[:,:,lmx_t]+pedge[:, :, lmx_t+2])

    return pedge

def compute_profile_constgrav(pedge, g0=9.80665):

    # molecular weights [kg/mol]
    mw_dry = 28.9647e-3  # dry air
    mw_h2o = 18.01528e-3 # water vapor

    outp = {}
    outp['pedge'] = pedge
    zmx = outp['pedge'].shape[-1]-1
    outp['pmid'] = 0.5*(outp['pedge'][:, :, 0:zmx]+outp['pedge'][:, :, 1:])

    outp['air_vcd'] = np.zeros_like(outp['pmid'])
    for iz in range(zmx):
#        outp['air_vcd'][:,:,iz] = (outp['pedge'][:,:,iz]-outp['pedge'][:,:,iz+1])/g0 * (1/(outp['mw'][:,:,iz]))
        outp['air_vcd'][:,:,iz] = (outp['pedge'][:,:,iz]-outp['pedge'][:,:,iz+1])*100/g0 * (1/mw_dry)
    outp['dry_air_vcd'] = outp['air_vcd']#*(1.0-xh2o)

    return outp

def group_data_in_gridbox(data, lat2d, lon2d, lat_bnds, lon_bnds, fill_value=1e30):
    """
    Bin and average 2D data into lat/lon gridboxes, handling _FillValue properly.

    Parameters:
        data       : 2D array (nlat, nlon)
        lat2d      : 2D array (nlat, nlon)
        lon2d      : 2D array (nlat, nlon)
        lat_bnds   : 1D array of latitude bin edges
        lon_bnds   : 1D array of longitude bin edges
        fill_value : Value representing missing data (default 1e36)

    Returns:
        avg_data : 2D array of mean values per bin
        count    : 2D array of counts per bin
    """
    # Flatten arrays
    flat_data = data.flatten()
    flat_lat = lat2d.flatten()
    flat_lon = lon2d.flatten()

    # Convert longitudes to -180–180 if needed
    flat_lon = np.where(flat_lon > 180, flat_lon - 360, flat_lon)

    # Mask out fill values and NaNs
    valid = (flat_data < (fill_value+1)) & np.isfinite(flat_data) & (flat_data > 0)
    flat_data = flat_data[valid]
    flat_lat = flat_lat[valid]
    flat_lon = flat_lon[valid]

    # Bin averages
    avg_data, _, _, _ = binned_statistic_2d(
        flat_lat, flat_lon, flat_data,
        statistic='mean',
        bins=[lat_bnds, lon_bnds]
    )

    # Bin counts
    count, _, _, _ = binned_statistic_2d(
        flat_lat, flat_lon, np.ones_like(flat_data),
        statistic='count',
        bins=[lat_bnds, lon_bnds]
    )

    return avg_data, count



def get_coverage(data, lat2d, lon2d, lat_bnds, lon_bnds, fill_value=1e30):
    """
    Bin and average 2D data into lat/lon gridboxes, handling _FillValue properly.

    Parameters:
        data       : 2D array (nlat, nlon)
        lat2d      : 2D array (nlat, nlon)
        lon2d      : 2D array (nlat, nlon)
        lat_bnds   : 1D array of latitude bin edges
        lon_bnds   : 1D array of longitude bin edges
        fill_value : Value representing missing data (default 1e36)

    Returns:
        avg_data : 2D array of mean values per bin
        count    : 2D array of counts per bin
    """
    # Flatten arrays
    flat_data = data.flatten()
    flat_lat = lat2d.flatten()
    flat_lon = lon2d.flatten()

    # Convert longitudes to -180–180 if needed
    flat_lon = np.where(flat_lon > 180, flat_lon - 360, flat_lon)

    count_total, _, _, _ = binned_statistic_2d(
        flat_lat, flat_lon, np.ones_like(flat_data),
        statistic='count',
        bins=[lat_bnds, lon_bnds]
    )

    # Mask out fill values and NaNs
    valid = (flat_data < (fill_value+1)) & np.isfinite(flat_data) & (flat_data > 0)
    flat_data = flat_data[valid]
    flat_lat = flat_lat[valid]
    flat_lon = flat_lon[valid]


    count_masked, _, _, _ = binned_statistic_2d(
        flat_lat, flat_lon, np.ones_like(flat_data),
        statistic='count',
        bins=[lat_bnds, lon_bnds]
    )

    coverage = count_masked/count_total

    return coverage, valid

    
