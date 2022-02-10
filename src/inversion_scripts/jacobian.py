#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
from multiprocessing.dummy import connection
import numpy as np
import xarray as xr
import re
import pickle
import os
import pandas as pd 
import datetime
from shapely.geometry import Polygon
from utils import save_obj

    
def read_tropomi(filename):
    """
    Read TROPOMI data and save important variables to dictionary.

    Arguments
        filename [str]  : TROPOMI netcdf data file to read

    Returns
        dat      [dict] : Dictionary of important variables from TROPOMI:
                            - CH4
                            - Latitude
                            - Longitude
                            - QA value
                            - UTC time
                            - Time (utc time reshaped for orbit)
                            - Averaging kernel
                            - SWIR albedo
                            - NIR albedo
                            - Blended albedo
                            - CH4 prior profile
                            - Dry air subcolumns
                            - Latitude bounds
                            - Longitude bounds
                            - Vertical pressure profile
    """

    # Initialize dictionary for TROPOMI data
    dat = {}
    
    # Store methane, QA, lat, lon
    tropomi_data = xr.open_dataset(filename, group='PRODUCT')
    dat['methane'] = tropomi_data['methane_mixing_ratio_bias_corrected'].values[0,:,:]
    dat['qa_value'] = tropomi_data['qa_value'].values[0,:,:]
    dat['longitude'] = tropomi_data['longitude'].values[0,:,:]
    dat['latitude'] = tropomi_data['latitude'].values[0,:,:]
    
    # Store UTC time
    reference_time = tropomi_data['time'].values
    delta_time = tropomi_data['delta_time'][0].values
    strdate = []
    if delta_time.dtype == '<m8[ns]':
        strdate = reference_time+delta_time
    elif delta_time.dtype == '<M8[ns]':
        strdate = delta_time
    else:
        print(delta_time.dtype)
        pass
    dat['utctime'] = strdate
    
    # Store time for whole orbit
    times = np.zeros(shape=dat['longitude'].shape, dtype='datetime64[ns]')
    for k in range(dat['longitude'].shape[0]):
        times[k,:] = strdate[k]
    dat['time'] = times
    tropomi_data.close()
 
    # Store column averaging kernel, SWIR and NIR surface albedo
    tropomi_data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
    dat['column_AK'] = tropomi_data['column_averaging_kernel'].values[0,:,:,::-1]
    dat['swir_albedo'] = tropomi_data['surface_albedo_SWIR'].values[0,:,:]
    dat['nir_albedo'] = tropomi_data['surface_albedo_NIR'].values[0,:,:]
    dat['blended_albedo'] = 2.4*dat['nir_albedo'] - 1.13*dat['swir_albedo']
    tropomi_data.close()
    
    # Store methane prior profile, dry air subcolumns
    tropomi_data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
    dat['methane_profile_apriori'] = tropomi_data['methane_profile_apriori'].values[0,:,:,::-1] # mol m-2
    dat['dry_air_subcolumns'] = tropomi_data['dry_air_subcolumns'].values[0,:,:,::-1]           # mol m-2
    
    # Also get pressure interval and surface pressure for use below
    pressure_interval = tropomi_data['pressure_interval'].values[0,:,:]/100                     # Pa -> hPa
    surface_pressure = tropomi_data['surface_pressure'].values[0,:,:]/100                       # Pa -> hPa
    tropomi_data.close()

    # Store latitude and longitude bounds for pixels
    tropomi_data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
    dat['longitude_bounds'] = tropomi_data['longitude_bounds'].values[0,:,:,:]
    dat['latitude_bounds'] = tropomi_data['latitude_bounds'].values[0,:,:,:]
    tropomi_data.close()
    
    # Store vertical pressure profile
    n1 = dat['methane'].shape[0] # length of along-track dimension (scanline) of retrieval field
    n2 = dat['methane'].shape[1] # length of across-track dimension (ground_pixel) of retrieval field
    pressures = np.zeros([n1,n2,12+1], dtype=np.float32)
    pressures.fill(np.nan)
    for i in range(12+1):
        pressures[:,:,i] = surface_pressure - i*pressure_interval
    dat['pressures'] = pressures
    
    return dat


def read_geoschem(date, gc_cache, build_jacobian=False, sens_cache=None):
    """
    Read GEOS-Chem data and save important variables to dictionary.

    Arguments
        date           [str]   : Date of interest
        gc_cache       [str]   : Path to GEOS-Chem output data
        build_jacobian [log]   : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        sens_cache     [str]   : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data
    
    Returns
        dat            [dict]  : Dictionary of important variables from GEOS-Chem:
                                    - CH4
                                    - Latitude
                                    - Longitude
                                    - PEDGE
    """
    
    # Assemble file paths to GEOS-Chem output collections for input data
    file_species = f'GEOSChem.SpeciesConc.{date}00z.nc4'        
    file_pedge = f'GEOSChem.LevelEdgeDiags.{date}00z.nc4'    
    
    # Read lat, lon, CH4 from the SpeciecConc collection
    filename = f'{gc_cache}/{file_species}'
    gc_data = xr.open_dataset(filename)
    LON = gc_data['lon'].values
    LAT = gc_data['lat'].values
    CH4 = gc_data['SpeciesConc_CH4'].values[0,:,:,:]
    CH4 = CH4*1e9 # Convert to ppb
    CH4 = np.einsum('lij->jil', CH4)
    gc_data.close()

    # Read PEDGE from the LevelEdgeDiags collection
    filename = f'{gc_cache}/{file_pedge}'
    gc_data = xr.open_dataset(filename)
    PEDGE = gc_data['Met_PEDGE'].values[0,:,:,:]
    PEDGE = np.einsum('lij->jil', PEDGE)
    gc_data.close()
    
    # Store GEOS-Chem data in dictionary
    dat = {}
    dat['lon'] = LON
    dat['lat'] = LAT
    dat['PEDGE'] = PEDGE
    dat['CH4'] = CH4
    
    # If need to construct Jacobian, read sensitivity data from GEOS-Chem perturbation simulations
    if build_jacobian:
        sensitivity_data = xr.open_dataset(f'{sens_cache}/Sensi_{date}.nc')
        sensitivities = sensitivity_data['Sensitivities'].values
        # Reshape so the data have dimensions (lon, lat, lev, grid_element)
        sensitivities = np.einsum('klji->ijlk', sensitivities)
        sensitivity_data.close()
        dat['Sensitivities'] = sensitivities

    return dat


def read_all_geoschem(all_strdate, gc_cache, build_jacobian=False, sens_cache=None):
    """ 
    Call readgeoschem() for multiple dates in a loop. 

    Arguments
        all_strdate    [list, str] : Multiple date strings
        gc_cache       [str]       : Path to GEOS-Chem output data
        build_jacobian [log]       : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        sens_cache     [str]       : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

    Returns
        dat            [dict]      : Dictionary of dictionaries. Each sub-dictionary is returned by read_geoschem()
    """

    dat={}
    for strdate in all_strdate:
        dat[strdate] = read_geoschem(strdate, gc_cache, build_jacobian, sens_cache) 
    
    return dat


def merge_pressure_grids(p_sat, p_gc):
    """
    Merge TROPOMI & GEOS-Chem vertical pressure grids

    Arguments
        p_sat   [float]    : Pressure edges from TROPOMI (13 edges)     <--- 13-1 = 12 pressure layers
        p_gc    [float]    : Pressure edges from GEOS-Chem (48 edges)   <--- 48-1 = 47 pressure layers

    Returns
        merged  [dict]     : Merged grid dictionary
                                - p_merge       : merged pressure-edge grid 
                                - data_type     : for each pressure edge in the merged grid, is it from GEOS-Chem or TROPOMI?
                                - edge_index    : indexes of pressure edges
                                - first_gc_edge : index of first GEOS-Chem pressure edge in the merged grid
    """

    # Combine p_sat and p_gc into merged vertical pressure grid
    p_merge = np.zeros(len(p_sat)+len(p_gc))
    p_merge.fill(np.nan)
    data_type = np.zeros(len(p_sat)+len(p_gc), dtype=int) 
    data_type.fill(-99)
    edge_index = []
    i=0; j=0; k=0
    while ((i < len(p_sat)) or (j < len(p_gc))):
        if i == len(p_sat):
            p_merge[k] = p_gc[j]
            data_type[k] = 2      # geos-chem edge
            j=j+1; k=k+1
            continue
        if j == len(p_gc):
            p_merge[k] = p_sat[i]
            data_type[k] = 1      # tropomi edge
            edge_index.append(k)        
            i=i+1; k=k+1   
            continue
        if p_sat[i] >= p_gc[j]:
            p_merge[k] = p_sat[i]
            data_type[k] = 1      # tropomi edge   
            edge_index.append(k)        
            i=i+1; k=k+1
        else:
            p_merge[k] = p_gc[j]
            data_type[k] = 2      # geos-chem edge
            j=j+1; k=k+1
    
    # Find the first GEOS-Chem pressure edge
    first_gc_edge = -99
    for i in range(len(p_sat)+len(p_gc)-1):
        if data_type[i] == 2:
            first_gc_edge = i
            break    
    
    # Save data to dictionary
    merged = {}
    merged['p_merge'] = p_merge
    merged['data_type'] = data_type
    merged['edge_index'] = edge_index
    merged['first_gc_edge'] = first_gc_edge
    
    return merged


def remap(gc_CH4, data_type, p_merge, edge_index, first_gc_edge):
    """
    Remap GEOS-Chem methane to the TROPOMI vertical grid.

    Arguments
        gc_CH4        [float]   : Methane from GEOS-Chem
        p_merge       [float]   : Merged TROPOMI + GEOS-Chem pressure levels, from merge_pressure_grids()
        data_type     [int]     : Labels for pressure edges of merged grid. 1=TROPOMI, 2=GEOS-Chem, from merge_pressure_grids()
        edge_index    [int]     : Indexes of pressure edges, from merge_pressure_grids()
        first_gc_edge [int]     : Index of first GEOS-Chem pressure edge in merged grid, from merge_pressure_grids()

    Returns
        sat_CH4       [float]   : GEOS-Chem methane in TROPOMI pressure coordinates
    """
    
    # Define CH4 concentrations in the layers of the merged pressure grid
    CH4 = np.zeros(len(p_merge)-1,)
    CH4.fill(np.nan)
    k=0
    for i in range(first_gc_edge, len(p_merge)-1):
        CH4[i] = gc_CH4[k]
        if data_type[i+1] == 2:
            k=k+1
    if first_gc_edge > 0:
        CH4[:first_gc_edge] = CH4[first_gc_edge]
    
    # Calculate the pressure-weighted mean methane for each TROPOMI layer
    delta_p = p_merge[:-1] - p_merge[1:]
    sat_CH4 = np.zeros(12)
    sat_CH4.fill(np.nan)
    for i in range(len(edge_index)-1):
        start = edge_index[i]
        end = edge_index[i+1]
        sum_conc_times_pressure_weights = sum(CH4[start:end] * delta_p[start:end])
        sum_pressure_weights = sum(delta_p[start:end])
        sat_CH4[i] = sum_conc_times_pressure_weights/sum_pressure_weights # pressure-weighted average
    
    return sat_CH4


def remap_sensitivities(sens_lonlat, data_type, p_merge, edge_index, first_gc_edge):
    """
    Remap GEOS-Chem sensitivity data (from perturbation simulations) to the TROPOMI vertical grid.

    Arguments
        sens_lonlat   [float]   : Sensitivity data from GEOS-Chem perturbation runs, for a specific lon/lat; has dims (lev, grid_element)
        p_merge       [float]   : Combined TROPOMI + GEOS-Chem pressure levels, from merge_pressure_grids()
        data_type     [int]     : Labels for pressure edges of merged grid. 1=TROPOMI, 2=GEOS-Chem, from merge_pressure_grids()
        edge_index    [int]     : Indexes of pressure edges, from merge_pressure_grids()
        first_gc_edge [int]     : Index of first GEOS-Chem pressure edge in merged grid, from merge_pressure_grids()

    Returns
        sat_deltaCH4  [float]   : GEOS-Chem methane sensitivities in TROPOMI pressure coordinates
    """
    
    # Define DeltaCH4 in the layers of the merged pressure grid, for all perturbed state vector elements
    n_elem = sens_lonlat.shape[1]
    deltaCH4 = np.zeros((len(p_merge)-1, n_elem))
    deltaCH4.fill(np.nan)
    k=0
    for i in range(first_gc_edge, len(p_merge)-1):
        deltaCH4[i,:] = sens_lonlat[k,:]
        if data_type[i+1] == 2:
            k=k+1
    if first_gc_edge > 0:
        deltaCH4[:first_gc_edge,:] = deltaCH4[first_gc_edge,:]
    
    # Calculate the weighted mean DeltaCH4 for each layer, for all perturbed state vector elements
    delta_p = p_merge[:-1] - p_merge[1:]
    delta_ps = np.transpose(np.tile(delta_p, (n_elem,1)))
    sat_deltaCH4 = np.zeros((12, n_elem))
    sat_deltaCH4.fill(np.nan)
    for i in range(len(edge_index)-1):
        start = edge_index[i]
        end = edge_index[i+1]
        sum_conc_times_pressure_weights = np.sum(deltaCH4[start:end,:] * delta_ps[start:end,:], 0)
        sum_pressure_weights = np.sum(delta_p[start:end])
        sat_deltaCH4[i,:] = sum_conc_times_pressure_weights/sum_pressure_weights # pressure-weighted average
    
    return sat_deltaCH4


def nearest_loc(query_location, reference_grid, tolerance=0.5):
    """ Find the index of the nearest grid location to a query location, with some tolerance. """

    distances = np.abs(reference_grid - query_location)
    ind = distances.argmin()
    if distances[ind] >= tolerance:
        return np.nan
    else:
        return ind


def apply_tropomi_operator(filename, n_elements, gc_startdate, gc_enddate, xlim, ylim, gc_cache, build_jacobian, sens_cache):
    """
    Apply the tropomi operator to map GEOS-Chem methane data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        n_elements     [int]        : Number of state vector elements
        gc_startdate   [datetime64] : First day of inversion period, for GEOS-Chem and TROPOMI
        gc_enddate     [datetime64] : Last day of inversion period, for GEOS-Chem and TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        gc_cache       [str]        : Path to GEOS-Chem output data
        build_jacobian [log]        : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        sens_cache     [str]        : If build_jacobian=True, this is the path to the GEOS-Chem sensitivity data

    Returns
        output         [dict]       : Dictionary with one or two fields:
 	   	 		                        - obs_GC : GEOS-Chem and TROPOMI methane data
                                                    - TROPOMI methane
                                                    - GEOS-Chem methane
                                                    - TROPOMI lat, lon
                                                    - TROPOMI lat index, lon index
				                      If build_jacobian=True, also include:
				                        - K      : Jacobian matrix
    """
    
    # Read TROPOMI data
    TROPOMI = read_tropomi(filename)
    
    # We're only going to consider data within lat/lon/time bounds, with QA > 0.5, and with safe surface albedo values
    sat_ind = np.where((TROPOMI['longitude'] >  xlim[0])      & (TROPOMI['longitude'] <  xlim[1])     & 
                       (TROPOMI['latitude']  >  ylim[0])      & (TROPOMI['latitude']  <  ylim[1])     & 
                       (TROPOMI['time'] >= gc_startdate)      & (TROPOMI['time'] <= gc_enddate)       &
                       (TROPOMI['qa_value']  >= 0.5)          &
                       (TROPOMI['swir_albedo'] > 0.05)        & (TROPOMI['blended_albedo'] < 0.85))

    # Number of TROPOMI observations
    n_obs = len(sat_ind[0])
    print('Found', n_obs, 'TROPOMI observations.')

    # If need to build Jacobian from GEOS-Chem perturbation simulation sensitivity data:
    if build_jacobian:
        # Initialize Jacobian K
        jacobian_K = np.zeros([n_obs,n_elements], dtype=np.float32)
        jacobian_K.fill(np.nan)
    
    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # For each TROPOMI observation
    for k in range(n_obs):  
        # Get the date and hour
        iSat = sat_ind[0][k] # lat index
        jSat = sat_ind[1][k] # lon index
        time = pd.to_datetime(str(TROPOMI['utctime'][iSat]))
        strdate = time.round('60min').strftime('%Y%m%d_%H')        
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_gc = read_all_geoschem(all_strdate, gc_cache, build_jacobian, sens_cache)
    
    # Initialize array with n_obs rows and 6 columns. Columns are TROPOMI CH4, GEOSChem CH4, longitude, latitude, II, JJ
    obs_GC = np.zeros([n_obs,6], dtype=np.float32)
    obs_GC.fill(np.nan)

    # For each TROPOMI observation: 
    for k in range(n_obs):
        
        # Get GEOS-Chem data for the date of the observation:
        iSat = sat_ind[0][k]
        jSat = sat_ind[1][k]
        p_sat = TROPOMI['pressures'][iSat,jSat,:]
        dry_air_subcolumns = TROPOMI['dry_air_subcolumns'][iSat,jSat,:] # mol m-2
        apriori = TROPOMI['methane_profile_apriori'][iSat,jSat,:]       # mol m-2
        avkern = TROPOMI['column_AK'][iSat,jSat,:]
        time = pd.to_datetime(str(TROPOMI['utctime'][iSat]))
        strdate = time.round('60min').strftime('%Y%m%d_%H')        
        GEOSCHEM = all_date_gc[strdate]
                
        # Find GEOS-Chem lats & lons closest to the corners of the TROPOMI pixel
        longitude_bounds = TROPOMI['longitude_bounds'][iSat,jSat,:]
        latitude_bounds = TROPOMI['latitude_bounds'][iSat,jSat,:]
        corners_lon_index = []
        corners_lat_index = []
        for l in range(4):
            iGC = nearest_loc(longitude_bounds[l], GEOSCHEM['lon'])
            jGC = nearest_loc(latitude_bounds[l], GEOSCHEM['lat'])
            corners_lon_index.append(iGC)
            corners_lat_index.append(jGC)
        # If the tolerance in nearest_loc() is not satisfied, skip the observation 
        if np.nan in corners_lon_index+corners_lat_index:
            continue
        # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the TROPOMI corners
        ij_GC = [(x,y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
        gc_coords = [(GEOSCHEM['lon'][i], GEOSCHEM['lat'][j]) for i,j in ij_GC]
        
        # Compute the overlapping area between the TROPOMI pixel and GEOS-Chem grid cells it touches
        overlap_area = np.zeros(len(gc_coords))
        dlon = GEOSCHEM['lon'][1] - GEOSCHEM['lon'][0]
        dlat = GEOSCHEM['lat'][1] - GEOSCHEM['lat'][0]
        # Polygon representing TROPOMI pixel
        polygon_tropomi = Polygon(np.column_stack((longitude_bounds,latitude_bounds)))
        # For each GEOS-Chem grid cell that touches the TROPOMI pixel: 
        for gridcellIndex in range(len(gc_coords)):
            # Define polygon representing the GEOS-Chem grid cell
            coords = gc_coords[gridcellIndex]
            geoschem_corners_lon = [coords[0]-dlon/2, coords[0]+dlon/2, coords[0]+dlon/2, coords[0]-dlon/2]
            geoschem_corners_lat = [coords[1]-dlat/2, coords[1]-dlat/2, coords[1]+dlat/2, coords[1]+dlat/2]
            polygon_geoschem = Polygon(np.column_stack((geoschem_corners_lon, geoschem_corners_lat)))
            # Calculate overlapping area as the intersection of the two polygons
            if polygon_geoschem.intersects(polygon_tropomi):
                  overlap_area[gridcellIndex] = polygon_tropomi.intersection(polygon_geoschem).area
        
        # If there is no overlap between GEOS-Chem and TROPOMI, skip to next observation:
        if sum(overlap_area) == 0:
            continue                  

        # =======================================================
        #       Map GEOS-Chem to TROPOMI observation space
        # =======================================================
        
        # Otherwise, initialize tropomi virtual xch4 and virtual sensitivity as zero
        area_weighted_virtual_tropomi = 0             # virtual tropomi xch4
        area_weighted_virtual_tropomi_sensitivity = 0 # virtual tropomi sensitivity
        
        # For each GEOS-Chem grid cell that touches the TROPOMI pixel: 
        for gridcellIndex in range(len(gc_coords)):
            
            # Get GEOS-Chem lat/lon indices for the cell
            iGC,jGC = ij_GC[gridcellIndex]
            
            # Get GEOS-Chem pressure edges for the cell
            p_gc = GEOSCHEM['PEDGE'][iGC,jGC,:]
            
            # Get GEOS-Chem methane for the cell
            gc_CH4 = GEOSCHEM['CH4'][iGC,jGC,:]
                
            # Get merged GEOS-Chem/TROPOMI pressure grid for the cell
            merged = merge_pressure_grids(p_sat, p_gc)
            
            # Remap GEOS-Chem methane to TROPOMI pressure levels
            sat_CH4 = remap(gc_CH4, merged['data_type'], merged['p_merge'], merged['edge_index'], merged['first_gc_edge']) # ppb
            
            # Convert ppb to mol m-2
            sat_CH4_molm2 = sat_CH4 * 1e-9 * dry_air_subcolumns   # mol m-2
            
            # Derive the column-averaged XCH4 that TROPOMI would see over this ground cell
            # using eq. 46 from TROPOMI Methane ATBD, Hasekamp et al. 2019
            virtual_tropomi_gridcellIndex = sum(apriori + avkern * (sat_CH4_molm2 - apriori)) / sum(dry_air_subcolumns) * 1e9 # ppb
            
            # Weight by overlapping area (to be divided out later) and add to sum
            area_weighted_virtual_tropomi += overlap_area[gridcellIndex] * virtual_tropomi_gridcellIndex # ppb m2

            # If building Jacobian matrix from GEOS-Chem perturbation simulation sensitivity data:
            if build_jacobian:
                
                # Get GEOS-Chem perturbation sensitivities at this lat/lon, for all vertical levels and state vector elements
                sens_lonlat = GEOSCHEM['Sensitivities'][iGC,jGC,:,:]
                
                # Map the sensitivities to TROPOMI pressure levels
                sat_deltaCH4 = remap_sensitivities(sens_lonlat, merged['data_type'], merged['p_merge'], merged['edge_index'], merged['first_gc_edge']) # mixing ratio, unitless
                
                # Tile the TROPOMI averaging kernel
                avkern_tiled = np.transpose(np.tile(avkern, (n_elements,1)))
                
                # Tile the TROPOMI dry air subcolumns
                dry_air_subcolumns_tiled = np.transpose(np.tile(dry_air_subcolumns, (n_elements,1)))   # mol m-2
                
                # Derive the change in column-averaged XCH4 that TROPOMI would see over this ground cell
                tropomi_sensitivity_gridcellIndex = np.sum(avkern_tiled*sat_deltaCH4*dry_air_subcolumns_tiled, 0) / sum(dry_air_subcolumns) # mixing ratio, unitless
                
                # Weight by overlapping area (to be divided out later) and add to sum
                area_weighted_virtual_tropomi_sensitivity += overlap_area[gridcellIndex] * tropomi_sensitivity_gridcellIndex # m2

        # Compute virtual TROPOMI observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_tropomi = area_weighted_virtual_tropomi / sum(overlap_area)
 
        # Save actual and virtual TROPOMI data                    
        obs_GC[k,0] = TROPOMI['methane'][iSat,jSat]     # Actual TROPOMI methane column observation
        obs_GC[k,1] = virtual_tropomi                   # Virtual TROPOMI methane column observation
        obs_GC[k,2] = TROPOMI['longitude'][iSat,jSat]   # TROPOMI longitude
        obs_GC[k,3] = TROPOMI['latitude'][iSat,jSat]    # TROPOMI latitude 
        obs_GC[k,4] = iSat                              # TROPOMI index of longitude
        obs_GC[k,5] = jSat                              # TROPOMI index of latitude

        if build_jacobian:
            # Compute TROPOMI sensitivity as weighted mean by overlapping area
            # i.e., need to divide out area [m2] from the previous step
            jacobian_K[k,:] = area_weighted_virtual_tropomi_sensitivity / sum(overlap_area)
    
    # Output 
    output = {}

    # Always return the coincident TROPOMI and GEOS-Chem data
    output['obs_GC'] = obs_GC

    # Optionally return the Jacobian
    if build_jacobian:
        output['K'] = jacobian_K
    
    return output




if __name__ == '__main__':
    import sys

    startday = sys.argv[1]
    endday = sys.argv[2]
    lonmin = float(sys.argv[3])
    lonmax = float(sys.argv[4])
    latmin = float(sys.argv[5])
    latmax = float(sys.argv[6])
    n_elements = int(sys.argv[7])
    tropomi_cache = sys.argv[8]
    isPost = sys.argv[9]
 
    # Reformat start and end days for datetime in configuration
    start = f'{startday[0:4]}-{startday[4:6]}-{startday[6:8]} 00:00:00'
    end = f'{endday[0:4]}-{endday[4:6]}-{endday[6:8]} 23:59:59'

    # Configuration
    workdir = '.'
    sens_cache = f'{workdir}/Sensi'
    if isPost.lower() == 'false':
        build_jacobian = True
        gc_cache = f'{workdir}/data_GC'
        outputdir = f'{workdir}/data_converted'
    else: # if sampling posterior simulation
        build_jacobian = False
        gc_cache = f'{workdir}/data_GC_posterior'
        outputdir = f'{workdir}/data_converted_posterior'
    xlim = [lonmin,lonmax]
    ylim = [latmin,latmax]
    gc_startdate = np.datetime64(datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S'))
    gc_enddate = np.datetime64(datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(days=1))
    print('Start:', start)
    print('End:', end)
    
    # Get TROPOMI data filenames for the desired date range
    allfiles = glob.glob(f'{tropomi_cache}/*.nc')
    sat_files = []
    for index in range(len(allfiles)):
        filename = allfiles[index]
        shortname = re.split('\/', filename)[-1]
        shortname = re.split('\.', shortname)[0]
        strdate = re.split('\.|_+|T',shortname)[4]
        strdate = datetime.datetime.strptime(strdate, '%Y%m%d')
        if ((strdate >= gc_startdate) and (strdate <= gc_enddate)):
            sat_files.append(filename)
    sat_files.sort()
    print('Found', len(sat_files), 'TROPOMI data files.')

    # Map GEOS-Chem to TROPOMI observation space
    # Also return Jacobian matrix if build_jacobian=True
    for filename in sat_files:
        
        # Check if TROPOMI file has already been processed
        print('========================')
        shortname = re.split('\/', filename)[-1]
        print(shortname)
        date = re.split('\.',shortname)[0]
        
        # If not yet processed, run apply_tropomi_operator()
        if not os.path.isfile(f'{outputdir}/{date}_GCtoTROPOMI.pkl'):
            print('Applying TROPOMI operator...')
            output = apply_tropomi_operator(filename, n_elements, gc_startdate, gc_enddate, xlim, ylim, gc_cache, build_jacobian, sens_cache)
        if output['obs_GC'].shape[0] > 0:
            print('Saving .pkl file')
            save_obj(output, f'{outputdir}/{date}_GCtoTROPOMI.pkl')

    print(f'Wrote files to {outputdir}')
