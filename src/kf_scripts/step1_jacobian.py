#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import numpy as np
import xarray as xr
import re
import pickle
import os
import pandas as pd 
import datetime
from shapely.geometry import Polygon

# Notes:
# ======
# - Manual step 1: Pre-define the GEOS-Chem sensitivities from the perturbation simulations.
#   Once all perturbation runs are finished, subtract output of the base run from all of the 
#   perturbation outputs, to obtain one 4D sensitivity array for each simulation time step 
#   (one 4D array per hour, in Lu's setup).
# - Manual step 2: Build the lat_ratio.csv file used to define the lat_ratio variable.
# - Not sure why we need both local times and UTC times.
# - Not clear on how longitudes get transformed into time offsets for computing local times.
# - Lu optimizes scaling factors, so he increases the sensitivities by a factor of two, due to 
#   1.5 scaling perturbation. This doesn't work if we are optimizing absolute emissions.
# - Not sure about some variables in the cal_weights() and remap()/remap2() functions:
# 	- data_type
# 	- location
#	- first_2
# - When Lu computes virtual TROPOMI column from GC data using the TROPOMI prior and averaging
#   kernel, he does it as the weighted mean mixing ratio [ppb] of the relevant GC ground cells.
#   Zhen does it as the weighted mean of number of molecules instead. This requires saving out
#   an additional GC diagnostic variable -- something like the mass column in addition to PEDGE.
# - Need to double-check units of Jacobian [mixing ratio, unitless] vs units of virtual TROPOMI
#   column [ppb] in use_AK_to_GC().
# - Is TROPP used for anything?

# ==================================================================================================
#
#                                      Define functions
#
# ==================================================================================================

def save_obj(obj, name):
    """ Save something with Pickle. """
    
    with open(name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# --------------------------------------------------------------------------------------------------

def load_obj(name):
    """ Load something with Pickle. """

    with open( name, 'rb') as f:
        return pickle.load(f)


# --------------------------------------------------------------------------------------------------
    
def read_tropomi_operational(filename):
    """
    Read TROPOMI data and save important variables to dictionary.

    Arguments
        filename [str]  : TROPOMI netcdf data file to read

    Returns
        met      [dict] : Dictionary of important variables from TROPOMI:
	    		    - CH4
 			    - Latitude
			    - Longitude
       		            - QA value
 			    - UTC time
 			    - Local time
			    - Averaging kernel
			    - CH4 prior profile
			    - Dry air subcolumns
			    - Latitude bounds
 			    - Longitude bounds
			    - Vertical pressure profile
    """

    # Initialize dictionary for TROPOMI data
    met = {}
    
    # Store methane, QA, lat, lon
    data=xr.open_dataset(filename, group='PRODUCT')
    met['methane'] = data['methane_mixing_ratio_bias_corrected'].values[0,:,:] # 3245x215
    met['qa_value'] = data['qa_value'].values[0,:,:]                           # 3245x215
    met['longitude'] = data['longitude'].values[0,:,:]                         # 3245x215
    met['latitude'] = data['latitude'].values[0,:,:]                           # 3245x215 
    
    # Store UTC time [****why is this necessary? time_utc is a standard tropomi variable]
    referencetime = data['time'].values
    delta_time = data['delta_time'][0].values                                  # 3245x1
    strdate = []
    if delta_time.dtype == '<m8[ns]':
        strdate = referencetime+delta_time
    elif delta_time.dtype == '<M8[ns]':
        strdate = delta_time
    else:
        print(delta_time.dtype)
        pass
    met['utctime'] = strdate
    
    # Store local time
    # [****This assumes 24 equal-sized (15-degree) time zones, but is that valid?]
    timeshift = np.array(met['longitude']/15*60, dtype=int)                             # Convert to minutes
    localtimes = np.zeros(shape=timeshift.shape, dtype='datetime64[ns]')
    for kk in range(timeshift.shape[0]):
        localtimes[kk,:] = strdate[kk]
    met['localtime'] = localtimes
    data.close()
 
    # Store averaging kernel 
    data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
    met['column_AK'] = data['column_averaging_kernel'].values[0,:,:,::-1]
    data.close()
    
    # Store methane prior profile, dry air subcolumns
    data = xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA')
    met['methane_profile_apriori'] = data['methane_profile_apriori'].values[0,:,:,::-1] # mol m-2
    met['dry_air_subcolumns'] = data['dry_air_subcolumns'].values[0,:,:,::-1]           # mol m-2
    
    # Also get pressure interval and surface pressure for use below
    pressure_interval = data['pressure_interval'].values[0,:,:]/100                     # Pa -> hPa
    surface_pressure = data['surface_pressure'].values[0,:,:]/100                       # Pa -> hPa
    data.close()

    # Store lat, lon bounds for pixels
    data=xr.open_dataset(filename, group='PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
    met['longitude_bounds'] = data['longitude_bounds'].values[0,:,:,:]
    met['latitude_bounds'] = data['latitude_bounds'].values[0,:,:,:]
    data.close()
    
    # Store vertical pressure profile
    N1 = met['methane'].shape[0]
    N2 = met['methane'].shape[1]
    pressures = np.zeros([N1,N2,13], dtype=np.float)
    pressures.fill(np.nan)
    for i in range(12+1):
        pressures[:,:,i] = surface_pressure - i*pressure_interval
    met['pressures'] = pressures
    
    # Done
    return met    


# --------------------------------------------------------------------------------------------------

def read_tropomi(filename):
    """
    The same as read_tropomi_operational() above, but for Hannah's processed research product files.
    The fields of the output dictionary have different shapes than from read_tropomi_operational()
    """

    met = {}
    data = xr.open_dataset(filename)
    data.close()
    met['methane'] = data['xch4_corrected'].values# 51975
    met['qa_value'] = data['qa_value'].values# 51975
    met['longitude'] = data['longitude_center'].values# 51975
    met['latitude'] = data['latitude_center'].values# 51975
    met['surface_albedo'] = data['surface_albedo'].values
    met['aerosol_optical_thickness'] = data['aerosol_optical_thickness'].values
    met['surface_altitude'] = data['surface_altitude'].values
    #met['surface_altitude_stdv'] = data['surface_altitude_stdv'].values
    dates = pd.DataFrame(data['time'].values[:,:-1],
                         columns=['year', 'month', 'day', 'hour', 'minute', 'second'])
    dates = pd.to_datetime(dates) #.dt.strftime('%Y%m%dT%H%M%S')
    met['utctime'] = dates
    met['localtime'] = dates
    met['column_AK'] = data['xch4_column_averaging_kernel'].values[:,::-1]#51975x12    
    met['methane_profile_apriori'] = data['ch4_profile_apriori'].values[:,::-1]
    pressure_interval = data['dp'].values #hPa
    surface_pressure = data['surface_pressure'].values #hPa
    met['dry_air_subcolumns'] = data['dry_air_subcolumns'].values[:,::-1]
    met['longitude_bounds'] = data['longitude_corners'].values
    met['latitude_bounds'] = data['latitude_corners'].values
    N1 = met['methane'].shape[0]
    pressures = np.zeros([N1,13], dtype=np.float)
    pressures.fill(np.nan)
    for i in range(12+1):
        pressures[:,i] = surface_pressure - i*pressure_interval
    met['pressures'] = pressures
    
    return met


# --------------------------------------------------------------------------------------------------

def read_GC(date, use_Sensi=False, Sensi_datadir=None, correct_strato=False, lat_mid=None, lat_ratio=None):
    """
    Read GEOS-Chem data and save important variables to dictionary.

    Arguments
        date           [str]   : date of interest
        use_Sensi      [log]   : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        Sensi_datadir  [str]   : If use_Sensi=True, this is the path to the GC sensitivity data
        correct_strato [log]   : Are we doing a latitudinal correction of GEOS-Chem stratospheric bias? 
        lat_mid        [float] : If correct_strato=True, this is the center latitude of each grid box
        lat_ratio      [float] : If correct_strato=True, this is the ratio for correcting GC stratospheric methane to match ACE-FTS
    
    Returns
        met           [dict]  : Dictionary of important variables from GEOS-Chem:
			          - CH4
			          - Latitude
			          - Longitude
			          - PEDGE
			          - TROPP, if correct_strato=True
                                  - CH4_adjusted, if correct_strato=True
    """
    
    # Assemble file paths to GC output collections for input data
    month = int(date[4:6])
    file_species = "GEOSChem.SpeciesConc."+date+"00z.nc4"        
    file_pedge = "GEOSChem.LevelEdgeDiags."+date+"00z.nc4"    
    file_troppause = "GEOSChem.StateMet."+date+"00z.nc4"    
    
    # Read lat, lon, CH4 from the SpeciecConc collection
    filename = GC_datadir+'/'+file_species
    data = xr.open_dataset(filename)
    LON = data['lon'].values
    LAT = data['lat'].values
    CH4 = data['SpeciesConc_CH4'].values[0,:,:,:];
    CH4 = CH4*1e9                                    # Convert to ppb
    CH4 = np.einsum('lij->jil', CH4)
    data.close()

    # Read PEDGE from the LevelEdgeDiags collection
    filename = GC_datadir+'/'+file_pedge
    data = xr.open_dataset(filename)
    PEDGE = data['Met_PEDGE'].values[0,:,:,:];
    PEDGE = np.einsum('lij->jil', PEDGE)
    data.close()
    
    # If want to do latitudinal correction of stratospheric bias in GEOS-Chem:
    if correct_strato:
        
        # Read tropopause level from StateMet collection
        filename = GC_datadir+'/'+file_troppause
        data = xr.open_dataset(filename)
        TROPP = data['Met_TropLev'].values[0,:,:];
        TROPP = np.einsum('ij->ji', TROPP)
        data.close()

        CH4_adjusted = CH4.copy()
        for i in range(len(LON)):
            for j in range(len(LAT)):
                l = int(TROPP[i,j])
                ind = np.where(lat_mid == LAT[j])[0][0]    # Find the location of lat in lat_mid        
                CH4_adjusted[i,j,l:] = CH4[i,j,l:]*lat_ratio[ind,month-1]
    
    # Store GC data in dictionary
    met = {}
    met['lon'] = LON
    met['lat'] = LAT
    met['PEDGE'] = PEDGE
    met['CH4'] = CH4
    if correct_strato:
        met['TROPP'] = TROPP
        met['CH4_adjusted'] = CH4_adjusted
    
    # If need to construct Jacobian, read sensitivity data from GEOS-Chem perturbation simulations
    if use_Sensi:
        filename = Sensi_datadir+'/Sensi_'+date+'.nc'
        data = xr.open_dataset(filename)
        Sensi = data['Sensitivities'].values
        Sensi = np.einsum('klji->ijlk', Sensi)
        data.close()
        # Now adjust the Sensitivity
        Sensi = Sensi*2              # Because we perturb the emissions by 50% [****Only makes sense if optimize scale factors]
        met['Sensi'] = Sensi

    # Done
    return met


# --------------------------------------------------------------------------------------------------

def read_all_GC(all_strdate, use_Sensi=False, Sensi_datadir=None, correct_strato=False, lat_mid=None, lat_ratio=None):
    """ 
    Call read_GC() for multiple dates in a loop. 

    Arguments
    Same as read_GC(), except instead of 'date' argument, use 
        allstr_date [list, str] : date strings 

    Returns
        met         [dict]      : Dictionary of dictionaries. Each sub-dictionary is returned by read_GC()
    """

    met={}
    for strdate in all_strdate:
        met[strdate] = read_GC(strdate, use_Sensi, Sensi_datadir, correct_strato, lat_mid, lat_ratio) 
    
    return met


# --------------------------------------------------------------------------------------------------

def cal_weights(Sat_p, GC_p):
    """
    Calculate pressure weights for TROPOMI & GC
    [****Need more detailed description. What exactly does the output dictionary look like?]

    Arguments
        Sat_p   [float]    : Pressure edge from TROPOMI (13)          <--- 13-1 = 12 pressure levels
        GC_p    [float]    : Pressure edge from GEOS-Chem (48)        <--- 48-1 = 47 pressure levels

    Returns
        weights [dict]     : Pressure weights
"""

    # Step 1: Combine Sat_p and GC_p 
    #         [****Possible to do this in fewer lines?]
    Com_p = np.zeros(len(Sat_p)+len(GC_p));                     # <--- 61-1 = 60 pressure levels
    Com_p.fill(np.nan)
    data_type = np.zeros(len(Sat_p)+len(GC_p), dtype=int); 
    data_type.fill(-99)
    location = []
    i=0;j=0;k=0
    while ((i < len(Sat_p)) or (j < len(GC_p))):
        if i == len(Sat_p):
            Com_p[k] = GC_p[j]
            data_type[k] = 2
            j=j+1; k=k+1
            continue
        if j == len(GC_p):
            Com_p[k] = Sat_p[i]
            data_type[k] = 1  
            location.append(k)        
            i=i+1; k=k+1   
            continue
        if Sat_p[i] >= GC_p[j]:
            Com_p[k] = Sat_p[i]
            data_type[k] = 1        
            location.append(k)        
            i=i+1; k=k+1
        else:
            Com_p[k] = GC_p[j]
            data_type[k] = 2
            j=j+1; k=k+1
    
    # Step 2: Find the first level of GC
    first_2 = -99
    for i in range(len(Sat_p)+len(GC_p)-1):
        if data_type[i] == 2:
            first_2 = i
            break    
    
    # Save data to dictionary
    weights = {}
    weights['data_type'] = data_type
    weights['Com_p'] = Com_p
    weights['location'] = location
    weights['first_2'] = first_2
    
    # Done
    return weights


# --------------------------------------------------------------------------------------------------

def remap(GC_CH4, data_type, Com_p, location, first_2):
    """
    Remap GEOS-Chem methane to TROPOMI vertical grid.

    Arguments
        GC_CH4    [float]   : Methane from GEOS-Chem (60) <---- 60 pressure levels?
        data_type [int]     : ???
        Com_p     [float]   : Combined TROPOMI + GEOS-Chem pressure levels
        location  [int]     : ???
        first_2   [int]     : ???

    Returns
        Sat_CH4   [float]   : GC methane in TROPOMI pressure coordinates
    """
    
    # ****Step 3: ???
    conc = np.zeros(len(Com_p)-1,); conc.fill(np.nan)
    k=0
    for i in range(first_2, len(Com_p)-1):
        conc[i] = GC_CH4[k]
        if data_type[i+1] == 2:
            k=k+1
    if first_2 > 0:
        conc[:first_2] = conc[first_2]
    
    # Step 4: Calculate the weighted mean methane for each layer
    delta_p = Com_p[:-1] - Com_p[1:]
    Sat_CH4 = np.zeros(12); Sat_CH4.fill(np.nan)
    for i in range(len(location)-1):
        start = location[i]
        end = location[i+1]
        fenzi = sum(conc[start:end]*delta_p[start:end])
        fenmu = sum(delta_p[start:end])
        Sat_CH4[i] = fenzi/fenmu   
    
    # Done
    return Sat_CH4


# --------------------------------------------------------------------------------------------------

def remap2(Sensi, data_type, Com_p, location, first_2):
    """
    Remap GEOS-Chem sensitivity data (from perturbation simulations) to TROPOMI vertical grid.

    Arguments
        Sensi     [float]   : 4D sensitivity data from GC perturbation runs, with dims (lon, lat, alt, cluster)****check dim order
        data_type [int]     : ???
        Com_p     [float]   : Combined TROPOMI + GEOS-Chem pressure levels
        location  [int]     : ???
        first_2   [int]     : ???

    Returns
        Sat_CH4   [float]   : GC methane in TROPOMI pressure coordinates
    """
    
    # ****Step 3: ???
    MM = Sensi.shape[1]
    conc = np.zeros((len(Com_p)-1, MM))
    conc.fill(np.nan)
    k=0
    for i in range(first_2, len(Com_p)-1):
        conc[i,:] = Sensi[k,:]
        if data_type[i+1] == 2:
            k=k+1
    if first_2 > 0:
        conc[:first_2,:] = conc[first_2,:]
    
    # Step 4: Calculate the weighted mean methane for each layer
    delta_p = Com_p[:-1] - Com_p[1:]
    delta_ps = np.transpose(np.tile(delta_p, (MM,1)))
    Sat_CH4 = np.zeros((12, MM)); Sat_CH4.fill(np.nan)
    for i in range(len(location)-1):
        start = location[i]
        end = location[i+1]
        fenzi = np.sum(conc[start:end,:]*delta_ps[start:end,:],0)
        fenmu = np.sum(delta_p[start:end])
        Sat_CH4[i,:] = fenzi/fenmu
    
    # Done
    return Sat_CH4


# --------------------------------------------------------------------------------------------------

def nearest_loc(loc_query, loc_grid, tolerance=0.5):
    """ Find the index of the nearest grid location to a query location, with some tolerance. """

    distances = np.abs(loc_grid - loc_query)
    ind = distances.argmin()
    if distances[ind] >= tolerance:
        return np.nan
    else:
        return ind


# --------------------------------------------------------------------------------------------------

def use_AK_to_GC_operational(filename, n_clust, GC_startdate, GC_enddate, xlim, ylim, use_Sensi, Sensi_datadir, correct_strato=False, lat_mid=None, lat_ratio=None):
    """
    THIS VERSION OF THE FUNCTION USES THE OPERATIONAL TROPOMI PRODUCT DATA

    Map GEOS-Chem data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        n_clust        [int]        : Number of clusters / state vector elements
        GC_startdate   [datetime64] : First day of inversion period, for GC and TROPOMI
        GC_enddate     [datetime64] : Last day of inversion period, for GC and TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        use_Sensi      [log]        : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        Sensi_datadir  [str]        : If use_Sensi=True, this is the path to the GC sensitivity data
        correct_strato [log]        : Are we doing a latitudinal correction of GEOS-Chem stratospheric bias?
        lat_mid        [float]      : If correct_strato=True, this is the center latitude of each grid box
        lat_ratio      [float]      : If correct_strato=True, this is the ratio for correcting GC stratospheric methane to match ACE-FTS

    Returns
        result         [dict]       : Dictionary with one or two fields:
 	   	 		        - obs_GC : GEOS-Chem and TROPOMI methane data
	 					     - TROPOMI methane
						     - GC methane
						     - TROPOMI lat, lon
						     - TROPOMI lat index, lon index
				      If use_Sensi=True, also include:
				        - KK     : Jacobian matrix
    """
    
    # Read TROPOMI data
    TROPOMI = read_tropomi(filename)
    
    # We're only going to consider data within lat/lon/time bounds, and with QA > 0.5
    sat_ind = np.where((TROPOMI['longitude'] >  xlim[0])      & (TROPOMI['longitude'] <  xlim[1])     & 
                       (TROPOMI['latitude']  >  ylim[0])      & (TROPOMI['latitude']  <  ylim[1])     & 
                       (TROPOMI['localtime'] >= GC_startdate) & (TROPOMI['localtime'] <= GC_enddate)  &
                       (TROPOMI['qa_value']  >= 0.5))
    # [****Why are we using local times here? Shouldn't we be using UTC?]

    # Number of TROPOMI observations
    NN = len(sat_ind[0])
    print('Found', NN, 'TROPOMI observations.')

    # If need to build Jacobian from GC perturbation simulation sensitivity data:
    if use_Sensi:
        # Initialize Jacobian K
        temp_KK = np.zeros([NN,n_clust], dtype=np.float32)
        temp_KK.fill(np.nan)
    
    # Initialize array with NN rows, 6 columns: TROPOMI-CH4, GC-CH4, longitude, latitude, II, JJ
    temp_obs_GC = np.zeros([NN,6], dtype=np.float32)
    temp_obs_GC.fill(np.nan)
    
    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # For each TROPOMI observation
    for iNN in range(NN):
        
        # Get the date and hour
        iSat = sat_ind[0][iNN] # lat index
        jSat = sat_ind[1][iNN] # lon index
        timeshift = int(TROPOMI['longitude'][iSat,jSat]/15*60)
        timeshift = 0 # Now I use UTC time instead of local time [****Why?]
        localtime = TROPOMI['utctime'][iSat] + np.timedelta64(timeshift, 'm') # local time
        localtime = pd.to_datetime(str(localtime))
        strdate = localtime.round('60min').strftime("%Y%m%d_%H")        
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_GC = read_all_GC(all_strdate, use_Sensi, Sensi_datadir, correct_strato, lat_mid, lat_ratio)
   
    # For each TROPOMI observation: 
    for iNN in range(NN):
        
        # Get GC data for the date of the observation:
        iSat = sat_ind[0][iNN]
        jSat = sat_ind[1][iNN]
        Sat_p = TROPOMI['pressures'][iSat,jSat,:]
        dry_air_subcolumns = TROPOMI['dry_air_subcolumns'][iSat,jSat,:]       # mol m-2
        priori = TROPOMI['methane_profile_apriori'][iSat,jSat,:]              # mol m-2
        AK = TROPOMI['column_AK'][iSat,jSat,:]
        timeshift = int(TROPOMI['longitude'][iSat,jSat]/15*60)
        timeshift = 0
        localtime = TROPOMI['utctime'][iSat]+np.timedelta64(timeshift,'m')    # local time
        localtime = pd.to_datetime(str(localtime))
        strdate = localtime.round('60min').strftime("%Y%m%d_%H")        
        GC = all_date_GC[strdate]
                
        # Find GC lats & lons closest to the corners of the TROPOMI pixel
        longitude_bounds = TROPOMI['longitude_bounds'][iSat,jSat,:]
        latitude_bounds = TROPOMI['latitude_bounds'][iSat,jSat,:]
        corners_lon = [];
        corners_lat = []
        for k in range(4):
            iGC = nearest_loc(longitude_bounds[k], GC['lon'])
            jGC = nearest_loc(latitude_bounds[k], GC['lat'])
            corners_lon.append(iGC)
            corners_lat.append(jGC)
        GC_ij = [(x,y) for x in set(corners_lon) for y in set(corners_lat)]
        GC_grids = [(GC['lon'][i], GC['lat'][j]) for i,j in GC_ij]
        
        # Compute the overlapping area between the TROPOMI pixel and GC grid cells it touches
        overlap_area = np.zeros(len(GC_grids))
        dlon = GC['lon'][1]-GC['lon'][0]
        dlat = GC['lat'][1]-GC['lat'][0]
        # Polygon representing TROPOMI pixel
        p0 = Polygon(np.column_stack((longitude_bounds,latitude_bounds)))
        # For each GC grid cell that touches the TROPOMI pixel: 
        for ipixel in range(len(GC_grids)):
            # Define polygon representing the GC grid cell
            item = GC_grids[ipixel]
            ap1 = [item[0]-dlon/2, item[0]+dlon/2, item[0]+dlon/2, item[0]-dlon/2]
            ap2 = [item[1]-dlat/2, item[1]-dlat/2, item[1]+dlat/2, item[1]+dlat/2]        
            p2 = Polygon(np.column_stack((ap1, ap2)))
            # Calculate overlapping area as the intersection of the two polygons
            if p2.intersects(p0):
                  overlap_area[ipixel] = p0.intersection(p2).area
        
        # If there is no overlap between GC and TROPOMI, skip to next observation:
        if sum(overlap_area) == 0:
            continue                  

        # =======================================================
        #       Map GEOS-Chem to TROPOMI observation space
        # =======================================================
        
        # Otherwise, initialize tropomi virtual xch4 and virtual sensitivity as zero
        GC_base_posteri = 0   # virtual tropomi xch4
        GC_base_sensi = 0     # virtual tropomi sensitivity
        
        # For each GC ground cell that touches the TROPOMI pixel: 
        for ipixel in range(len(GC_grids)):
            
            # Get GC lat/lon indices for the cell
            iGC,jGC = GC_ij[ipixel]
            # Get GC pressure edges for the cell
            GC_p = GC['PEDGE'][iGC,jGC,:]
            # Get GC methane for the cell
            if correct_strato:
                GC_CH4 = GC['CH4_adjusted'][iGC,jGC,:]
            else:
                GC_CH4 = GC['CH4'][iGC,jGC,:]
            # Calculate pressure weights for the cell
            ww = cal_weights(Sat_p, GC_p)
            # Map GC methane to TROPOMI pressure levels
            Sat_CH4 = remap(GC_CH4, ww['data_type'], ww['Com_p'], ww['location'], ww['first_2'])          # ppb
            # Convert ppb to mol m-2
            Sat_CH4_2 = Sat_CH4 * 1e-9 * dry_air_subcolumns                                               # mol m-2
            # Derive the column-averaged XCH4 that TROPOMI would see over this ground cell
            tropomi_sees_ipixel = sum(priori + AK * (Sat_CH4_2-priori)) / sum(dry_air_subcolumns) * 1e9   # ppb
            # Weight by overlapping area (to be divided out later) and add to sum
            GC_base_posteri += overlap_area[ipixel] * tropomi_sees_ipixel                                 # ppb m2

            # If building Jacobian matrix from GC perturbation simulation sensitivity data:
            if use_Sensi:            
                # Get GC perturbation sensitivities at this lat/lon, for all vertical levels and clusters
                Sensi = GC['Sensi'][iGC,jGC,:,:]
                # Map the sensitivities to TROPOMI pressure levels
                Sat_CH4 = remap2(Sensi, ww['data_type'], ww['Com_p'], ww['location'], ww['first_2'])      # mixing ratio, unitless 
                # Tile the TROPOMI averaging kernel
                AKs = np.transpose(np.tile(AK, (n_clust,1)))
                # Tile the TROPOMI dry air subcolumns
                dry_air_subcolumns_s = np.transpose(np.tile(dry_air_subcolumns, (n_clust,1)))             # mol m-2
                # Derive the change in column-averaged XCH4 that TROPOMI would see over this ground cell
                ap = np.sum(AKs*Sat_CH4*dry_air_subcolumns_s, 0) / sum(dry_air_subcolumns)                # mixing ratio, unitless
                # Weight by overlapping area (to be divided out later) and add to sum
                GC_base_sensi += overlap_area[ipixel] * ap                                                # m2

        # Compute virtual TROPOMI observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_tropomi = GC_base_posteri/sum(overlap_area)
 
        # Save real and virtual TROPOMI data                        
        temp_obs_GC[iNN,0] = TROPOMI['methane'][iSat,jSat]     # TROPOMI methane
        temp_obs_GC[iNN,1] = virtual_tropomi                   # GC virtual TROPOMI methane
        temp_obs_GC[iNN,2] = TROPOMI['longitude'][iSat,jSat]   # TROPOMI longitude
        temp_obs_GC[iNN,3] = TROPOMI['latitude'][iSat,jSat]    # TROPOMI latitude 
        temp_obs_GC[iNN,4] = iSat                              # TROPOMI index of longitude
        temp_obs_GC[iNN,5] = jSat                              # TROPOMI index of lattitude
        if use_Sensi:
            # Compute TROPOMI sensitivity as weighted mean by overlapping area
            # i.e., need to divide out area [m2] from the previous step
            temp_KK[iNN,:] = GC_base_sensi/sum(overlap_area)
    
    # Output 
    result = {}
    result['obs_GC'] = temp_obs_GC  # Always return the coincident TROPOMI and GC data
    if use_Sensi:
        result['KK'] = temp_KK      # Optionally return the Jacobian
    
    # Done    
    return result


# --------------------------------------------------------------------------------------------------

def use_AK_to_GC(filename, n_clust, GC_startdate, GC_enddate, xlim, ylim, use_Sensi, Sensi_datadir, correct_strato=False, lat_mid=None, lat_ratio=None):
    """
    THIS VERSION OF THE FUNCTION USES HANNAH'S PROCESSED RESEARCH PRODUCT DATA

    Map GEOS-Chem data to TROPOMI observation space.

    Arguments
        filename       [str]        : TROPOMI netcdf data file to read
        n_clust        [int]        : Number of clusters / state vector elements
        GC_startdate   [datetime64] : First day of inversion period, for GC and TROPOMI
        GC_enddate     [datetime64] : Last day of inversion period, for GC and TROPOMI
        xlim           [float]      : Longitude bounds for simulation domain
        ylim           [float]      : Latitude bounds for simulation domain
        use_Sensi      [log]        : Are we trying to map GEOS-Chem sensitivities to TROPOMI observation space?
        Sensi_datadir  [str]        : If use_Sensi=True, this is the path to the GC sensitivity data
        correct_strato [log]        : Are we doing a latitudinal correction of GEOS-Chem stratospheric bias?
        lat_mid        [float]      : If correct_strato=True, this is the center latitude of each grid box
        lat_ratio      [float]      : If correct_strato=True, this is the ratio for correcting GC stratospheric methane to match ACE-FTS

    Returns
        result         [dict]       : Dictionary with one or two fields:
 	   	 		        - obs_GC : GEOS-Chem and TROPOMI methane data
	 					     - TROPOMI methane
						     - GC methane
						     - TROPOMI lat, lon
						     - TROPOMI lat index, lon index
				      If use_Sensi=True, also include:
				        - KK     : Jacobian matrix
    """
    
    # Read TROPOMI data
    TROPOMI = read_tropomi(filename)
    
    # We're only going to consider data within lat/lon/time bounds, and with QA > 0.5
    sat_ind = np.where((TROPOMI['longitude'] >  xlim[0])      & (TROPOMI['longitude'] <  xlim[1])     & 
                       (TROPOMI['latitude']  >  ylim[0])      & (TROPOMI['latitude']  <  ylim[1])     & 
                       (TROPOMI['localtime'] >= GC_startdate) & (TROPOMI['localtime'] <= GC_enddate)  &
                       (TROPOMI['qa_value']  >= 0.5)          & (TROPOMI['methane'] <= 3000)) # added this methane condition because Lu did...
    # [****localtime is same as utctime here. Should better streamline code.]

    # Number of TROPOMI observations
    NN = len(sat_ind[0])
    print('Found', NN, 'TROPOMI observations.')

    # If need to build Jacobian from GC perturbation simulation sensitivity data:
    if use_Sensi:
        # Initialize Jacobian K
        temp_KK = np.zeros([NN,n_clust], dtype=np.float32)
        temp_KK.fill(np.nan)
    
    # Initialize array with NN rows, 5 columns: TROPOMI-CH4, GC-CH4, longitude, latitude, II
    temp_obs_GC = np.zeros([NN,5], dtype=np.float32)
    temp_obs_GC.fill(np.nan)
    
    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # For each TROPOMI observation
    for iNN in range(NN):
        
        # Get the date and hour
        iSat = sat_ind[0][iNN] # satellite observation index
        localtime = TROPOMI['utctime'][iSat]
        localtime = pd.to_datetime(str(localtime))
        strdate = localtime.round('60min').strftime("%Y%m%d_%H")        
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_GC = read_all_GC(all_strdate, use_Sensi, Sensi_datadir, correct_strato, lat_mid, lat_ratio)
   
    # For each TROPOMI observation: 
    for iNN in range(NN):
        
        # Get GC data for the date of the observation:
        iSat = sat_ind[0][iNN]
        Sat_p = TROPOMI['pressures'][iSat,:]
        dry_air_subcolumns = TROPOMI['dry_air_subcolumns'][iSat,:]       # mol m-2
        priori = TROPOMI['methane_profile_apriori'][iSat,:]              # mol m-2
        AK = TROPOMI['column_AK'][iSat,:]
        localtime = TROPOMI['utctime'][iSat]                             # local time
        localtime = pd.to_datetime(str(localtime))
        strdate = localtime.round('60min').strftime("%Y%m%d_%H")        
        GC = all_date_GC[strdate]
                
        # Find GC lats & lons closest to the corners of the TROPOMI pixel
        longitude_bounds = TROPOMI['longitude_bounds'][iSat,:]
        latitude_bounds = TROPOMI['latitude_bounds'][iSat,:]
        corners_lon = [];
        corners_lat = []
        for k in range(4):
            iGC = nearest_loc(longitude_bounds[k], GC['lon'])
            jGC = nearest_loc(latitude_bounds[k], GC['lat'])
            corners_lon.append(iGC)
            corners_lat.append(jGC)
        GC_ij = [(x,y) for x in set(corners_lon) for y in set(corners_lat)]
        GC_grids = [(GC['lon'][i], GC['lat'][j]) for i,j in GC_ij]
        
        # Compute the overlapping area between the TROPOMI pixel and GC grid cells it touches
        overlap_area = np.zeros(len(GC_grids))
        dlon = GC['lon'][1]-GC['lon'][0]
        dlat = GC['lat'][1]-GC['lat'][0]
        # Polygon representing TROPOMI pixel
        p0 = Polygon(np.column_stack((longitude_bounds,latitude_bounds)))
        # For each GC grid cell that touches the TROPOMI pixel: 
        for ipixel in range(len(GC_grids)):
            # Define polygon representing the GC grid cell
            item = GC_grids[ipixel]
            ap1 = [item[0]-dlon/2, item[0]+dlon/2, item[0]+dlon/2, item[0]-dlon/2]
            ap2 = [item[1]-dlat/2, item[1]-dlat/2, item[1]+dlat/2, item[1]+dlat/2]        
            p2 = Polygon(np.column_stack((ap1, ap2)))
            # Calculate overlapping area as the intersection of the two polygons
            if p2.intersects(p0):
                  overlap_area[ipixel] = p0.intersection(p2).area
        
        # If there is no overlap between GC and TROPOMI, skip to next observation:
        if sum(overlap_area) == 0:
            continue                  

        # =======================================================
        #       Map GEOS-Chem to TROPOMI observation space
        # =======================================================
        
        # Otherwise, initialize tropomi virtual xch4 and virtual sensitivity as zero
        GC_base_posteri = 0   # virtual tropomi xch4
        GC_base_sensi = 0     # virtual tropomi sensitivity
        
        # For each GC ground cell that touches the TROPOMI pixel: 
        for ipixel in range(len(GC_grids)):
            
            # Get GC lat/lon indices for the cell
            iGC,jGC = GC_ij[ipixel]
            # Get GC pressure edges for the cell
            GC_p = GC['PEDGE'][iGC,jGC,:]
            # Get GC methane for the cell
            if correct_strato:
                GC_CH4 = GC['CH4_adjusted'][iGC,jGC,:]
            else:
                GC_CH4 = GC['CH4'][iGC,jGC,:]
            # Calculate pressure weights for the cell
            ww = cal_weights(Sat_p, GC_p)
            # Map GC methane to TROPOMI pressure levels
            Sat_CH4 = remap(GC_CH4, ww['data_type'], ww['Com_p'], ww['location'], ww['first_2'])          # ppb
            # Convert ppb to mol m-2
            Sat_CH4_2 = Sat_CH4 * 1e-9 * dry_air_subcolumns                                               # mol m-2
            # Derive the column-averaged XCH4 that TROPOMI would see over this ground cell
            tropomi_sees_ipixel = sum(priori + AK * (Sat_CH4_2-priori)) / sum(dry_air_subcolumns) * 1e9   # ppb
            # Weight by overlapping area (to be divided out later) and add to sum
            GC_base_posteri += overlap_area[ipixel] * tropomi_sees_ipixel                                 # ppb m2

            # If building Jacobian matrix from GC perturbation simulation sensitivity data:
            if use_Sensi:            
                # Get GC perturbation sensitivities at this lat/lon, for all vertical levels and clusters
                Sensi = GC['Sensi'][iGC,jGC,:,:]
                # Map the sensitivities to TROPOMI pressure levels
                Sat_CH4 = remap2(Sensi, ww['data_type'], ww['Com_p'], ww['location'], ww['first_2'])      # mixing ratio, unitless 
                # Tile the TROPOMI averaging kernel
                AKs = np.transpose(np.tile(AK, (n_clust,1)))
                # Tile the TROPOMI dry air subcolumns
                dry_air_subcolumns_s = np.transpose(np.tile(dry_air_subcolumns, (n_clust,1)))             # mol m-2
                # Derive the change in column-averaged XCH4 that TROPOMI would see over this ground cell
                ap = np.sum(AKs*Sat_CH4*dry_air_subcolumns_s, 0) / sum(dry_air_subcolumns)                # mixing ratio, unitless
                # Weight by overlapping area (to be divided out later) and add to sum
                GC_base_sensi += overlap_area[ipixel] * ap                                                # m2

        # Compute virtual TROPOMI observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_tropomi = GC_base_posteri/sum(overlap_area)
 
        # Save real and virtual TROPOMI data                        
        temp_obs_GC[iNN,0] = TROPOMI['methane'][iSat]     # TROPOMI methane
        temp_obs_GC[iNN,1] = virtual_tropomi              # GC virtual TROPOMI methane
        temp_obs_GC[iNN,2] = TROPOMI['longitude'][iSat]   # TROPOMI longitude
        temp_obs_GC[iNN,3] = TROPOMI['latitude'][iSat]    # TROPOMI latitude 
        temp_obs_GC[iNN,4] = iSat                         # TROPOMI index of observation
        if use_Sensi:
            # Compute TROPOMI sensitivity as weighted mean by overlapping area
            # i.e., need to divide out area [m2] from the previous step
            temp_KK[iNN,:] = GC_base_sensi/sum(overlap_area)
    
    # Output 
    result = {}
    result['obs_GC'] = temp_obs_GC  # Always return the coincident TROPOMI and GC data
    if use_Sensi:
        result['KK'] = temp_KK      # Optionally return the Jacobian
    
    # Done    
    return result


# ==================================================================================================
#
#                                      Run the code
#
# ==================================================================================================

if __name__ == '__main__':
    import sys

    startday = sys.argv[1]
    endday = sys.argv[2]
    workdir = sys.argv[3]
    outputdir = sys.argv[4]

    outputdir = outputdir+'/'
 
    # Reformat start and end days for datetime in configuration
    start = startday[0:4]+'-'+startday[4:6]+'-'+startday[6:8]+' 00:00:00'
    end = endday[0:4]+'-'+endday[4:6]+'-'+endday[6:8]+' 23:59:59'

    # Configuration
    use_Sensi = True
    correct_strato = False 
    #workdir = "/n/holyscratch01/jacob_lab/dvaron"
    Sensi_datadir = f'{workdir}/Sensi'
    #Sat_datadir = '/n/holyscratch01/jacob_lab/dvaron/data_TROPOMI/' # always use same sat data dir
    #Sat_datadir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/constant_prior/inversion/data_TROPOMI/'
    #Sat_datadir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/blended_prior/inversion/data_TROPOMI_lu_uptomarch2020/'
    #Sat_datadir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/blended_prior/inversion/data_TROPOMI_hannah/'
    #Sat_datadir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/uninformative_prior/inversion/data_TROPOMI_lu_uptomarch2020/'
    #Sat_datadir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/uninformative_prior/inversion/data_TROPOMI_hannah/'
    #Sat_datadir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/nudged_prior/inversion/data_TROPOMI_lu_uptomarch2020/'
    Sat_datadir = '/n/jacob_lab/Lab/seasasfs02/dvaron/GEOSChem/sensitivity_permian/nudged_prior/inversion/data_TROPOMI_hannah/'
    GC_datadir = f'{workdir}/data_GC'
    #outputdir = workdir+"data_converted/"
    n_clust = 235+8
    xlim = [-111,-95]
    ylim = [24,39]
    GC_startdate = np.datetime64(datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S'))
    GC_enddate = np.datetime64(datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(days=1))
    print('Start:', start)
    print('End:', end)

    # Move to Step1 directory
    #os.chdir(workdir+"Step1_convert_GC")

    # Get TROPOMI data filenames for the desired date range
    allfiles = glob.glob(Sat_datadir+'*.nc')
    Sat_files = []
    for index in range(len(allfiles)):
        filename = allfiles[index]
        shortname = re.split('\/', filename)[-1]
        shortname = re.split('\.', shortname)[0]
        strdate = re.split('\.|_+|T',shortname)[4]
        strdate2 = datetime.datetime.strptime(strdate, '%Y%m%d')
        if ((strdate2 >= GC_startdate) and (strdate2 <= GC_enddate)):
            Sat_files.append(filename)
    Sat_files.sort()
    print('Found',len(Sat_files),'TROPOMI data files.')

    # Map GEOS-Chem to TROPOMI observation space; also return Jacobian matrix if use_Sensi=True
    for filename in Sat_files:
        # Check if TROPOMI file has already been processed
        print('========================')
        temp = re.split('\/', filename)[-1]
        print(temp)
        date = re.split('\.',temp)[0]
        # If not yet processed, run use_AK_to_GC()
        if ~os.path.isfile(outputdir+date+'_GCtoTROPOMI.pkl'):
            if correct_strato:
                df = pd.read_csv("./lat_ratio.csv", index_col=0)
                lat_mid = df.index
                lat_ratio = df.values
                result = use_AK_to_GC(filename, n_clust, GC_startdate, GC_enddate, xlim, ylim, use_Sensi, Sensi_datadir, correct_strato, lat_mid, lat_ratio)
            else:
                print('Running use_AK_to_GC().')
                result = use_AK_to_GC(filename, n_clust, GC_startdate, GC_enddate, xlim, ylim, use_Sensi, Sensi_datadir)
        save_obj(result, outputdir+date+'_GCtoTROPOMI.pkl')
