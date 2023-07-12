import numpy as np
import xarray as xr
from sklearn.cluster import KMeans
from bs4 import BeautifulSoup
import requests
import os
from zipfile import ZipFile
import sys
import fiona
import regionmask
import xesmf as xe
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
from requests.exceptions import Timeout
import yaml
import ruamel.yaml


def get_nested_grid_bounds(land_cover_pth):
    """
    Get the lat/lon bounds of the nested grid window for the inversion.
    The land cover file path specifies the window.
    """

    land_cover = xr.load_dataset(land_cover_pth)
    minLat_allowed = np.min(land_cover.lat.values)
    maxLat_allowed = np.max(land_cover.lat.values)
    minLon_allowed = np.min(land_cover.lon.values)
    maxLon_allowed = np.max(land_cover.lon.values)

    return minLat_allowed, maxLat_allowed, minLon_allowed, maxLon_allowed


def check_nested_grid_compatibility(lat_min, lat_max, lon_min, lon_max, land_cover_pth):
    """
    Check whether input lat/lon bounds are compatible with (contained within) the nested grid window.
    The land cover file path specifies the window.
    """

    (
        minLat_allowed,
        maxLat_allowed,
        minLon_allowed,
        maxLon_allowed,
    ) = get_nested_grid_bounds(land_cover_pth)

    if (
        lat_min < minLat_allowed
        or lat_max > maxLat_allowed
        or lon_min < minLon_allowed
        or lon_max > maxLon_allowed
    ):
        compatible = False
    else:
        compatible = True

    return compatible


def make_rectilinear_state_vector_file(
    config_path, land_cover_pth, hemco_diag_pth, save_pth,
):
    """
    Generates the state vector file for an analytical inversion.

    Arguments
        config_path    [str]   : Path to configuration file
        land_cover_pth [str]   : Path to land cover file
        hemco_diag_pth [str]   : Path to initial HEMCO diagnostics file
        save_pth       [str]   : Where to save the state vector file

    Returns
        ds_statevector []     : xarray dataset containing state vector field formatted for HEMCO

    Notes
        - Land cover file looks like 'GEOSFP.20200101.CN.025x03125.NA.nc' (or 0.5-deg equivalent)
        - HEMCO diags file needs to be global, is used to include offshore emissions in state vector
        - Land cover file and HEMCO diags file need to have the same grid resolution
    """

    # Get config
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    lat_min = config["LatMin"]
    lat_max = config["LatMax"]
    lon_min = config["LonMin"]
    lon_max = config["LonMax"]
    buffer_deg = config["BufferDeg"]
    land_threshold = config["LandThreshold"]
    emis_threshold = config["OffshoreEmisThreshold"]
    k_buffer_clust = config["nBufferClusters"]

    # Load land cover data and HEMCO diagnostics
    lc = xr.load_dataset(land_cover_pth)
    hd = xr.load_dataset(hemco_diag_pth)

    # Require hemco diags on same global grid as land cover map
    hd["lon"] = hd["lon"] - 0.03125  # initially offset by 0.03125 degrees

    # Select / group fields together
    lc = (lc["FRLAKE"] + lc["FRLAND"] + lc["FRLANDIC"]).drop("time").squeeze()
    hd = (hd["EmisCH4_Oil"] + hd["EmisCH4_Gas"]).drop("time").squeeze()

    # Check compatibility of region of interest with nesting window
    compatible = check_nested_grid_compatibility(
        lat_min, lat_max, lon_min, lon_max, land_cover_pth
    )
    if not compatible:
        raise ValueError(
            "Region of interest not contained within selected NestedRegion; see config.yml."
        )

    # Define bounds of inversion domain
    (
        minLat_allowed,
        maxLat_allowed,
        minLon_allowed,
        maxLon_allowed,
    ) = get_nested_grid_bounds(land_cover_pth)
    lon_min_inv_domain = np.max([lon_min - buffer_deg, minLon_allowed])
    lon_max_inv_domain = np.min([lon_max + buffer_deg, maxLon_allowed])
    lat_min_inv_domain = np.max([lat_min - buffer_deg, minLat_allowed])
    lat_max_inv_domain = np.min([lat_max + buffer_deg, maxLat_allowed])

    # Subset inversion domain for land cover and hemco diagnostics fields
    lc = lc.isel(lon=lc.lon >= lon_min_inv_domain, lat=lc.lat >= lat_min_inv_domain)
    lc = lc.isel(lon=lc.lon <= lon_max_inv_domain, lat=lc.lat <= lat_max_inv_domain)
    hd = hd.isel(lon=hd.lon >= lon_min_inv_domain, lat=hd.lat >= lat_min_inv_domain)
    hd = hd.isel(lon=hd.lon <= lon_max_inv_domain, lat=hd.lat <= lat_max_inv_domain)

    # Initialize state vector from land cover, replacing all values with NaN (to be filled later)
    statevector = lc.where(lc == -9999.0)

    # Set pixels in buffer areas to 0
    statevector[:, (statevector.lon < lon_min) | (statevector.lon > lon_max)] = 0
    statevector[(statevector.lat < lat_min) | (statevector.lat > lat_max), :] = 0

    # Also set pixels over water to 0, unless there are offshore emissions
    if land_threshold:
        # Where there is neither land nor emissions, replace with 0
        land = lc.where((lc > land_threshold) | (hd > emis_threshold))
        statevector.values[land.isnull().values] = -9999

    # Fill in the remaining NaNs with state vector element values
    statevector.values[statevector.isnull().values] = np.arange(
        1, statevector.isnull().sum() + 1
    )[::-1]

    # Assign buffer pixels (the remaining 0's) to state vector
    # -------------------------------------------------------------------------
    buffer_area = statevector.values == 0

    # Get image coordinates of all pixels in buffer area
    irows = np.arange(buffer_area.shape[0])
    icols = np.arange(buffer_area.shape[1])
    irows = np.transpose(np.tile(irows, (len(icols), 1)))
    icols = np.tile(icols, (len(irows), 1))
    irows_good = irows[buffer_area > 0]
    icols_good = icols[buffer_area > 0]
    coords = [[icols_good[j], irows_good[j]] for j in range(len(irows_good))]

    # K-means
    X = np.array(coords)
    kmeans = KMeans(n_clusters=k_buffer_clust, random_state=0).fit(X)

    # Assign pixels to state vector
    highres_statevector_max = np.nanmax(statevector.values)
    n_rows = statevector.shape[0]
    n_cols = statevector.shape[1]
    for r in range(n_rows):
        for c in range(n_cols):
            if statevector[r, c].values == 0:
                statevector[r, c] = (
                    kmeans.predict([[c, r]])[0] + 1 + highres_statevector_max
                )
    # -------------------------------------------------------------------------

    # Make dataset
    da_statevector = statevector.copy()
    ds_statevector = da_statevector.to_dataset(name="StateVector")

    # Add attribute metadata
    ds_statevector.lat.attrs["units"] = "degrees_north"
    ds_statevector.lat.attrs["long_name"] = "Latitude"
    ds_statevector.lon.attrs["units"] = "degrees_east"
    ds_statevector.lon.attrs["long_name"] = "Longitude"
    ds_statevector.StateVector.attrs["units"] = "none"
    ds_statevector.StateVector.attrs["missing_value"] = -9999
    ds_statevector.StateVector.attrs["_FillValue"] = -9999

    # Save
    if save_pth is not None:
        print("Saving file {}".format(save_pth))
        ds_statevector.to_netcdf(
            save_pth,
            encoding={
                v: {"zlib": True, "complevel": 9} for v in ds_statevector.data_vars
            },
        )

    return ds_statevector



def get_country_code(ROI):
    """
    Matches the given country to its country code by searching through online database

    Arguments
        ROI            [str]  : Name of the country specified in config file

    Returns
        countryCode    [str]  : Three-letter country code (e.g. USA) of the ROI
    """
    geodata_url = "https://gadm.org/download_country.html"  # URL of the database for country shapefiles

    try:
        # Send a GET request to the geodata URL with a timeout of 5 seconds
        page = requests.get(geodata_url, timeout=5)

        # Create a BeautifulSoup parser object to parse the HTML content
        parser = BeautifulSoup(page.content, 'html.parser')

        found = False  # Flag to track if the ROI is found in the database

        # Iterate over each 'option' tag in the parsed HTML content
        for country in parser.find_all('option'):
            if country.string == ROI:  # Check if the country name matches the given ROI
                countryInfo = str(country)
                # Extract the country code from the HTML string
                countryCode = countryInfo.split("=\"")[1].split("_")[0]
                print(f"{ROI} was found in the database with a country code: {countryCode}")  # Print the country code
                found = True

        if not found:
            # If the ROI is not found in the database, print an error message with a link to the download page
            print("Error: The ROI entered does not match any of the regions in our database. "
                  "Please visit https://gadm.org/download_country.html")
            return None

        return countryCode  # Return the country code

    except Timeout:
        # Handle timeout exception
        print("Timeout error: The request to the web server timed out.")
        return None


def get_shapefile(countryCode):
    """
    Fetches the shapefile corresponding to the ROI and downloads it in the current directory

    Arguments
        countryCode    [str]  : Three-letter country code (e.g. USA) of the ROI

    Returns
        shp_name        [str] : Path to the dowloaded shapefile
    """
    shp_url = f"https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_{countryCode}_shp.zip"  # URL for country shapefile

    try:
        # Send a GET request to the shapefile URL with a timeout of 30 seconds
        rshp = requests.get(shp_url, allow_redirects=True, timeout=30)

        file = f"{countryCode}_SHP.zip"
        # Save the shapefile zip content to a file
        open(file, 'wb').write(rshp.content)

        path = f"/{countryCode}_SHP"
        path = str(os.getcwd()) + path
        print(path)

        # Create a directory to extract the shapefile contents
        os.mkdir(path)

        # Extract the shapefile contents from the zip file to the specified directory
        with ZipFile(file, 'r') as zObject:
            zObject.extractall(path=path)

        shp_name = f"{path}/gadm41_{countryCode}_0.shp"
        return shp_name

    except Timeout:
        # Handle timeout exception
        print("Timeout error: The request to download the shapefile timed out.")
        return None


def make_statevector(config, shapefile_path, land_cover_pth, hemco_diag_pth, save_pth):
    """
    Generates the state vector file for an analytical inversion based on the given shapefile.

    Arguments
        config_path    {dict}  : Config file
        shapefile_path [str]   : Path to the shapefile
        land_cover_pth [str]   : Path to land cover file
        hemco_diag_pth [str]   : Path to initial HEMCO diagnostics file
        save_pth       [str]   : Where to save the state vector file

    Returns
        ds_statevector []     : xarray dataset containing state vector field formatted for HEMCO

    Notes
        - Land cover file looks like 'GEOSFP.20200101.CN.025x03125.NA.nc' (or 0.5-deg equivalent)
        - HEMCO diags file needs to be global, is used to include offshore emissions in state vector
        - Land cover file and HEMCO diags file need to have the same grid resolution
        - Shapefile, including the buffer area specified in the config file, must be fully contained within the nested region (e.g. ME)
    """
    shape = gpd.read_file(shapefile_path)

    # Define min/max latitude and longitude for inversion domain
    bufferdeg = config['BufferDeg']
    bounds = shape.bounds
    lon_min = bounds['minx'][0] - bufferdeg
    lon_max = bounds['maxx'][0] + bufferdeg
    lat_min = bounds['miny'][0] - bufferdeg
    lat_max = bounds['maxy'][0] + bufferdeg

    # Create lat/lon grid (can increase dlon, dlat if the regridding in next steps is slow)
    dlon = 0.01
    dlat = 0.01
    lon = np.arange(lon_min, lon_max+dlon, dlon)
    lat = np.arange(lat_min, lat_max+dlat, dlat)

    # Inversion spatial resolution
    if config['Res'] == '0.25x0.3125':
        lat_res = 0.25
        lon_res = 0.3125
    elif config['Res'] == '0.5x0.625':
        lat_res = 0.5
        lon_res = 0.625
    
    # Make mask from shape file
    mask = regionmask.mask_geopandas(shape, lon, lat) + 1 # Add 1 so the mask values are 1 instead of 0 

    # Need to regrid to the grid HEMCO expects
    reference_lat_grid = np.arange(-90 , 90+lat_res , lat_res)
    reference_lon_grid = np.arange(-180, 180+lon_res, lon_res)

    # Find closest reference coordinates to selected lat/lon bounds
    lat_min = reference_lat_grid[np.abs(reference_lat_grid - lat_min).argmin()]
    lon_min = reference_lon_grid[np.abs(reference_lon_grid - lon_min).argmin()]
    lat_max = reference_lat_grid[np.abs(reference_lat_grid - lat_max).argmin()]
    lon_max = reference_lon_grid[np.abs(reference_lon_grid - lon_max).argmin()]

    # Create an xESMF regridder object to resample the mask on the grid HEMCO expects
    new_lat_grid = np.arange(lat_min, lat_max+lat_res, lat_res)
    new_lon_grid = np.arange(lon_min, lon_max+lon_res, lon_res)
    ds_out = xr.Dataset({'lat': (['lat'], new_lat_grid),
                        'lon': (['lon'], new_lon_grid),
                        }
                    )
    #ds_out = ds_out.rename({'longitude': 'lon', 'latitude': 'lat'})

    regridder = xe.Regridder(mask, ds_out, 'nearest_s2d')  

    mask = regridder(mask)


    state_vector = mask.copy()

    # Set pixels in buffer areas to 0
    state_vector.values[mask.isnull()] = 0
    state_vector.values[mask.isnull()] = 0

    lc = xr.load_dataset(land_cover_pth)
    hd = xr.load_dataset(hemco_diag_pth)
    hd["lon"] = hd["lon"] - 0.03125  # initially offset by 0.03125 degrees

    # Group together
    lc = (lc['FRLAKE'] + lc['FRLAND'] + lc['FRLANDIC']).drop('time').squeeze()
    hd = (hd["EmisCH4_Oil"] + hd["EmisCH4_Gas"]).drop("time").squeeze()

    # Subset the area of interest
    lc = lc.isel(lon=lc.lon>=lon_min, lat=lc.lat>=lat_min)
    lc = lc.isel(lon=lc.lon<=lon_max, lat=lc.lat<=lat_max)
    hd = hd.isel(lon=hd.lon >= lon_min, lat=hd.lat >= lat_min)
    hd = hd.isel(lon=hd.lon <= lon_max, lat=hd.lat <= lat_max)

    # Set pixels over water to 0, unless there are offshore emissions
    if config['LandThreshold']:
        # Where there is neither land nor emissions, replace with 0
        land = lc.where((lc > config['LandThreshold']) | (hd > config['OffshoreEmisThreshold']))
        state_vector.values[land.isnull().values] = 0
    
    # Enumerate state vector elements
    n_lat = len(new_lat_grid)
    n_lon = len(new_lon_grid)
    count = 1
    for r in range(n_lat):
        for c in range(n_lon):
            if state_vector[r,c] == 1:
                state_vector[r,c] = count
                count += 1

    # Now set pixels over water to NaN, unless there are offshore emissions
    if config['LandThreshold']:
        # Where there is no land or offshore emissions, replace with NaN
        state_vector = state_vector.where((lc > config['LandThreshold']) | (hd > config['OffshoreEmisThreshold']))


    # Add buffer elements
    buffer_elements = np.abs((mask > 0) - 1).values

    # Get image coordinates of buffer elements
    irows = np.arange(buffer_elements.shape[0])
    icols = np.arange(buffer_elements.shape[1])
    irows = np.transpose(np.tile(irows,(len(icols),1)))
    icols = np.tile(icols,(len(irows),1)) * (buffer_elements > 0)

    # Select image coordinates of buffer elements
    irows_buffer = irows[buffer_elements > 0]
    icols_buffer = icols[buffer_elements > 0]
    coords = [[icols_buffer[j], irows_buffer[j]] for j in range(len(irows_buffer))]

    # Kmeans based on image coordinates
    X = np.array(coords)
    kmeans = KMeans(n_clusters=config['nBufferClusters'], random_state=0).fit(X)

    # Assign labels
    for r in range(n_lat):
        for c in range(n_lon):
            if state_vector[r,c].values == 0:
                state_vector[r,c] = kmeans.predict([[c,r]])[0] + count

    # Add units attribute
    state_vector.attrs['units'] = 'none'

    if config['LandThreshold']:
        state_vector.values[state_vector.isnull()] = -9999
    
    # Create dataset
    ds_state_vector = xr.Dataset(
                    {
                        "StateVector": (["lat", "lon"], state_vector),
                    },
                    coords={
                        "lon": ("lon", new_lon_grid),
                        "lat": ("lat", new_lat_grid),
                    },
                    
                    )

    # Add attribute metadata
    ds_state_vector.lat.attrs['units'] = 'degrees_north'
    ds_state_vector.lat.attrs['long_name'] = 'Latitude'
    ds_state_vector.lon.attrs['units'] = 'degrees_east'
    ds_state_vector.lon.attrs['long_name'] = 'Longitude'
    ds_state_vector.StateVector.attrs['units'] = 'none'
    ds_state_vector.StateVector.attrs['missing_value'] = -9999
    ds_state_vector.StateVector.attrs['_FillValue'] = -9999

    ds_state_vector.to_netcdf(
        save_pth,
        encoding={
            v: {"zlib": True, "complevel": 9} for v in ds_state_vector.data_vars
        },
    )
    return ds_state_vector

def make_country_state_vector_file(config_path, land_cover_pth, hemco_diag_pth, save_pth):
    # Load the configuration file
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    
    # Extract the region of interest (ROI) from the configuration
    roi = config["RegionOfInterest"]
    
    # Get the country code based on the ROI
    code = get_country_code(roi)
    
    # Gets the shapefile for the corrresponding country and saves it in the current directory
    shapefile_path = get_shapefile(code)
    
    # Update the configuration dict with the shapefile path
    config["ShapeFile"] = shapefile_path

    #update the original configuration file using the ruamel.yaml library
    ruamelyaml = ruamel.yaml.YAML()
    ruamelyaml.preserve_quotes = True
    with open(config_path) as fp:
        changed_config = ruamelyaml.load(fp)
    changed_config['ShapeFile'] = config["ShapeFile"]
    with open(config_path, 'w') as fp:
        ruamelyaml.dump(changed_config, fp)
    
    # Print the path where the shapefile is saved
    print(f"Shapefile saved at {shapefile_path}")
    
    # Call the function to create the state vector using the provided configuration and paths
    make_statevector(config, shapefile_path, land_cover_pth, hemco_diag_pth, save_pth)


if __name__ == "__main__":
    import sys

    config_path = sys.argv[1]
    config = yaml.load(open(config_path), Loader=yaml.FullLoader)
    land_cover_pth = sys.argv[2]
    hemco_diag_pth = sys.argv[3]
    save_pth = sys.argv[4]
    # If the Region of Interest is specified, then the make_country_state_vector_file function is called. Else, the make_rectilinear_state_vector_file function is called
    try:
        if (config["RegionOfInterest"] is not None) and config["CreateAutomaticRectilinearStateVectorFile"] == False:
            make_country_state_vector_file(
            config_path, 
            land_cover_pth, 
            hemco_diag_pth, 
            save_pth
            )
        elif (config["RegionOfInterest"] is None) and config["CreateAutomaticRectilinearStateVectorFile"] == True:
            print("Creating Rectilinear StateVector File")
            make_rectilinear_state_vector_file(
            config_path, 
            land_cover_pth, 
            hemco_diag_pth, 
            save_pth,
            )
        else:
            print("ERROR: The method of creating a statevector was not specified correctly.")
    except KeyError:
        print("No specified ROI")
        make_rectilinear_state_vector_file(
        config_path, 
        land_cover_pth, 
        hemco_diag_pth, 
        save_pth,
        )
