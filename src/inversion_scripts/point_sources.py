# This file contains functions related to accessing point source data
import os
import datetime
import requests
import pandas as pd
import geopandas as gpd
from bs4 import BeautifulSoup
from shapely.geometry import Point
from dateutil.relativedelta import relativedelta

import numpy as np
from shapely.geometry import Polygon
from warnings import warn


class PointSources:
    '''
    Point source plume and emissions data filtered to ROI.

    Multiple point source datasets can be included and data
    will be concatenated. Current options include:
        - SRON plume data
        - CarbonMapper plume data
    
    
    Datasources must have for each point
    * emissions in kg/hr
    * number of detections
    * point geometry on WGS84 crs
        

    Parameters
    ----------
    geofilter: GeoFilter 
    datasources: list of point source data. Can
        include CarbonMapper, SRON. Others
        not yet implemented.

    Use
    ---
    ps = PointSources(geofilter, datasources)

    To return a list of [x,y] coordinates:
    ps.get_coords()
    
    To return a list of [x,y] coordinates
    on the state vector grid (i.e. multiple
    points in single grid are given as a
    single coordinate):
    ps.get_gridded_coords()

    '''
    def __init__(self, geofilter, datasources):

        self.geofilter = geofilter
        if (type(datasources) != list):
            datasources = [datasources]
        self.datasources = datasources
        self.points = self._get_points()
        self.grid_ds = self.grid_all_data()

    def _get_points(self):
        # concatenate and return points geometry
        # within the ROI for all point source datasets 
        in_points = []
        for d in self.datasources:
            # True if point in ROI, else False
            inout = self.geofilter.filter_points(d.gdf)
            if inout is not None:
                in_points.append(d.gdf[inout])
        if len(in_points) > 0:
            gdf_combined = pd.concat(in_points)
            # filter out nans emissions
            gdf_filtered = gdf_combined[~gdf_combined.emission.isna()].copy()
            return gdf_filtered.geometry
        else:
            return None
    
    
    def grid_datasource(self, datasource):
        '''
        put all point source data on the state vector
        grid. Need:
         - emissions [kg/hr]
         - n observations
         - n detections
         
        return a dataset
        '''
        
        # 1. filter points to roi
        inout = self.geofilter.filter_points(datasource.gdf)
        
        if inout.sum() > 0:
            gdf_inroi = datasource.gdf[inout]
        else:
            return None
        
        # 2. filter out nans emissions
        gdf_filtered = gdf_inroi[~gdf_inroi.emission.isna()].copy()
        
        # 3. get I,J indices
        lon_dist = gdf_filtered.geometry.x.values[:,None] - self.geofilter.lons[None,:]
        lat_dist = gdf_filtered.geometry.y.values[:,None] - self.geofilter.lats[None,:]
        ilon = np.argmin( np.abs(lon_dist), 1)
        ilat = np.argmin( np.abs(lat_dist), 1)
        gdf_filtered['I'] = ilon
        gdf_filtered['J'] = ilat
        
        # 4. average up to grid, sum count, etc
        dfgb = gdf_filtered.groupby(['I','J'])
        gdf_grid = dfgb.sum()[['detection_count','observation_count']]
        gdf_grid['emission'] = dfgb.mean()[['emission']]
        gdf_grid['persistence'] = gdf_grid.detection_count / gdf_grid.observation_count

        # 5. make it on a grid
        dalist = []
        for v in gdf_grid.columns:
            
            indat = np.full(
                (self.geofilter.lats.shape[0],
                 self.geofilter.lons.shape[0]),
                np.nan
            )
            
            indat[
                gdf_grid.reset_index().J,
                gdf_grid.reset_index().I
            ] = gdf_grid.reset_index()[v]
            
            vda = xr.DataArray(
                data = indat,
                coords = {
                    'lat':self.geofilter.lats,
                    'lon':self.geofilter.lons
                },
                dims = ['lat','lon'],
                name = v
            )
            dalist.append(vda)

        # vars of this are columns of gdf_grid
        ds = xr.merge(dalist)
        
        return ds
    
    
    def grid_all_data(self):
        d_list = []
        d_names = []
        for d in self.datasources:
            ds_i = self.grid_datasource(d)
            d_list.append(ds_i)
            if ds_i is not None:
                d_names.append(d.myname)
        if not all([i is None for i in d_list]):
            ds = xr.concat(d_list, 'observer')
            ds = ds.assign_coords({'observer': d_names})
            return ds
        else:
            msg = 'No point sources in ROI'
            warnings.warn(msg)
            return None
        
    
    def get_coords(self):
        '''
        Get y,x (lat,lon) coordinates for point sources
        from all datasets
         
        Returns: list of [y,x] coordinates in [lon,lat]
        '''
        if self.points is None:
            return []
        else:
            coords = np.column_stack((
                self.points.y.values,
                self.points.x.values
            )).tolist()
            return coords
        
    
    def get_gridded_coords(self):
        '''
        Get y,x (lat,lon) coordinates for point sources
        from all datasets
         
        Returns: list of [y,x] coordinates in [lon,lat]
        '''
        if self.grid_ds is None:
            return []
        else:
            dftmp = self.grid_ds.emission.to_dataframe().reset_index()
            coords = (
                dftmp[~dftmp.emission.isna()][['lat','lon']]
                .values.tolist()
            )
            return coords

class GeoFilter:
    '''
    Geometry of ROI and inversion domain

    Parameters
    ----------
    config: dict of IMI config.yml contents

    Use
    ---
    gf = GeoFilter(config)

    To check whether some other point geometry 
    is in ROI:
    gf.filter_points()

    '''
    def __init__(self, config):
        self.config = config
        self.geo = self._make_roi_geometry()
        self.svds = self._get_state_vector_file()
        self.lons = self.svds.lon.values
        self.lats = self.svds.lat.values

    def _get_state_vector_file(self):
        
        # infer lat lons from state vector file
        # this file should already exist
        svf = f'{config["OutputPath"]}/{config["RunName"]}/StateVector.nc'
        
        try:
            svds = xr.open_dataset(svf)
        except FileNotFoundError:
            msg = (
                f'State vector file {svf} '
                'cannot be found. Re-run IMI '
                'with "RunSetup" True.'
            )
            sys.exit(msg)
        
        return svds

        
        
    def _make_roi_geometry(self):
        
        custom_vectorfile = not self.config["CreateAutomaticRectilinearStateVectorFile"]
        
        if custom_vectorfile:
            shapefile_path = self.config["ShapeFile"]
            shp_geo = gpd.read_file(shapefile_path).geometry
            if shp_geo.shape[0] > 1:
                msg = (
                    'Shapefile has >1 shape, only using point sources in '
                    'first shape, please check results'
                )
                warnings.warn(msg)
            geo = shp_geo.iloc[0]
            
        else:
            lon0 = self.config['LonMin']
            lat0 = self.config['LatMin']
            lon1 = self.config['LonMax']
            lat1 = self.config['LatMax']
            geo = Polygon([
                [lon0, lat0],
                [lon1, lat0],
                [lon1, lat1],
                [lon0, lat1],
                [lon0, lat0]
            ])
            
        return geo
    
    def filter_points(self, points):
        '''
        Return boolean dataframe, True for points
        in ROI, false for points outside ROI
         
        Arguments
        ---------
        points: GeoDataFrame
         
        Returns: Series with booleans, index
            matches input GeoDataFrame
        '''
        if points.shape[0] == 0:
            print('No CarbonMapper points found')
            return None
            
        inout = points.within(self.geo)
        
        if not inout.any():
            print('No CarbonMapper points in ROI')
            return None
        
        else:
            return inout

    
class SRONPlumes:
    pass


class SRONPlumes:
    pass



class CarbonMapper:
    '''
    CarbonMapper point source data retrieved
    from the API. Converts data to GeoDataFrame
    with point geometries.

    Parameters
    ----------
    config: dict of IMI config.yml contents

    Use
    ---
    cmapper = CarbonMapper(config)

    Access the GeoDataFrame with data:
    
    cmapper.gdf
    

    '''
    def __init__(self, config):
        self.config = config
        self.gdf = self._get_data()
        self.myname = 'CarbonMapper'

    def _get_data(self):
        '''
        api call to carbonmapper
        '''

        # build query
        bbox = [
            self.config['LonMin'],
            self.config['LatMin'],
            self.config['LonMax'],
            self.config['LatMax'],
        ]
        q_opts = [
            'bbox='+'&bbox='.join([str(i) for i in bbox]),
            'plume_gas=CH4',
            'status=published',
            'eps=100', # minimal distance [m], plumes within this distance will be grouped
            'minpoints=1', # min num plumes in a group
        ]
        base_url = 'https://api.carbonmapper.org'
        endpoint = '/api/v1/catalog/sources.geojson'
        q_params = '?' + '&'.join(q_opts)
        url = base_url + endpoint + q_params

        # get data
        print("Fetching plumes from CarbonMapper...")
        response = requests.get(url)

        # Raise an exception if the API call returns an HTTP error status
        response.raise_for_status()  
        
        # Process the API response
        data = response.json()

        # geodataframe from geojson
        gdf = gpd.GeoDataFrame.from_features(data['features'])
        
        # keep only vars we want (subject to change...)
        keepv = ['detection_count','observation_count','emission','geometry']
        gdf = gdf.rename({
            'detection_date_count': 'detection_count',
            'observation_date_count': 'observation_count',
            'emission_auto': 'emission'
        }, axis=1)[keepv]

        # is a geodataframe with points as geometry
        return gdf



def read_point_source_csv(path):
    """
    Description:
        Read lat, lon coordinates either from a csv file and return as list
    arguments:
        path   String : path to a csv file
    Returns:                [[]] : list of [lat, lon] coordinates of floats
    """
    if not path.endswith(".csv"):
        raise Exception(
            "ForcedNativeResolutionElements expects either a .csv file or a list of lists."
        )
    coords_df = pd.read_csv(path)

    # check if lat and lon columns are present
    if not ("lat" in coords_df.columns and "lon" in coords_df.columns):
        raise Exception(
            "lat or lon columns are not present in the csv file."
            + " csv file must have lat and lon in header using lowercase."
        )
    # select lat and lon columns and convert to list of lists
    return coords_df[["lat", "lon"]].values.tolist()


def get_point_source_coordinates(config):
    """
    Description:
        Read point point sources from config file, csv file, or user-specified
        external point source dataset
    arguments:
        config        {} : dictionary of config variables
    Returns:        [[]] : list of [lat, lon] coordinates of floats
    """
    coords = []
    # first we attempt to read direct, user inputted point sources
    if "ForcedNativeResolutionElements" in config.keys():
        coord_var = config["ForcedNativeResolutionElements"]
        # handle path to csv file containg coordinates
        if isinstance(coord_var, str):
            coords = read_point_source_csv(coord_var)
        # handle list of lists
        elif isinstance(coord_var, list):
            coords = coord_var
        else:
            # Variable is neither a string nor a list
            print(
                "Warning: No ForcedNativeResolutionElements specified or invalid format."
            )

    # then we read point sources from external datasources
    plumes = []
    got_plumes = False
    if "PointSourceDatasets" in config.keys():
        if "SRON" in config["PointSourceDatasets"]:
            print("Fetching plumes from SRON database...")
            plumes += SRON_plumes(config)
            got_plumes = True
        if "CarbonMapper" in config["PointSourceDatasets"]:
            print("Fetching plumes from CarbonMapper database...")
            # append CarbonMapper
            gf = GeoFilter(config)
            cmapper = CarbonMapper(config)
            ps = PointSources(gf, cmapper)
            plumes += ps.get_gridded_coords()
            got_plumes = True
            inout = gf.filter_points(cmapper.gdf)
            cmapper.gdf[inout].to_csv('/n/holylfs05/LABS/jacob_lab/Users/jeast/proj/imi/point_source_sv/carbonmapper_sources_testing.csv')
    
        if not got_plumes:
            print(
                'No valid external point source datasets specified. Valid values are: "SRON", "CarbonMapper"'
            )

    # append point sources
    coords = coords + plumes
    
    return coords


def get_plumes(month, year):
    """
    Description:
        Scrapes the SRON database for weekly methane plumes, saving each week's data as csv files
    arguments:
        month        String : the month (number) for which to select the plumes
        year         String : the year (number) for which to select the plumes
    Returns:                pd.Dataframe() : pandas dataframe with all of the plumes detected for that month
    """
    sron_url = "https://earth.sron.nl/wp-content/uploads/"  # URL of the SRON database for weekly methane plumes
    url = sron_url + year + "/" + month.zfill(2)
    write_dir = "SRON_plumes"
    if not os.path.exists(write_dir):
        os.makedirs(write_dir)
    response = requests.get(url)
    parser = BeautifulSoup(response.content, "html.parser")
    plume = pd.DataFrame()
    for link in parser.find_all("a"):
        if ".csv" in link.get("href") and "SRON_Weekly_Methane_Plumes" in link.get(
            "href"
        ):  # filters through all .csv files containing "SRON_Weekly_Methane_Plumes"
            csvUrl = url + "/" + link.get("href")
            dates = csvUrl.split("_v")[1]
            try:
                rcsv = requests.get(csvUrl, allow_redirects=True)
                file_path = f"{write_dir}/SRON_{dates}"
                # downloads all of the plumes from that week in a CSV file in the current directory
                open(file_path, "wb").write(rcsv.content)
                # reads from the csv file into a pandas dataframe
                df = pd.read_csv(file_path)
                plume = pd.concat([plume, df], ignore_index=True)
            except Exception as err:
                print(
                    f"Warning: Unable to access data for csv file at {csvUrl}. "
                    + "The file may not exist or there may be a connection problem."
                    + f"\nError message: {err}"
                )
    return plume


def shapefile_filter(plumes, shapefile_path):
    """
    Description:
        Removes any plumes (coordinates) that are not within a given shapefile
    arguments:
        plumes   pd.Dataframe() : a pandas dataframe with columns 'lon' and 'lat'
        shapefile_path      String : a string with the path to the shapefile of the ROI
    Returns:                pd.Dataframe() : pandas dataframe only containing coordinates within the shapefile
    """
    shapefile = gpd.read_file(shapefile_path)
    for lon, lat in zip(plumes["lon"], plumes["lat"]):
        point = Point(lon, lat)
        is_within = shapefile.contains(point)
        # checks if it is within any of the polygons if multiple
        is_within_any = is_within.any()
        if not is_within_any:
            plumes = plumes[
                (plumes["lon"] != float(lon)) | (plumes["lat"] != float(lat))
            ]

    return plumes


def rectangular_filter(plumes, LatMax, LatMin, LonMax, LonMin):
    """
    Description:
        Removes any plumes (coordinate) not within a given set of coordinates
    arguments:
        plumes   pd.Dataframe() : a pandas dataframe with columns 'lon' and 'lat'
        LatMax          float : a float indicating the maximum latitude in the ROI
        LatMin          float : a float indicating the minimum latitude in the ROI
        LonMax          float : a float indicating the maximum longitude in the ROI
        LonMin          float : a float indicating the minimum longitude in the ROI
    Returns:                pd.Dataframe() : pandas dataframe
    """
    inLat = (plumes["lat"] > LatMin) & (plumes["lat"] < LatMax)
    inLon = (plumes["lon"] > LonMin) & (plumes["lon"] < LonMax)
    filtered_plumes = plumes[inLat & inLon]
    return filtered_plumes


def SRON_plumes(config):
    """
    Description:
        Selects all the recorded methane plumes on the SRON database for the selected time frame and region
    arguments:
        config              parsed YAML file
    Returns:                [[]] : list of [lat, lon] coordinates of floats
    """

    # variables from config file, specifying time period and region.
    plumes = pd.DataFrame()
    shapefile_path = config["ShapeFile"]
    startDate = datetime.datetime.strptime(str(config["StartDate"]), "%Y%m%d")
    endDate = datetime.datetime.strptime(str(config["EndDate"]), "%Y%m%d")
    if endDate.year < 2023:  # SRON plumes are only available beginning in 2023
        return []
    custom_vectorfile = not config["CreateAutomaticRectilinearStateVectorFile"]
    LatMax = config["LatMax"]
    LatMin = config["LatMin"]
    LonMax = config["LonMax"]
    LonMin = config["LonMin"]
    currentDate = startDate

    # calls the get_plumes function for every month in the selected time frame
    while currentDate <= endDate:
        if currentDate.year >= 2023:  # SRON plumes are only available beginning in 2023
            p = get_plumes(str(currentDate.month), str(currentDate.year))
            plumes = pd.concat([plumes, pd.DataFrame(p)], ignore_index=True)
        currentDate = currentDate + relativedelta(months=1)
    
    # if no plumes found in the time frame, returns an empty list
    if plumes.empty:
        return []
        
    # filters through the dataset to remove any plumes outside the ROI
    if custom_vectorfile:
        plumes = shapefile_filter(
            plumes, shapefile_path
        )  # calls function to filter through coordinates found in shapefile
    else:
        plumes = rectangular_filter(plumes, LatMax, LatMin, LonMax, LonMin)

    #plumes.to_csv('plumes_testing.csv')
    plumes_list = plumes[["lat", "lon"]].values.tolist()
    return plumes_list
