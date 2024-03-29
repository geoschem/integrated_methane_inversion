import os
import warnings
import requests
import pandas as pd
import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import Point
from shapely.geometry import Polygon


class PointSources:
    '''
    Point source plume and emissions data filtered to ROI.

    Multiple point source datasets can be included and data
    will be concatenated. Current options include:
        - SRON plume data
        - CarbonMapper plume data
    
    
    Datasrouces must have for each point
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
        self.grid_ds = self._grid_all_data()

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
            gdf_filtered = gdf_combined[~gdf_combined.emission_rate.isna()].copy()
            return gdf_filtered.geometry
        else:
            return None
    
    
    def _grid_datasource(self, datasource):
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
        gdf_filtered = gdf_inroi[~gdf_inroi.emission_rate.isna()].copy()
        
        # 3. get I,J indices
        lon_dist = gdf_filtered.geometry.x.values[:,None] - self.geofilter.lons[None,:]
        lat_dist = gdf_filtered.geometry.y.values[:,None] - self.geofilter.lats[None,:]
        ilon = np.argmin( np.abs(lon_dist), 1)
        ilat = np.argmin( np.abs(lat_dist), 1)
        gdf_filtered['I'] = ilon
        gdf_filtered['J'] = ilat
        
        # 4. average up to grid, sum count, etc
        dfgb = gdf_filtered.groupby(['I','J'])
        gdf_grid = (
            dfgb
            .count()[['emission_rate']]
            .rename({'emission_rate': 'plume_count'}, axis=1)
        )
        gdf_grid['emission_rate'] = dfgb.mean()[['emission_rate']]
        #gdf_grid['persistence'] = gdf_grid.detection_count / gdf_grid.observation_count

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
    
    
    def _grid_all_data(self):
        d_list = []
        d_names = []
        for d in self.datasources:
            ds_i = self._grid_datasource(d)
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
            dftmp = self.grid_ds.emission_rate.to_dataframe().reset_index()
            coords = (
                dftmp[~dftmp['emission_rate'].isna()][['lat','lon']]
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
        svf = f'{self.config["OutputPath"]}/{self.config["RunName"]}/StateVector.nc'
        
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




class PlumeObserver:

    def __init__(self, usecached, myname):
        self.usecached = usecached
        self.myname = myname
        
        if self.usecached:
            msg = (
                f'Using cached plume data for {self.myname}'
            )
            print(msg)
            infile = f'{self.myname}_plumes/plumes.geojson'
            try:
                self.gdf = gpd.read_file(infile)
            except:
                msg = (
                    'No cached file, fetching new data.'
                )
                warnings.warn(msg)
                self.gdf = self._get_data()
        
        else:
            self.gdf = self._get_data()
        
        self._cache_data()
        
    def _get_data(self):
        raise NotImplementedError
        
    def _cache_data(self):
        # dir for saving file
        write_dir = f'{self.myname}_plumes'
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        # write file to geojson
        self.gdf.to_file(
            f'{write_dir}/plumes.geojson',
            driver = 'GeoJSON' 
        )




class CarbonMapper(PlumeObserver):
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
    def __init__(self, config, usecached = False):
        self.myname = 'CarbonMapper'
        self.config = config
        super().__init__(usecached, self.myname)

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
        q_opts = lambda dlimit, offset: [
            'bbox='+'&bbox='.join([str(i) for i in bbox]),
            'plume_gas=CH4',
            f'limit={dlimit}',
            f'offset={offset}',
        ]
        base_url = 'https://api.carbonmapper.org'
        endpoint = '/api/v1/catalog/plumes/annotated'
        
        # get data
        print("Fetching plumes from CarbonMapper...")
        
        # first check how much data
        q_params = '?' + '&'.join(q_opts(1,0))
        url = base_url + endpoint + q_params
        response = requests.get(url)
        # Raise an exception if the API call returns an HTTP error status
        response.raise_for_status()  
        # Process the API response
        data = response.json()
        
        # total num data
        ndata = data['bbox_count']
        limit = 1000
        npages = ndata // limit + 1
        print(f'Found {ndata} Carbon Mapper plumes to download.')
        
        # loop through api pages
        rawdat = []
        print('Downloading',end='')
        for ioffset in range(npages):
            offset = ioffset * limit
            print('.',end='')
            q_params = '?' + '&'.join(q_opts(limit,offset))
            url = base_url + endpoint + q_params
            response = requests.get(url)
            # Raise an exception if the API call returns an HTTP error status
            response.raise_for_status()  
            # Process the API response
            data = response.json()
            rawdat += data['items']
        print()

        # geodataframe from json data
        df_raw = pd.json_normalize(rawdat)
        lonlats = np.array(df_raw['geometry_json.coordinates'].tolist())
        gdf = gpd.GeoDataFrame(
            df_raw,
            geometry = gpd.points_from_xy(lonlats[:,0], lonlats[:,1]),
            crs = "EPSG:4326"
        )
        
        # keep only vars we want (subject to change...)
        keepv = [
            'emission_rate','emission_rate_std',
            'time', 'instrument', 'platform',
            'geometry'
        ]
        gdf = gdf.rename({
            'scene_timestamp': 'time',
            'emission_auto': 'emission_rate',
            'emission_uncertainty_auto': 'emission_rate_std'
        }, axis=1)[keepv]
        
        gdf['time'] = pd.to_datetime(gdf['time'])

        # is a geodataframe with points as geometry
        return gdf


class IMEO(PlumeObserver):
    '''
    Methane plume data from IMEO portal.
    Converts data to GeoDataFrame with
    point geometries.

    Parameters
    ----------
    config: dict of IMI config.yml contents

    Use
    ---
    imeo = IMEO(config)

    Access the GeoDataFrame with data:
    
    imeo.gdf
    

    '''
    def __init__(self, config, usecached = False):
        self.myname = 'IMEO'
        self.config = config
        super().__init__(usecached, self.myname)

    def _get_data(self):
        '''
        download data from imeo
        '''
       
        # dir for saving file
        write_dir = "IMEO_plumes"
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
    
        # build url for geojson
        url = (
            'https://unepazeconomyadlsstorage.blob.core.windows.net/'
            'public/unep_methanedata_detected_plumes.geojson'
        )
            
        print("Fetching plumes from IMEO...")
        try:
            file_path = f"{write_dir}/{url.split('/')[-1]}"
            rgjson = requests.get(url, allow_redirects=True)
            open(file_path, "wb").write(rgjson.content)
            
            gdf = gpd.read_file(file_path)
            points = gpd.points_from_xy(gdf.lon, gdf.lat, crs = "EPSG:4326")
            gdf = (
                gdf
                .rename({'geometry': 'plume_geometry'}, axis=1)
                .set_geometry(points)
            )
            
            # data are in kg/hr
            # keep only vars we want (subject to change...)
            keepv = ['emission_rate','emission_rate_std','satellite','time','geometry']
            gdf = gdf.rename({
                'ch4_fluxrate': 'emission_rate',
                'ch4_fluxrate_std': 'emission_rate_std',
                'tile_date': 'time'
            }, axis=1)[keepv]
            gdf['time'] = pd.to_datetime(gdf['time'])
            
            return gdf
        
        except Exception as err:
            msg = (
                f"Unable to access data at {url}. "
                + "The file may not exist or there may be a connection problem."
                + f"\nError message: {err}"
            )
            warnings.warn(msg)
            
        

class SRON:
    pass
