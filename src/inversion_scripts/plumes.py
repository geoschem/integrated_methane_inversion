import os
import io
import sys
import warnings
import requests
import pandas as pd
import geopandas as gpd
import numpy as np
import xarray as xr
from bs4 import BeautifulSoup
from shapely.geometry import Point
from shapely.geometry import Polygon
from functools import partial


print = partial(print, flush=True)

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
        include CarbonMapper, SRON, IMEO. Others
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
        self.gdf = self._merge_dataframes()

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
    
    def _merge_dataframes(self):
        merged_df = pd.concat([s.gdf.copy() for s in self.datasources])
        merged_df['time'] = pd.to_datetime(merged_df['time'], utc=True)
        return merged_df
    
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
        
        if inout is None:
            return None
        elif inout.sum() > 0:
            gdf_inroi = datasource.gdf[inout].copy()
        else:
            msg = (
                'Unable to detect point sources, '
                'continuing without plume data.'
            )
            warnings.warn(msg)
            return None

        
        # 3. get I,J indices
        lon_dist = gdf_inroi.geometry.x.values[:,None] - self.geofilter.lons[None,:]
        lat_dist = gdf_inroi.geometry.y.values[:,None] - self.geofilter.lats[None,:]
        ilon = np.argmin( np.abs(lon_dist), 1)
        ilat = np.argmin( np.abs(lat_dist), 1)
        gdf_inroi['I'] = ilon
        gdf_inroi['J'] = ilat
        
        # 4. average up to grid, sum count, etc
        gdf_inroi['non_detect'] = gdf_inroi['emission_rate'].isna()
        gdf_inroi['plume_count'] = ~gdf_inroi['non_detect']

        keepv = ['plume_count','emission_rate', 'non_detect', 'I', 'J']
        dfgb = gdf_inroi[keepv].groupby(['I','J'])

        gdf_grid = dfgb.sum()[['plume_count']]
        gdf_grid['non_detect'] = dfgb.sum()[['non_detect']]
        gdf_grid['emission_rate'] = dfgb.agg(np.nanmean)[['emission_rate']]

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
        '''
        For all data source (e.g. IMEO, CarbonMapper, SRON)
        provided, grid point source information to the GEOS-Chem
        grid by calling _grid_datasource for each datasource.
        
        returns: xarray.dataset of gridded point source
                 information if it exists, else none. Dataset
                 has coordinate dimension "observer" which
                 is the name of the datasource (e.g. "IMEO")
        '''
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
        from all datasets. Non-nan points only.
         
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
        
    
    def get_gridded_coords(self, emission_rate_filter = 0, plume_count_filter = 0):
        '''
        Get y,x (lat,lon) coordinates for point sources
        from all datasets

        Arguments
        ---------
        emission_rate_filter: kg/hr filter, grid cells with
            mean emission less than this are not included.
            0 means all plumes included.

        plume_count_filter: grid cells with plume_count less
            than this are not included. 0 means no filtering
            applied and only emission_rate_filter included.
         
        Returns
        -------
        list of [y,x] coordinates in [lon,lat]

        '''
        if self.grid_ds is None:
            return []

        else:

            if plume_count_filter > 0:

                print(f'Point sources: {emission_rate_filter = }')
                print(f'Point sources: {plume_count_filter = }')

                criteria = lambda x: (
                    (x.emission_rate > emission_rate_filter) | 
                    (x.plume_count > plume_count_filter)
                )

            else:

                print(f'Point sources: {emission_rate_filter = }')
                print(f'Point sources: No plume count filter applied')

                criteria = lambda x: x.emission_rate > emission_rate_filter

            dftmp = (
                self.grid_ds
                .where(criteria)
                .emission_rate
                .mean('observer')
                .to_dataframe()
                .reset_index()
            )
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
    use_shapefile: bool, whether to define the ROI
        using the shapefile provided in the config
        file.

    Use
    ---
    gf = GeoFilter(config)

    To check whether some other point geometry 
    is in ROI:
    gf.filter_points()

    '''
    def __init__(self, config, use_shapefile = False):
        self.config = config
        self.use_shapefile = use_shapefile
        self.geo = self._make_roi_geometry()
        self.svds = self._get_state_vector_file()
        self.lons = self.svds.lon.values
        self.lats = self.svds.lat.values

    def _get_state_vector_file(self):
        
        # infer lat lons from state vector file
        # this file should already exist
        # expand any environment variables
        svf = f'{self.config["OutputPath"]}/{self.config["RunName"]}/StateVector.nc'
        svf = os.path.expandvars(svf)
        
        try:
            svds = xr.load_dataset(svf)
        except FileNotFoundError:
            msg = (
                f'State vector file {svf} '
                'cannot be found. Re-run IMI '
                'with "RunSetup" True.'
            )
            sys.exit(msg)
        
        return svds

        
        
    def _make_roi_geometry(self):
        
        # optionally use the shapefile
        # provided in the config file
        if self.use_shapefile:
            shapefile_path = self.config["ShapeFile"]
            shp_geo = gpd.read_file(shapefile_path).geometry
            if shp_geo.shape[0] > 1:
                msg = (
                    'Shapefile has >1 shape, only using point sources in '
                    'first shape, please check results'
                )
                warnings.warn(msg)
            geo = shp_geo.iloc[0]

            
        # else define ROI based on bounds
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
            print('No plumes found.')
            return None
            
        inout = points.within(self.geo)
        
        if not inout.any():
            print('No plumes in ROI')
            return None
        
        else:
            return inout




class PlumeObserver:

    def __init__(self, myname, config, usecached, cache=None):
        self.usecached = usecached
        self.myname = myname
        self.cache = cache
        self.config = config
        
        if self.usecached:
            msg = (
                f'Using cached plume data for {self.myname}'
            )
            print(msg)
            if self.cache is not None:
                infile = self.cache
            else:
                basedir = f'{self.config["OutputPath"]}/{self.config["RunName"]}'
                infile = f'{basedir}/{self.myname}_plumes/plumes.geojson'
            try:
                gdf = gpd.read_file(infile)
                gdf['time'] = pd.to_datetime(gdf['time'])
                self.gdf = gdf
            except:
                msg = (
                    f'Cannot open cached file at {infile}, fetching new data.'
                )
                warnings.warn(msg)
                self.gdf = self._get_data()
        
        else:
            self.gdf = self._get_data()

        if self.gdf is not None:
            if self.gdf.shape[0] > 0:
        
                # cache all data
                self._cache_data()
                
                # filter data according to 
                # options from subclass
                self._filter_data()
        
    def _filter_data(self):
        raise NotImplementedError
        
    def _get_data(self):
        raise NotImplementedError
        
    def _cache_data(self):
        # dir for saving file
        basedir = f'{self.config["OutputPath"]}/{self.config["RunName"]}'
        write_dir = f'{basedir}/{self.myname}_plumes'
        if not os.path.exists(write_dir):
            os.makedirs(write_dir)
        # write file to geojson
        self.gdf.to_file(
            f'{write_dir}/plumes.geojson',
            driver = 'GeoJSON' 
        )


class SRON(PlumeObserver):
    '''
    SRON point source data retrieved
    from the API. Converts data to GeoDataFrame
    with point geometries.

    Parameters
    ----------
    config: dict of IMI config.yml contents
    usecached: bool, whether to use saved data
    cache: None or path to existing file. Will search 
        for saved data in IMI output dir if no
        cache given.

    Use
    ---
    sron = SRON(config, cache=cache_file)

    Access the GeoDataFrame with data:
    
    sron.gdf
    

    '''
    def __init__(self, config, usecached = True, cache = None):
        self.myname = 'SRON'
        self.cache = cache
        super().__init__(self.myname, config, usecached, self.cache)
        
    def _filter_data(self):
        pass

    def _get_data(self):
        
        # get data
        print(f"Fetching plumes from {self.myname}...")
        
        # URL of the SRON database for weekly methane plumes
        sron_url = 'https://ftp.sron.nl/pub/memo/CSVs/'
        
        response = requests.get(sron_url)
        if not response.status_code == 200:
            msg = (
                f"Couldn't access SRON site {sron_url}.  "
                f"failed with error {response.status_code}: "
                f"{response.reason}.\n"
                "SRON plumes not retrieved.  "
                "Check SRON website or remove SRON from "
                "point source list in config file."
            )
            raise Exception(msg)
        
        parser = BeautifulSoup(response.content, "html.parser")
        
        dfplume = pd.DataFrame()
        
        # retrieve csvs, read into dataframe
        dfplume = pd.DataFrame()
        for link in parser.find_all("a"):
            file_stem = link.get('href')
            if file_stem.endswith('.csv'):
                try:
                    # read the raw text, make it a dataframe
                    sron_file = '/'.join([sron_url, file_stem])
                    rcsv = requests.get(sron_file, allow_redirects=True)
                    df = pd.read_csv(io.StringIO(rcsv.text), dtype = {'date':str})
                    dfplume = pd.concat([dfplume, df], ignore_index=True)
                except Exception as err:
                    print(
                        f"Warning: Unable to access data for csv file at {sron_file}. "
                        + "The file may not exist or there may be a connection problem."
                        + f"\nError message: {err}"
                    )

        
        df = dfplume.copy()
        
        # drop nan times
        # drops when there is tropomi error and
        # therefore no data reported
        df = df.loc[~df['time_UTC'].isna(),:]
        
        # SRON changed naming convention for uncertainty
        # columns, so need to capture both
        df['emission_rate_std'] = np.where(
            np.isnan(df['uncertainty_t/h']),
            df['uncertainty'],
            df['uncertainty_t/h']
        )

        gdf = gpd.GeoDataFrame(
            data = df,
            geometry = gpd.points_from_xy(df.lon, df.lat),
            crs = 'EPSG:4326'
        ).rename({
            'source_rate_t/h': 'emission_rate',
        }, axis=1)
        
        gdf['time'] = pd.to_datetime(
            gdf['date'] + 'T' + gdf['time_UTC']
        )
        gdf['instrument'] = 'TROPOMI'
        
        # t/h to kg/hr
        gdf.loc[:,'emission_rate'] *= 1000
        gdf.loc[:,'emission_rate_std'] *= 1000
        
        gdf = gdf[['emission_rate', 'emission_rate_std', 'time', 'instrument', 'geometry']]

        # for SRON, subset to time of interest
        gdf = (
            gdf
            .set_index('time')
            .sort_index()
            .loc[str(self.config['StartDate']):str(self.config['EndDate'])]
            .reset_index()
        )
        gdf['datasource'] = self.myname

        return gdf
    


class CarbonMapper(PlumeObserver):
    '''
    CarbonMapper point source data retrieved
    from the API. Converts data to GeoDataFrame
    with point geometries.

    Parameters
    ----------
    config: dict of IMI config.yml contents
    usecached: bool, whether to use saved data
    cache: None or path to existing file. Will search 
        for saved data in IMI output dir if no
        cache given.
    sat_only: bool, whether to filter out aircraft
        data. Option for CarbonMapper only.

    Use
    ---
    cmapper = CarbonMapper(config)

    Access the GeoDataFrame with data:
    
    cmapper.gdf
    

    '''
    def __init__(self, config, usecached = False, cache=None, sat_only = False):
        self.myname = 'CarbonMapper'
        self.cache = cache
        self.sat_only = sat_only
        super().__init__(self.myname, config, usecached, self.cache)
        
    def _filter_data(self):
        # optionally filter out aircraft data
        if self.sat_only:
            print('Dropping aircraft data from CarbonMapper dataset')
            aircraft_instruments = ['GAO', 'ang']
            self.gdf = self.gdf[~self.gdf.instrument.isin(aircraft_instruments)].copy()
        

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
        if None in bbox:
            # then use statevector file bounds
            ds = xr.load_dataset(self.config['StateVectorFile'])
            bbox = [
                int(ds.lon.min().values),
                int(ds.lat.min().values),
                int(ds.lon.max().values),
                int(ds.lat.max().values)
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
            try:
                response.raise_for_status() 
            except Exception as e:
                print(f'CarbonMapper API error: {e}')
                print('Continuing WITHOUT CarbonMapper plumes.')
                return None
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
        gdf['datasource'] = self.myname
        
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
    usecached: bool, whether to use saved data
    cache: None or path to existing file. Will search 
        for saved data in IMI output dir if no
        cache given.

    Use
    ---
    imeo = IMEO(config)

    Access the GeoDataFrame with data:
    
    imeo.gdf
    

    '''
    def __init__(self, config, usecached = False, cache = None):
        self.myname = 'IMEO'
        self.cache = cache
        super().__init__(self.myname, config, usecached, self.cache)
        
    def _filter_data(self):
        pass

    def _get_data(self):
        '''
        download data from imeo
        '''
       
        # dir for saving file
        basedir = f'{self.config["OutputPath"]}/{self.config["RunName"]}'
        write_dir = f'{basedir}/{self.myname}_plumes'
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
            keepv = ['emission_rate','emission_rate_std','instrument','time','geometry']
            gdf = gdf.rename({
                'ch4_fluxrate': 'emission_rate',
                'ch4_fluxrate_std': 'emission_rate_std',
                'tile_date': 'time',
                'satellite': 'instrument'
            }, axis=1)[keepv]
            gdf['time'] = pd.to_datetime(gdf['time'], format='ISO8601')
            gdf['datasource'] = self.myname
            
            return gdf
        
        except Exception as err:
            msg = (
                f"Unable to access data at {url}. "
                + "The file may not exist or there may be a connection problem."
                + f"\nError message: {err}"
            )
            warnings.warn(msg)
            
        
