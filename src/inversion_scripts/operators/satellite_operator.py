import os
from pathlib import Path
import importlib.util
import numpy as np
import xarray as xr
import pandas as pd
import datetime
import pygeohash as pgh
from shapely.geometry import Polygon
from src.inversion_scripts.utils import (
    read_and_filter_satellite,
    mixing_ratio_conv_factor,
    get_strdate,
    check_is_OH_element,
    check_is_BC_element,
)
from src.inversion_scripts.operators.operator_utilities import (
    get_gc_lat_lon,
    read_all_geoschem,
    merge_pressure_grids,
    remap,
    remapping_weights,
    get_gridcell_list,
    nearest_loc,
    get_overlap_area_CSgrid,
)
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="xarray")

def get_goopy_config_path():
    spec = importlib.util.find_spec("GOOPy")
    if spec is not None and spec.submodule_search_locations:
        config_path = Path(spec.submodule_search_locations[0]) / "config.yaml"
        if config_path.exists():
            return config_path

    repo_config_path = Path(__file__).resolve().parents[3] / "GOOPy" / "config.yaml"
    if repo_config_path.exists():
        return repo_config_path

    raise FileNotFoundError(
        "Could not locate GOOPy/config.yaml from the importable GOOPy package "
        "or the repository root. Make sure the IMI repository root is on PYTHONPATH."
    )

def obs_to_xarray_dataset(obs_mapped_to_gc, species, filename, config):
    """
    Convert superobservations from structured numpy array to xarray Dataset.
    
    Arguments
        obs_mapped_to_gc [numpy.ndarray] : Structured array of superobservations
        species          [str]           : Species name (CH4 or CO2)
        filename         [str]           : Original satellite filename (for metadata)
        config           [dict]          : Configuration dictionary
    
    Returns
        ds               [xr.Dataset]    : xarray Dataset with superobservations
    """
    n_obs = len(obs_mapped_to_gc)
    
    # Get dimensions from the first observation
    n_lev_p = len(obs_mapped_to_gc["p_sat"][0])
    n_lev_layer = len(obs_mapped_to_gc["layer"][0])
    
    # Create dimensions and coordinates
    coords = {
        'nobs': np.arange(n_obs),
        'layer': obs_mapped_to_gc["layer"][0],  # layer index from first observation
        'corner': np.arange(4),
    }
    
    # Map species name for netcdf variable name
    if species.upper() == 'CH4':
        satellite_column_name = 'methane_mixing_ratio_blended'
        prior_name = 'methane_profile_apriori'
    elif species.upper() == 'CO2':
        satellite_column_name = 'xco2_mixing_ratio_blended'
        prior_name = 'co2_profile_apriori'
    else:
        satellite_column_name = f'{species}_mixing_ratio_blended'
        prior_name = f'{species}_profile_apriori'

    time_utc = pd.to_datetime(
        obs_mapped_to_gc['time'],
        format='%Y%m%d_%H',
    ).strftime('%Y-%m-%dT%H:%M:%S').to_numpy(dtype=str)

    # IMI stores these layer-resolved fields surface-to-TOA. GOOPy's
    # TROPOMI_blended parser constructs pressure edges in the native file order
    # and then flips N_EDGES/N_CENTERS to get descending pressure. Write the
    # center fields in that native order so GOOPy's flip restores IMI ordering.
    avkern_for_goopy = np.array(
        [obs_mapped_to_gc['avkern'][i] for i in range(n_obs)],
        dtype=np.float32,
    )[:, ::-1]
    apriori_for_goopy = np.array(
        [obs_mapped_to_gc['apriori'][i] for i in range(n_obs)],
        dtype=np.float32,
    )[:, ::-1]
    dryair_for_goopy = np.array(
        [obs_mapped_to_gc['dry_air_subcolumns'][i] for i in range(n_obs)],
        dtype=np.float32,
    )[:, ::-1]
    
    # Create data variables matching TROPOMI_blended format from GOOPy config.
    # GOOPy uses latitude/longitude to choose the model grid cell. These
    # superobservations have already been assigned to GEOS-Chem cells by IMI,
    # so use the GC cell centers here rather than the weighted satellite
    # footprint centers.
    data_vars = {
        'latitude': (['nobs'], obs_mapped_to_gc['lat'].astype(np.float32)),
        'longitude': (['nobs'], obs_mapped_to_gc['lon'].astype(np.float32)),
        'satellite_latitude': (['nobs'], obs_mapped_to_gc['lat_sat'].astype(np.float32)),
        'satellite_longitude': (['nobs'], obs_mapped_to_gc['lon_sat'].astype(np.float32)),
        'time_utc': (['nobs'], time_utc),
        satellite_column_name: (['nobs'], obs_mapped_to_gc[species].astype(np.float32)),
        'surface_pressure': (['nobs'], obs_mapped_to_gc['surface_pressure'].astype(np.float32)),
        'surface_albedo_NIR': (['nobs'], obs_mapped_to_gc['nir_albedo'].astype(np.float32)),
        'surface_albedo_SWIR': (['nobs'], obs_mapped_to_gc['swir_albedo'].astype(np.float32)),
        'column_averaging_kernel': (
            ['nobs', 'layer'],
            avkern_for_goopy
        ),
        prior_name: (
            ['nobs', 'layer'],
            apriori_for_goopy
        ),
        'pressure_interval': (
            ['nobs'],
            (obs_mapped_to_gc['p_sat'][:, 0].astype(np.float32) - obs_mapped_to_gc['p_sat'][:, 1].astype(np.float32)) * 100.0
            if n_lev_p > 1 else np.ones(n_obs, dtype=np.float32)
        ), # NOTE: we multiply by 100 to convert from hPa to Pa, which is the unit expected by GOOPy and the TROPOMI_blended format
        'dry_air_subcolumns': (
            ['nobs', 'layer'],
            dryair_for_goopy
        ),
        'observation_count': (['nobs'], obs_mapped_to_gc['observation_count'].astype(np.float32)),
        'qa_value': (['nobs'], np.ones(n_obs, dtype=np.float32)),  # Set to 1.0 for all valid observations
        'gc_cell_lat': (['nobs'], obs_mapped_to_gc['lat'].astype(np.float32)),
        'gc_cell_lon': (['nobs'], obs_mapped_to_gc['lon'].astype(np.float32)),
        'gc_cell_index_lat': (['nobs'], obs_mapped_to_gc['jGC'].astype(np.int32)),
        'gc_cell_index_lon': (['nobs'], obs_mapped_to_gc['iGC'].astype(np.int32)),
    }
    
    # Create the xarray Dataset
    ds = xr.Dataset(data_vars, coords=coords)
    
    # Add attributes for metadata
    ds.attrs['title'] = 'Superobservations: satellite observations averaged to GEOS-Chem grid cells'
    ds.attrs['source_file'] = os.path.basename(filename)
    ds.attrs['species'] = species
    ds.attrs['history'] = f'Created by satellite_operator.py at {pd.Timestamp.now().isoformat()}'
    
    return ds


def save_superobservations(ds, filename, output_path):
    """
    Save superobservations xarray Dataset to netcdf file.
    
    Arguments
        ds          [xr.Dataset] : xarray Dataset of superobservations
        filename    [str]        : Original satellite filename (for output filename)
        output_path [str]        : Directory path to save the netcdf file
    
    Returns
        output_file [str]        : Path to the saved netcdf file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)
    
    # Create output filename based on input filename
    base_name = os.path.basename(filename)
    # Replace or append to indicate this is superobservations
    if base_name.endswith('.nc'):
        output_file = os.path.join(output_path, base_name.replace('.nc', '_superobservations.nc'))
    else:
        output_file = os.path.join(output_path, base_name + '_superobservations.nc')
    
    # Save to netcdf with compression
    ds.to_netcdf(output_file)
    # ds.to_netcdf(output_file, encoding={
    #     var: {'zlib': True, 'complevel': 4}
    #     for var in ds.data_vars
    # })
    
    print(f"Saved superobservations to {output_file}")
    return output_file


def apply_operator(operator, params, obs_mapped_to_gc, config):
    """
    Run the chosen operator based on selected instrument

    Arguments
        operator [str]    : Data conversion operator to use
        params   [dict]   : parameters to run the given operator
    Returns
        output   [dict]   : Dictionary with:
                            - obs_GC : GEOS-Chem and satellite column data
                            - satellite columns
                            - GEOS-Chem columns
                            - satellite lat, lon
                            - satellite lat index, lon index
    """
    return goopy_apply_operator(
        operator,
        params["filename"],
        params["species"],
        params["satellite_product"],
        params["satellite_cache"],
        params["n_elements"],
        params["gc_startdate"],
        params["gc_enddate"],
        params["xlim"],
        params["ylim"],
        params["gc_cache"],
        params["period_i"],
        obs_mapped_to_gc,
        config,
        params["use_water_obs"],
    )
    # if operator == "satellite_average":
    #     return apply_average_satellite_operator(
    #         params["filename"],
    #         params["species"],
    #         params["satellite_product"],
    #         params["satellite_cache"],
    #         params["n_elements"],
    #         params["gc_startdate"],
    #         params["gc_enddate"],
    #         params["xlim"],
    #         params["ylim"],
    #         params["gc_cache"],
    #         params["period_i"],
    #         obs_mapped_to_gc,
    #         config,
    #         params["use_water_obs"],
    #     )
    # elif operator == "satellite":
    #     return apply_satellite_operator(
    #         params["filename"],
    #         params["species"],
    #         params["satellite_product"],
    #         params["satellite_cache"],
    #         params["n_elements"],
    #         params["gc_startdate"],
    #         params["gc_enddate"],
    #         params["xlim"],
    #         params["ylim"],
    #         params["gc_cache"],
    #         params["period_i"],
    #         config,
    #         params["use_water_obs"],
    #     )
    # else:
    #     raise ValueError("Error: invalid operator selected.")


def superobservations(
    filename,
    species,
    satellite_product,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    period_i,
    config,
    use_water_obs=False
):
    """
    Compute superobservations for the given satellite file by averaging observations within each grid cell. 
    """
    # Read satellite data
    satellite, sat_ind = read_and_filter_satellite(
        filename, satellite_product, gc_startdate, gc_enddate,
        xlim, ylim, use_water_obs)

    # Number of satellite observations
    n_obs = len(sat_ind[0])
    if n_obs == 0:
        print(f"No satellite observations found in {filename}. Skipping.")
        return None
    print("Found", n_obs, "satellite observations.")

    # Define time threshold (hour 00 after the inversion period)
    date_after_inversion = str(gc_enddate + np.timedelta64(1, "D"))[:10].replace(
        "-", ""
    )
    time_threshold = f"{date_after_inversion}_00"

    # map satellite obs into gridcells and average the observations
    # into each gridcell. Only returns gridcells containing observations
    if config["UseGCHP"]:
        if config['STRETCH_GRID']:
            sf_formatted = f"{config['STRETCH_FACTOR']:.2f}".replace(".", "d")
            target_geohash = pgh.encode(config['TARGET_LAT'], config['TARGET_LON'])
            gridspec_path = f"c{config['CS_RES']}_s{sf_formatted}_t{target_geohash}_gridspec.nc"
        else:
            gridspec_path = f"c{config['CS_RES']}_gridspec.nc"
        GC_shape = (6, config['CS_RES'], config['CS_RES'])
        CSgridDir = f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}/CS_grids"

        obs_mapped_to_gc = average_satellite_observations_to_CSgrid(
            satellite, species, filename, sat_ind, time_threshold, CSgridDir, gridspec_path, GC_shape
        )
    else:
        # get the lat/lons of gc gridcells
        gc_lat_lon = get_gc_lat_lon(gc_cache, gc_startdate)
        obs_mapped_to_gc = average_satellite_observations(
            satellite, species, gc_lat_lon, sat_ind, time_threshold
        )
        GC_shape = (len(gc_lat_lon['lat']), len(gc_lat_lon['lon']))
    
    # Create xarray dataset from obs_mapped_to_gc
    ds = obs_to_xarray_dataset(obs_mapped_to_gc, species, filename, config)
    
    # Save superobservations to netcdf files
    output_dir = os.path.join(os.path.expandvars(config['OutputPath']), config['RunName'], f'{filename}_superobservations')
    save_superobservations(ds, filename, output_dir)

    return obs_mapped_to_gc, output_dir

def goopy_apply_operator(
    operator,
    filename,
    species,
    satellite_product,
    satellite_cache,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    period_i,
    obs_mapped_to_gc,
    config,
    use_water_obs=False
):
    import sys
    import yaml
    import tempfile
    import glob
    import importlib
    
    # Read the full GOOPy config
    goopy_config_path = get_goopy_config_path()
    with open(goopy_config_path, 'r') as f:
        goopy_config = yaml.safe_load(f)
    goopy_package_path = goopy_config_path.parent
    goopy_repo_path = goopy_package_path.parent
    
    # Update LOCAL_SETTINGS with the values needed for this run
    save_dir = f'{gc_cache}/../goopy_output'
    if operator == "satellite_average":
        satellite_files = sorted(glob.glob(os.path.join(satellite_cache, '*.nc')))
        if len(satellite_files) != 1:
            raise ValueError(
                f"Expected exactly one superobservation file in {satellite_cache}, "
                f"found {len(satellite_files)}."
            )
        goopy_obs_file = satellite_files[0]
    elif operator == "satellite":
        goopy_obs_file = os.path.join(satellite_cache, os.path.basename(filename))
        if not os.path.exists(goopy_obs_file):
            raise FileNotFoundError(
                f"Expected satellite file {goopy_obs_file} does not exist."
            )
    else:
        raise ValueError(f"Error: invalid operator selected: {operator}")

    goopy_config['LOCAL_SETTINGS'].update({
        'SAVE_INTERPOLATION': 'False',
        'SATELLITE_NAME': 'TROPOMI_blended',  # TODO: make this not TROPOMI-specific
        'OBS_DIR': os.path.dirname(goopy_obs_file),
        'OBS_FILE_FORMAT': os.path.basename(goopy_obs_file),
        'MODEL_LEVEL_EDGE_DIR': gc_cache,
        # 'LEVEL_EDGE_FILE_FORMAT': 'GEOSChem.LevelEdgeDiags.*.nc4',
        'LEVEL_EDGE_FILE_FORMAT': 'GEOSChem.StateMetLevEdge.*.nc4',
        'MODEL_CONCENTRATION_DIR': gc_cache,
        'CONCENTRATION_FILE_FORMAT': 'GEOSChem.SpeciesConc.*.nc4',
        'SAVE_DIR': save_dir,
    })
    
    # Write to a temporary file
    temp_config = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
    yaml.dump(goopy_config, temp_config)
    temp_config.close()
    
    # Set the config file path for GOOPy and import
    sys.argv[1] = temp_config.name
    for import_path in (goopy_repo_path, goopy_package_path):
        import_path = str(import_path)
        if import_path not in sys.path:
            sys.path.insert(0, import_path)
    GOOPy_main = importlib.import_module("GOOPy.main")
    GOOPy_main = importlib.reload(GOOPy_main)
    
    GOOPy_main.apply_operator(GOOPy_main.config)
    
    # Clean up the temporary file
    os.unlink(temp_config.name)

    goopy_output_file = os.path.join(
        save_dir,
        os.path.splitext(os.path.basename(goopy_obs_file))[0] + '_operator.nc'
    )
    if not os.path.exists(goopy_output_file):
        available_outputs = sorted(
            f for f in os.listdir(save_dir) if f.endswith('_operator.nc')
        ) if os.path.isdir(save_dir) else []
        raise FileNotFoundError(
            f"Expected GOOPy output {goopy_output_file} was not created. "
            f"Available GOOPy outputs: {available_outputs}"
        )

    # debug_row = int(os.environ.get("IMI_DEBUG_ROW", "-1"))

    # print("GOOPy output operator:", operator)
    # print("GOOPy output file:", goopy_output_file)
    # with xr.open_dataset(goopy_output_file) as ds:
    #     if operator == "satellite_average" and debug_row >= 0:
    #         print("\n=== GOOPy AVERAGED OUTPUT ROW DEBUG ===")
    #         print("operator:", operator)
    #         print("file:", goopy_output_file)
    #         print("row:", debug_row)
    #         print("dims:", dict(ds.sizes))
    #         print("data vars:", list(ds.data_vars))

    #         for var_name, da in ds.data_vars.items():
    #             if "N_OBS" in da.dims:
    #                 row = da.isel(N_OBS=debug_row).values
    #             elif "nobs" in da.dims:
    #                 row = da.isel(nobs=debug_row).values
    #             else:
    #                 continue

    #             print(f"\n{var_name}")
    #             print("  dims:", da.dims)
    #             print("  shape:", da.shape)
    #             print("  row value:", row)

    #             arr = np.asarray(row)
    #             if arr.size == 1 and np.isfinite(arr).all():
    #                 print("  scalar as ppb if mol/mol:", float(arr) * 1e9)

    #         print("=== END GOOPy AVERAGED OUTPUT ROW DEBUG ===\n")

    with xr.open_dataset(goopy_output_file) as ds:
        # virtual_satellite = ds['SATELLITE_COLUMN'].values.astype(np.float32)
        virtual_satellite = ds['MODEL_COLUMN_CH4'].values.astype(np.float32)

    if operator == "satellite":
        return format_goopy_satellite_output(
            filename,
            species,
            satellite_product,
            gc_startdate,
            gc_enddate,
            xlim,
            ylim,
            virtual_satellite,
            config,
            use_water_obs,
        )

    # mostly copied from apply_average_satellite_operator() in this file, but with some edits to read in the GOOPy output instead of re-computing the operator components
    n_gridcells = len(obs_mapped_to_gc)
    gc_lat_lon = get_gc_lat_lon(gc_cache, gc_startdate)
    GC_shape = (len(gc_lat_lon['lat']), len(gc_lat_lon['lon']))

    # Initialize array with n_gridcells rows and 5 columns. Columns are
    # satellite species, GEOSChem species, longitude, latitude, observation counts
    obs_GC = np.empty([n_gridcells, 5], dtype=np.float32)
    obs_GC.fill(np.nan)

    if config['UseGCHP']:
        GC_index = np.ravel_multi_index((obs_mapped_to_gc["nfi"],
                                         obs_mapped_to_gc["Ydimi"],
                                         obs_mapped_to_gc["Xdimi"]), GC_shape)
    else:
        GC_index = np.ravel_multi_index((obs_mapped_to_gc["jGC"], #lat
                                         obs_mapped_to_gc["iGC"] # lon
                                         ), GC_shape)

    all_strdate = [gridcell["time"] for gridcell in obs_mapped_to_gc]
    all_strdate = list(set(all_strdate))

    # Read GEOS-Chem data for simulated truth in OSSE simulation
    if config["EnableOSSE"]:
        osse_gc_cache = "./data_geoschem_osse"

        # check if the osse_gc_cache exists
        assert os.path.exists(osse_gc_cache), (
            f"OSSE GEOS-Chem cache directory {osse_gc_cache} does not exist. "
            "Please run the OSSE simulation first."
        )

    if len(virtual_satellite) != n_gridcells:
        raise ValueError(
            f"GOOPy output has {len(virtual_satellite)} satellite columns, "
            f"but {n_gridcells} superobservations were expected."
        )

    # TODO: delete after testing
    debug_row = int(os.environ.get("IMI_DEBUG_ROW", "-1"))
    # print(f"HELLO\n\n\n\n\n\n\n\n\n\n{n_gridcells=}, {len(virtual_satellite)=}, {debug_row=}\n\n\n\n\n\n\n\n\n\n")

    debug_file = os.environ.get("IMI_DEBUG_FILE_SUBSTR")

    if debug_file and debug_file not in os.path.basename(filename):
        pass
    elif 0 <= debug_row < n_gridcells:
        np.set_printoptions(precision=8, suppress=False)

        one = obs_mapped_to_gc[debug_row:debug_row + 1]
        native_virtual = get_virtual_satellite(
            one["time"][0],
            gc_cache,
            one,
            n_elements,
            config,
        )

        print("\n=== IMI -> GOOPy SUPEROBS DEBUG ===")
        print("filename:", os.path.basename(filename))
        print("n_gridcells:", n_gridcells)
        print("native IMI virtual column ppb:", native_virtual[0] * 1e9)
        print("GOOPy virtual column ppb:", virtual_satellite[debug_row] * 1e9)
        print("GOOPy - native IMI ppb:", virtual_satellite[debug_row] * 1e9 - native_virtual[0] * 1e9)

        dryair = obs_mapped_to_gc["dry_air_subcolumns"][debug_row]
        apriori = obs_mapped_to_gc["apriori"][debug_row]
        avkern = obs_mapped_to_gc["avkern"][debug_row]
        p_sat = obs_mapped_to_gc["p_sat"][debug_row]

        print("\n--- IMI arrays in GOOPy-comparable form ---")
        print("IMI pressure edges:", p_sat)
        print("IMI pressure diffs:", np.diff(p_sat))
        print("IMI pressure descending?", np.all(np.diff(p_sat) < 0))

        print("IMI pressure weights:", dryair / np.sum(dryair))
        print("IMI prior profile vmr:", apriori / dryair)
        print("IMI averaging kernel:", avkern)

        print("IMI pressure weights reversed:", (dryair / np.sum(dryair))[::-1])
        print("IMI prior profile vmr reversed:", (apriori / dryair)[::-1])
        print("IMI averaging kernel reversed:", avkern[::-1])
        print("------------------------------------------\n")

        print(f"debug row: {debug_row}")
        print("iGC, jGC:", obs_mapped_to_gc["iGC"][debug_row], obs_mapped_to_gc["jGC"][debug_row])
        print("GC lat/lon:", obs_mapped_to_gc["lat"][debug_row], obs_mapped_to_gc["lon"][debug_row])
        print("sat lat/lon:", obs_mapped_to_gc["lat_sat"][debug_row], obs_mapped_to_gc["lon_sat"][debug_row])
        print("time:", obs_mapped_to_gc["time"][debug_row])
        print("satellite column:", obs_mapped_to_gc[species][debug_row])
        print("GOOPy virtual column ppb:", virtual_satellite[debug_row] * 1e9)
        print("p_sat:", obs_mapped_to_gc["p_sat"][debug_row])
        print("dry_air_subcolumns:", obs_mapped_to_gc["dry_air_subcolumns"][debug_row])
        print("apriori:", obs_mapped_to_gc["apriori"][debug_row])
        print("avkern:", obs_mapped_to_gc["avkern"][debug_row])
        print("observation_count:", obs_mapped_to_gc["observation_count"][debug_row])
        print("===================================\n")


    obs_GC[:, 1] = virtual_satellite * 1e9  # convert from mol/mol to ppb
    obs_GC[:, 2] = obs_mapped_to_gc["lon_sat"]
    obs_GC[:, 3] = obs_mapped_to_gc["lat_sat"]
    obs_GC[:, 4] = obs_mapped_to_gc["observation_count"]

    for strdate in all_strdate:
        gridcell_dict = obs_mapped_to_gc[obs_mapped_to_gc["time"] == strdate]
        sel_idx = np.where(obs_mapped_to_gc["time"] == strdate)[0]
        # virtual_satellite = get_virtual_satellite(
        #     strdate, gc_cache, gridcell_dict, n_elements, config
        # )
        if config["EnableOSSE"]:
            synthetic_virtual_satellite = get_virtual_satellite(
                strdate, osse_gc_cache, gridcell_dict, n_elements, config
            ) * 1e9  # convert to ppb

        # Save actual and virtual satellite data
        if config["EnableOSSE"]:
            # Synthetic observations if using OSSE, add random noise later
            obs_GC[sel_idx, 0] = synthetic_virtual_satellite
        else:
            # Actual satellite species column observation
            obs_GC[sel_idx, 0] = gridcell_dict[species]

    # add random noise to synthetic observations if using OSSE
    if config["EnableOSSE"]:
        noise = np.random.normal(
            loc=0.0,
            scale=float(config["ObsErrorOSSE"]),
            size=obs_GC[:,0].shape,
        )
        obs_GC[:,0] += noise

    # Output
    output = {}

    # Always return the coincident satellite and GEOS-Chem data
    output["obs_GC"] = obs_GC
    output["GC_index"] = GC_index

    return output


def format_goopy_satellite_output(
    filename,
    species,
    satellite_product,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    virtual_satellite,
    config,
    use_water_obs=False,
):
    """
    Format GOOPy output for the unaveraged satellite visualization path.

    GOOPy operates on the whole input file, while IMI visualization uses the
    same spatial/quality-filtered pixels as the native satellite operator. If
    GOOPy returned full-file output, subset it to those filtered pixels.
    """
    satellite, sat_ind = read_and_filter_satellite(
        filename, satellite_product, gc_startdate, gc_enddate,
        xlim, ylim, use_water_obs
    )

    n_obs = len(sat_ind[0])
    if n_obs == 0:
        return None

    satellite_shape = satellite["longitude"].shape
    sat_flat_ind = np.ravel_multi_index(sat_ind, satellite_shape)

    if len(virtual_satellite) == np.product(satellite_shape):
        virtual_satellite = virtual_satellite[sat_flat_ind]
    elif len(virtual_satellite) != n_obs:
        raise ValueError(
            f"GOOPy output has {len(virtual_satellite)} satellite columns, "
            f"but the filtered satellite data has {n_obs} observations and "
            f"the full satellite grid has {np.product(satellite_shape)} pixels."
        )

    obs_GC = np.empty([n_obs, 6], dtype=np.float32)
    obs_GC.fill(np.nan)

    i_sat = sat_ind[0]
    j_sat = sat_ind[1]

    obs_GC[:, 0] = satellite[species][i_sat, j_sat]
    obs_GC[:, 1] = virtual_satellite * 1e9  # convert from mol/mol to ppb
    obs_GC[:, 2] = satellite["longitude"][i_sat, j_sat]
    obs_GC[:, 3] = satellite["latitude"][i_sat, j_sat]
    obs_GC[:, 4] = i_sat
    obs_GC[:, 5] = j_sat

    output = {}
    output["obs_GC"] = obs_GC

    return output


def apply_average_satellite_operator(
    filename,
    species,
    satellite_product,
    satellite_cache,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    period_i,
    obs_mapped_to_gc,
    config,
    use_water_obs=False
):
    """
    Apply the averaging satellite operator to map GEOS-Chem data to satellite observation space.

    Arguments
        filename          [str]        : satellite netcdf data file to read
        species           [str]        : The species (CH4 or CO2) to use
        satellite_product [str]        : "BlendedTROPOMI", "TROPOMI", or "Other", specifying the data used in the inversion.
        n_elements        [int]        : Number of state vector elements
        gc_startdate      [datetime64] : First day of inversion period, for GEOS-Chem and satellite
        gc_enddate        [datetime64] : Last day of inversion period, for GEOS-Chem and satellite
        xlim              [float]      : Longitude bounds for simulation domain
        ylim              [float]      : Latitude bounds for simulation domain
        gc_cache          [str]        : Path to GEOS-Chem output data
        period_i       [int]        : kalman filter period
        obs_mapped_to_gc [numpy.ndarray] : structured array of grid-cell-averaged satellite observations mapped to GC gridcells
        config         [dict]       : dict of the config file
        use_water_obs  [bool]       : if True, use observations over water

    Returns
        output            [dict]       : Dictionary with:
                                        - obs_GC : GEOS-Chem and satellite data
                                        - satellite gas
                                        - GEOS-Chem gas
                                        - satellite lat, lon
                                        - satellite lat index, lon index
    """
    n_gridcells = len(obs_mapped_to_gc)
    gc_lat_lon = get_gc_lat_lon(gc_cache, gc_startdate)
    GC_shape = (len(gc_lat_lon['lat']), len(gc_lat_lon['lon']))

    # Initialize array with n_gridcells rows and 5 columns. Columns are
    # satellite species, GEOSChem species, longitude, latitude, observation counts
    obs_GC = np.empty([n_gridcells, 5], dtype=np.float32)
    obs_GC.fill(np.nan)

    if config['UseGCHP']:
        GC_index = np.ravel_multi_index((obs_mapped_to_gc["nfi"],
                                         obs_mapped_to_gc["Ydimi"],
                                         obs_mapped_to_gc["Xdimi"]), GC_shape)
    else:
        GC_index = np.ravel_multi_index((obs_mapped_to_gc["jGC"], #lat
                                         obs_mapped_to_gc["iGC"] # lon
                                         ), GC_shape)

    all_strdate = [gridcell["time"] for gridcell in obs_mapped_to_gc]
    all_strdate = list(set(all_strdate))

    # Read GEOS-Chem data for simulated truth in OSSE simulation
    if config["EnableOSSE"]:
        osse_gc_cache = "./data_geoschem_osse"

        # check if the osse_gc_cache exists
        assert os.path.exists(osse_gc_cache), (
            f"OSSE GEOS-Chem cache directory {osse_gc_cache} does not exist. "
            "Please run the OSSE simulation first."
        )

    for strdate in all_strdate:
        gridcell_dict = obs_mapped_to_gc[obs_mapped_to_gc["time"] == strdate]
        sel_idx = np.where(obs_mapped_to_gc["time"] == strdate)[0]
        virtual_satellite = get_virtual_satellite(
            strdate, gc_cache, gridcell_dict, n_elements, config
        )
        if config["EnableOSSE"]:
            synthetic_virtual_satellite = get_virtual_satellite(
                strdate, osse_gc_cache, gridcell_dict, n_elements, config
            ) * 1e9  # convert to ppb

        # Save actual and virtual satellite data
        if config["EnableOSSE"]:
            # Synthetic observations if using OSSE, add random noise later
            obs_GC[sel_idx, 0] = synthetic_virtual_satellite
        else:
            # Actual satellite species column observation
            obs_GC[sel_idx, 0] = gridcell_dict[species]

        obs_GC[sel_idx, 1] = virtual_satellite * 1e9  # Virtual satellite column observation and convert to ppb
        obs_GC[sel_idx, 2] = gridcell_dict["lon_sat"]  # satellite longitude
        obs_GC[sel_idx, 3] = gridcell_dict["lat_sat"]  # satellite latitude
        obs_GC[sel_idx, 4] = gridcell_dict["observation_count"]  # observation counts

    # add random noise to synthetic observations if using OSSE
    if config["EnableOSSE"]:
        noise = np.random.normal(
            loc=0.0,
            scale=float(config["ObsErrorOSSE"]),
            size=obs_GC[:,0].shape,
        )
        obs_GC[:,0] += noise

    # Output
    output = {}

    # Always return the coincident satellite and GEOS-Chem data
    output["obs_GC"] = obs_GC
    output["GC_index"] = GC_index

    return output


def apply_satellite_operator(
    filename,
    species,
    satellite_product,
    satellite_cache,
    n_elements,
    gc_startdate,
    gc_enddate,
    xlim,
    ylim,
    gc_cache,
    period_i,
    config,
    use_water_obs=False,
):
    """
    Apply the satellite operator to map GEOS-Chem species data to satellite observation space.

    Arguments
        filename           [str]        : Satellite netcdf data file to read
        species            [str]        : The species (CH4 or CO2) to use
        satellite_product  [str]        : "BlendedTROPOMI", "TROPOMI", or "Other", specifying the data used in the inversion.
        n_elements         [int]        : Number of state vector elements
        gc_startdate       [datetime64] : First day of inversion period, for GEOS-Chem and satellite
        gc_enddate         [datetime64] : Last day of inversion period, for GEOS-Chem and satellite
        xlim               [float]      : Longitude bounds for simulation domain
        ylim               [float]      : Latitude bounds for simulation domain
        gc_cache           [str]        : Path to GEOS-Chem output data
        period_i       [int]            : kalman filter period
        config         [dict]          : dict of the config file
        use_water_obs  [bool]          : if True, use observations over water

    Returns
        output             [dict]       : Dictionary with one or two fields:
                                                    - obs_GC : GEOS-Chem and satellite species data
                                                    - satellite species
                                                    - GEOS-Chem species
                                                    - satellite lat, lon
                                                    - satellite lat index, lon index
                                                      If build_jacobian=True, also include:
                                                        - K      : Jacobian matrix
    """

    # Read satellite data
    satellite, sat_ind = read_and_filter_satellite(
        filename, satellite_product, gc_startdate, gc_enddate,
        xlim, ylim, use_water_obs)

    # Number of satellite observations
    n_obs = len(sat_ind[0])
    if n_obs == 0:
        return None
    # print("Found", n_obs, "satellite observations.")

    # Initialize a list to store the dates we want to look at
    all_strdate = []

    # Define time threshold (hour 00 after the inversion period)
    date_after_inversion = str(gc_enddate + np.timedelta64(1, "D"))[:10].replace(
        "-", ""
    )
    time_threshold = f"{date_after_inversion}_00"

    # For each satellite observation
    for k in range(n_obs):
        # Get the date and hour
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index
        time = pd.to_datetime(str(satellite["time"][iSat,jSat]))
        strdate = get_strdate(time, time_threshold)
        all_strdate.append(strdate)
    all_strdate = list(set(all_strdate))

    # Read GEOS_Chem data for the dates of interest
    all_date_gc = read_all_geoschem(
        all_strdate, gc_cache, config
    )

    # Initialize array with n_obs rows and 6 columns. Columns are satellite
    # mixing ratio, GEOSChem mixing ratio, longitude, latitude, II, JJ
    obs_GC = np.zeros([n_obs, 6], dtype=np.float32)
    obs_GC.fill(np.nan)

    if config['UseGCHP']:
        if config['STRETCH_GRID']:
            sf_formatted = f"{config['STRETCH_FACTOR']:.2f}".replace(".", "d")
            target_geohash = pgh.encode(config['TARGET_LAT'], config['TARGET_LON'])
            gridspec_path = f"c{config['CS_RES']}_s{sf_formatted}_t{target_geohash}_gridspec.nc"
        else:
            gridspec_path = f"c{config['CS_RES']}_gridspec.nc"
        GC_shape = (6, config['CS_RES'], config['CS_RES'])
        CSgridDir = f"{os.path.expandvars(config['OutputPath']) }/{config['RunName']}/CS_grids"

        overlap_area_all = get_overlap_area_CSgrid(satellite, filename, sat_ind, CSgridDir,
                            gridspec_path, GC_shape) # (n_dst, n_valid_obs)

    # For each satellite observation:
    for k in range(n_obs):

        # Get GEOS-Chem data for the date of the observation:
        iSat = sat_ind[0][k]
        jSat = sat_ind[1][k]
        longitude_bounds = satellite["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = satellite["latitude_bounds"][iSat, jSat, :]

        p_sat = satellite["pressures"][iSat, jSat, :]
        dry_air_subcolumns = satellite["dry_air_subcolumns"][iSat, jSat, :]  # mol m-2
        apriori = satellite["profile_apriori"][iSat, jSat, :]  # mol m-2
        avkern = satellite["column_AK"][iSat, jSat, :]
        time = pd.to_datetime(str(satellite["time"][iSat,jSat]))
        strdate = get_strdate(time, time_threshold)
        GEOSCHEM = all_date_gc[strdate]

        if config['UseGCHP']:
            overlap_area_csr = overlap_area_all[:, k]
            gc_coords = overlap_area_csr.nonzero()[0]       # row indices of non-zero entries
            overlap_area = overlap_area_csr.data
        else:
            # Polygon representing satellite pixel
            polygon_satellite = Polygon(np.column_stack((longitude_bounds, latitude_bounds)))

            dlon = np.median(np.diff(GEOSCHEM["lon"]))  # GEOS-Chem lon resolution
            dlat = np.median(np.diff(GEOSCHEM["lat"]))  # GEOS-Chem lon resolution

            # Find GEOS-Chem lats & lons closest to the corners of the satellite pixel
            corners_lon_index = []
            corners_lat_index = []
            for l in range(4):
                iGC = nearest_loc(
                    longitude_bounds[l], GEOSCHEM["lon"], tolerance=max(dlon, 0.5)
                )
                jGC = nearest_loc(
                    latitude_bounds[l], GEOSCHEM["lat"], tolerance=max(dlat, 0.5)
                )
                corners_lon_index.append(iGC)
                corners_lat_index.append(jGC)
            # If the tolerance in nearest_loc() is not satisfied, skip the observation
            if np.nan in corners_lon_index + corners_lat_index:
                continue
            # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the satellite corners
            ij_GC = [(x, y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
            gc_coords = [(GEOSCHEM["lon"][i], GEOSCHEM["lat"][j]) for i, j in ij_GC]

            # Compute the overlapping area between the satellite pixel and GEOS-Chem grid cells it touches
            overlap_area = np.zeros(len(gc_coords), dtype=np.float32)
            # For each GEOS-Chem grid cell that touches the satellite pixel:
            for gridcellIndex in range(len(gc_coords)):
                # Define polygon representing the GEOS-Chem grid cell
                coords = gc_coords[gridcellIndex]
                geoschem_corners_lon = [
                    coords[0] - dlon / 2,
                    coords[0] + dlon / 2,
                    coords[0] + dlon / 2,
                    coords[0] - dlon / 2,
                ]
                geoschem_corners_lat = [
                    coords[1] - dlat / 2,
                    coords[1] - dlat / 2,
                    coords[1] + dlat / 2,
                    coords[1] + dlat / 2,
                ]
                # If this is a global 2.0 x 2.5 grid, extend the eastern-most grid cells to 180 degrees
                if (dlon == 2.5) and (coords[0] == 177.5):
                    for i in [1, 2]:
                        geoschem_corners_lon[i] += dlon / 2

                polygon_geoschem = Polygon(
                    np.column_stack((geoschem_corners_lon, geoschem_corners_lat))
                )

                # Calculate overlapping area as the intersection of the two polygons
                if polygon_geoschem.intersects(polygon_satellite):
                    overlap_area[gridcellIndex] = polygon_satellite.intersection(
                        polygon_geoschem
                    ).area

        # If there is no overlap between GEOS-Chem and satellite, skip to next observation:
        if np.sum(overlap_area) == 0:
            continue

        # =======================================================
        #       Map GEOS-Chem to satellite observation space
        # =======================================================

        # Otherwise, initialize satellite virtual mixing ratios and virtual
        #  sensitivity as zero
        area_weighted_virtual_satellite = 0  # virtual satellite mixing ratio

        # For each GEOS-Chem grid cell that touches the satellite pixel:
        for gridcellIndex in range(len(gc_coords)):
            if overlap_area[gridcellIndex] == 0:
                continue
            if config['UseGCHP']:
                # Get GEOS-Chem lat/lon indices for the cell
                f, j, x = np.unravel_index(gc_coords[gridcellIndex], GC_shape)

                # Get GEOS-Chem pressure edges for the cell
                p_gc = GEOSCHEM["PEDGE"][f, j, x, :]

                # Get GEOS-Chem mixing ratios for the cell
                gc_species = GEOSCHEM[species][f, j, x, :]
            else:
                # Get GEOS-Chem lat/lon indices for the cell
                iGC, jGC = ij_GC[gridcellIndex]

                # Get GEOS-Chem pressure edges for the cell
                p_gc = GEOSCHEM["PEDGE"][iGC, jGC, :]

                # Get GEOS-Chem mixing ratios for the cell
                gc_species = GEOSCHEM[species][iGC, jGC, :]

            # Get merged GEOS-Chem/satellite pressure grid for the cell
            merged = merge_pressure_grids(p_sat, p_gc)

            # Remap GEOS-Chem mixing ratios to satellite pressure levels
            sat_species = remap(
                gc_species,
                merged["data_type"],
                merged["p_merge"],
                merged["edge_index"],
                merged["first_gc_edge"],
            )  # ppb

            # Convert volumetric mixing ratio to mol m-2
            sat_species_molm2 = sat_species * 1/mixing_ratio_conv_factor(species) * dry_air_subcolumns  # mol m-2

            # Derive the column-averaged mixing ratio that satellite would
            # see over this ground cell
            virtual_satellite_gridcellIndex = apply_averaging_kernel(
                apriori, avkern, sat_species_molm2, dry_air_subcolumns, species
            ) # Volumetric mixing ratio

            # Weight by overlapping area (to be divided out later) and add to sum
            area_weighted_virtual_satellite += (
                overlap_area[gridcellIndex] * virtual_satellite_gridcellIndex
            )  # ppb m2

        # Compute virtual satellite observation as weighted mean by overlapping area
        # i.e., need to divide out area [m2] from the previous step
        virtual_satellite = area_weighted_virtual_satellite / sum(overlap_area)

        # Save actual and virtual satellite data
        obs_GC[k, 0] = satellite[species][
            iSat, jSat
        ]  # Actual satellite mixing ratio column observation
        obs_GC[k, 1] = virtual_satellite  # Virtual satellite mixing ratio column observation
        obs_GC[k, 2] = satellite["longitude"][iSat, jSat]  # satellite longitude
        obs_GC[k, 3] = satellite["latitude"][iSat, jSat]  # satellite latitude
        obs_GC[k, 4] = iSat  # satellite index of longitude
        obs_GC[k, 5] = jSat  # satellite index of latitude

    # Output
    output = {}

    # Always return the coincident satellite and GEOS-Chem data
    output["obs_GC"] = obs_GC

    return output


def average_satellite_observations(
        satellite, species, gc_lat_lon, sat_ind, time_threshold
    ):
    """
    Map satellite observations into appropriate gc gridcells. Then average all
    observations within a gridcell for processing. Use area weighting if
    observation overlaps multiple gridcells.

    Arguments
        satellite      [dict]   : Dict of satellite data
        species        [str]    : Name of species analyzed (CO2 or CH4)
        gc_lat_lon     [list]   : list of dictionaries containing  gc gridcell info
        sat_ind        [int]    : index list of satellite data that passes filters

    Returns
        numpy.ndarray
            Structured array of grid-cell-averaged values with fields:

            - iGC, jGC : int
                GEOS-Chem longitude and latitude indices
            - lat, lon : float
                Grid cell center coordinates
            - lat_sat, lon_sat : float
                Weighted-average satellite footprint center
            - species : float
                Weighted-average species column
            - time : str
                Averaged observation time (string from `get_strdate`)
            - p_sat : float[n_lev]
                Weighted-average vertical pressure profile
            - dry_air_subcolumns : float[n_lev]
                Weighted-average dry-air subcolumn profile
            - apriori : float[n_lev]
                Weighted-average a priori species profile
            - avkern : float[n_lev]
                Weighted-average averaging kernel
            - observation_count : float
                Effective number of contributing observations (fractional if
                split across multiple grid cells)

    """
    n_obs = len(sat_ind[0])
    # print("Found", n_obs, "satellite observations.")
    gc_lats = gc_lat_lon["lat"]
    gc_lons = gc_lat_lon["lon"]
    dlon = np.median(np.diff(gc_lat_lon["lon"]))  # GEOS-Chem lon resolution
    dlat = np.median(np.diff(gc_lat_lon["lat"]))  # GEOS-Chem lon resolution
    gridcell_dicts = get_gridcell_list(gc_lons, gc_lats, species)

    for k in range(n_obs):
        iSat = sat_ind[0][k]  # lat index
        jSat = sat_ind[1][k]  # lon index

        # Find GEOS-Chem lats & lons closest to the corners of the satellite pixel
        longitude_bounds = satellite["longitude_bounds"][iSat, jSat, :]
        latitude_bounds = satellite["latitude_bounds"][iSat, jSat, :]
        corners_lon_index = []
        corners_lat_index = []

        for l in range(4):
            iGC = nearest_loc(longitude_bounds[l], gc_lons, tolerance=max(dlon, 0.5))
            jGC = nearest_loc(latitude_bounds[l], gc_lats, tolerance=max(dlat, 0.5))
            corners_lon_index.append(iGC)
            corners_lat_index.append(jGC)

        # If the tolerance in nearest_loc() is not satisfied, skip the observation
        if np.nan in corners_lon_index + corners_lat_index:
            continue

        # Get lat/lon indexes and coordinates of GEOS-Chem grid cells closest to the satellite corners
        ij_GC = [(x, y) for x in set(corners_lon_index) for y in set(corners_lat_index)]
        gc_coords = [(gc_lons[i], gc_lats[j]) for i, j in ij_GC]

        # Compute the overlapping area between the satellite pixel and GEOS-Chem grid cells it touches
        overlap_area = np.zeros(len(gc_coords))

        # Polygon representing satellite pixel
        polygon_satellite = Polygon(np.column_stack((longitude_bounds, latitude_bounds)))
        for gridcellIndex in range(len(gc_coords)):
            # Define polygon representing the GEOS-Chem grid cell
            coords = gc_coords[gridcellIndex]
            geoschem_corners_lon = [
                coords[0] - dlon / 2,
                coords[0] + dlon / 2,
                coords[0] + dlon / 2,
                coords[0] - dlon / 2,
            ]
            geoschem_corners_lat = [
                coords[1] - dlat / 2,
                coords[1] - dlat / 2,
                coords[1] + dlat / 2,
                coords[1] + dlat / 2,
            ]
            polygon_geoschem = Polygon(
                np.column_stack((geoschem_corners_lon, geoschem_corners_lat))
            )

            # Calculate overlapping area as the intersection of the two polygons
            if polygon_geoschem.intersects(polygon_satellite):
                overlap_area[gridcellIndex] = polygon_satellite.intersection(
                    polygon_geoschem
                ).area
        # If there is no overlap between GEOS-Chem and satellite, skip to next observation:
        total_overlap_area = sum(overlap_area)
        if total_overlap_area == 0:
            continue

        # iterate through any gridcells with observation overlap
        # weight each observation if observation extent overlaps with multiple
        # gridcells
        for index, overlap in enumerate(overlap_area):
            if not overlap == 0:
                # get the matching dictionary for the gridcell with the overlap
                gridcell_dict = gridcell_dicts[ij_GC[index][0]][ij_GC[index][1]]
                gridcell_dict["lat_sat"].append(satellite["latitude"][iSat, jSat])
                gridcell_dict["lon_sat"].append(satellite["longitude"][iSat, jSat])
                gridcell_dict["overlap_area"].append(overlap)
                gridcell_dict["p_sat"].append(satellite["pressures"][iSat, jSat, :])
                gridcell_dict["surface_pressure"].append(
                    satellite["surface_pressure"][iSat, jSat]
                )
                gridcell_dict["nir_albedo"].append(satellite["nir_albedo"][iSat, jSat])
                gridcell_dict["swir_albedo"].append(satellite["swir_albedo"][iSat, jSat])
                gridcell_dict["dry_air_subcolumns"].append(
                    satellite["dry_air_subcolumns"][iSat, jSat, :]
                )
                gridcell_dict["apriori"].append(
                    satellite["profile_apriori"][iSat, jSat, :]
                )
                gridcell_dict["avkern"].append(satellite["column_AK"][iSat, jSat, :])
                gridcell_dict[
                    "time"
                ].append(  # convert times to epoch time to make taking the mean easier
                    int(pd.to_datetime(str(satellite["time"][iSat,jSat])).strftime("%s"))
                )
                gridcell_dict[species].append(
                    satellite[species][iSat, jSat]
                )  # Actual satellite mixing ratio column observation
                # record weights for averaging later
                gridcell_dict["observation_weights"].append(
                    overlap / total_overlap_area
                )
                # increment the observation count based on overlap area
                gridcell_dict["observation_count"] += overlap / total_overlap_area

    # filter out gridcells without any observations
    gridcell_dicts = [
        item for item in gridcell_dicts.flatten() if item["observation_count"] > 0
    ]
    # weighted average observation values for each gridcell
    for gridcell_dict in gridcell_dicts:
        gridcell_dict["lat_sat"] = np.average(
            gridcell_dict["lat_sat"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["lon_sat"] = np.average(
            gridcell_dict["lon_sat"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["overlap_area"] = np.average(
            gridcell_dict["overlap_area"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict[species] = np.average(
            gridcell_dict[species],
            weights=gridcell_dict["observation_weights"],
        )
        # take mean of epoch times and then convert gc filename time string
        time = pd.to_datetime(
            datetime.datetime.fromtimestamp(int(np.mean(gridcell_dict["time"])))
        )
        gridcell_dict["time"] = get_strdate(time, time_threshold)
        # for multi-dimensional arrays, we only take the average across the 0 axis
        gridcell_dict["p_sat"] = np.average(
            gridcell_dict["p_sat"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["surface_pressure"] = np.average(
            gridcell_dict["surface_pressure"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["nir_albedo"] = np.average(
            gridcell_dict["nir_albedo"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["swir_albedo"] = np.average(
            gridcell_dict["swir_albedo"],
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["dry_air_subcolumns"] = np.average(
            gridcell_dict["dry_air_subcolumns"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["apriori"] = np.average(
            gridcell_dict["apriori"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )
        gridcell_dict["avkern"] = np.average(
            gridcell_dict["avkern"],
            axis=0,
            weights=gridcell_dict["observation_weights"],
        )

    if not gridcell_dicts:
        # nothing to return
        return np.zeros(0, dtype=[
            ("iGC","i4"), ("jGC","i4"),
            ("lat_sat","f4"), ("lon_sat","f4"),
            (species,"f4"), ("time","U13"),
            ("p_sat","f4",(0,)),
            ("surface_pressure","f4"),
            ("nir_albedo","f4"),
            ("swir_albedo","f4"),
            ("dry_air_subcolumns","f4",(0,)),
            ("apriori","f4",(0,)),
            ("avkern","f4",(0,)),
            ("layer","f4",(0,)),
            ("observation_count","f4"),
            ("lat","f4"), ("lon","f4"),
        ])

    # infer vertical sizes from the first item
    n_lev_p       = len(gridcell_dicts[0]["p_sat"])
    n_lev_dryair  = len(gridcell_dicts[0]["dry_air_subcolumns"])
    n_lev_apriori = len(gridcell_dicts[0]["apriori"])
    n_lev_avkern  = len(gridcell_dicts[0]["avkern"])
    n_layers  = len(gridcell_dicts[0]["layer"])

    dtype_latlon = [
        ("iGC","i4"), ("jGC","i4"),
        ("lat_sat","f4"), ("lon_sat","f4"),
        (species,"f4"), ("time","U13"),
        ("p_sat","f4",(n_lev_p,)),
        ("surface_pressure","f4"),
        ("nir_albedo","f4"),
        ("swir_albedo","f4"),
        ("dry_air_subcolumns","f4",(n_lev_dryair,)),
        ("apriori","f4",(n_lev_apriori,)),
        ("avkern","f4",(n_lev_avkern,)),
        ("layer","f4",(n_layers,)),
        ("observation_count","f4"),
        ("lat","f4"), ("lon","f4"),
    ]

    arr = np.zeros(len(gridcell_dicts), dtype=dtype_latlon)

    for idx, cell in enumerate(gridcell_dicts):
        arr["iGC"][idx]   = cell["iGC"]
        arr["jGC"][idx]   = cell["jGC"]
        arr["lat_sat"][idx] = np.float32(cell["lat_sat"])
        arr["lon_sat"][idx] = np.float32(cell["lon_sat"])
        arr[species][idx] = np.float32(cell[species])
        arr["time"][idx]    = cell["time"]  # already a short string from get_strdate
        arr["p_sat"][idx]   = np.asarray(cell["p_sat"], dtype=np.float32)
        arr["surface_pressure"][idx] = np.float32(cell["surface_pressure"])
        arr["nir_albedo"][idx] = np.float32(cell["nir_albedo"])
        arr["swir_albedo"][idx] = np.float32(cell["swir_albedo"])
        arr["dry_air_subcolumns"][idx] = np.asarray(cell["dry_air_subcolumns"], dtype=np.float32)
        arr["apriori"][idx] = np.asarray(cell["apriori"], dtype=np.float32)
        arr["avkern"][idx]  = np.asarray(cell["avkern"], dtype=np.float32)
        arr["observation_count"][idx] = np.float32(cell["observation_count"])
        arr["lat"][idx] = np.float32(cell["lat"])
        arr["lon"][idx] = np.float32(cell["lon"])
        # arr["layer"][idx] = np.asarray(cell["layer"], dtype=np.float32)
        arr["layer"][idx] = np.arange(n_layers)

    return arr

# TODO: update this to match what average_sat_observations does with layers
def average_satellite_observations_to_CSgrid(
        satellite, species, filename, sat_ind, time_threshold,
        CSgridDir, gridspec_path, GC_shape):
    """
    Map satellite observations into appropriate GCHP gridcells. Then average all
    observations within a gridcell for processing.
    Use area weighting (conservative regridding weights x destination area)
    if observation overlaps multiple gridcells.

    Parameters
    ----------
    satellite : dict
        Dictionary of satellite satellite data.
    filename : str
        Path to the satellite dataset file.
    sat_ind : int
        Indices of valid satellite observations (after filtering).
    time_threshold : str or float
        Threshold used to bin/average times into representative strings.
    CSgridDir : str
        Path to cubed-sphere grid definitions.
    gridspec_path : str
        Path to the GCHP gridspec file.
    GC_shape : tuple
        Shape of the GEOS-Chem cubed-sphere grid.

    Returns
    -------
    output_dicts : np.ndarray of dict-like records
        Array of per-gridcell averaged observation data with fields:
            - nfi, Ydimi, Xdimi : cubed-sphere indices
            - lat_sat, lon_sat  : averaged satellite coordinates
            - species           : averaged species column
            - time              : averaged time string
            - p_sat             : averaged pressure grid
            - dry_air_subcolumns: averaged dry air columns
            - apriori           : averaged prior species profiles
            - avkern            : averaged averaging kernels
            - observation_count : number of obs contributing to the average
    """
    Sat_shape = satellite["longitude"].shape
    # flatten indices to valid obs
    sat_ind_flat = np.ravel_multi_index(sat_ind, Sat_shape)
    sat_mask = np.zeros(np.product(Sat_shape), dtype=bool)
    sat_mask[sat_ind_flat] = True

    overlap_area = get_overlap_area_CSgrid(satellite, filename, sat_ind, CSgridDir,
                            gridspec_path, GC_shape) # (n_dst, n_valid_obs)

    # Observation-centric weights
    # observation weights w_ij = overlap_ij / sum_j(overlap_ij)
    # where w_ij is the observation weight for destination grid cell j from the observation cell i
    # overlap_ij is the overlap area between destination cell j and the observation cell i
    # total overlap area for each observation is the sum of the overlapping area
    # between observation cell i and all its overlapping destination grid cell j
    total_overlap_area = np.asarray(overlap_area.sum(axis=0)).ravel() # (n_valid_obs)
    inv_total_overlap_area = np.zeros_like(total_overlap_area, dtype=np.float32)
    valid_obs = total_overlap_area > 0
    inv_total_overlap_area[valid_obs] = 1.0 / total_overlap_area[valid_obs]

    # obs_weight: each obs contributes fully across overlapping cells
    # element-wise multiplication:
    obs_weights = overlap_area.multiply(inv_total_overlap_area).tocsr() # (n_dst × n_valid_obs)

    sum_obs_weights = np.asarray(obs_weights.sum(axis=1)).ravel() # (n_dst)
    GC_mask = sum_obs_weights > 0
    GC_indices = np.where(GC_mask)[0]  # only valid destination cells
    n_valid_GC = len(GC_indices)

    # ---------- Scalars ----------
    sat_lon_flat = satellite["longitude"].ravel()[sat_mask]
    sat_lat_flat = satellite["latitude"].ravel()[sat_mask]
    sat_species_flat = satellite[species].ravel()[sat_mask]
    surface_pressure_flat = satellite["surface_pressure"].ravel()[sat_mask]
    nir_albedo_flat = satellite["nir_albedo"].ravel()[sat_mask]
    swir_albedo_flat = satellite["swir_albedo"].ravel()[sat_mask]
    sat_time_flat = pd.to_datetime(satellite["time"].ravel()[sat_mask])

    # For each destination grid cell j which overlap with at least one observation cell i,
    # weighted_obs_j = sum_i(obs_weights_ij * obs_i) / sum_i(obs_weights_ij)
    sat_lon_avg = (obs_weights[GC_indices, :] @ sat_lon_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    sat_lat_avg = (obs_weights[GC_indices, :] @ sat_lat_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    sat_species_avg = (obs_weights[GC_indices, :] @ sat_species_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    surface_pressure_avg = (obs_weights[GC_indices, :] @ surface_pressure_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    nir_albedo_avg = (obs_weights[GC_indices, :] @ nir_albedo_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    swir_albedo_avg = (obs_weights[GC_indices, :] @ swir_albedo_flat).astype(np.float32) / sum_obs_weights[GC_indices]
    sat_time_avg = (obs_weights[GC_indices, :] @ sat_time_flat.astype(np.int64)) / sum_obs_weights[GC_indices]

    sat_time_avg = pd.to_datetime(sat_time_avg)  # convert back to datetime
    sat_time_str = np.array([get_strdate(t, time_threshold) for t in sat_time_avg])

    # --- Multi-dimensional arrays ---
    def weighted_avg_profiles(var_flat):
        # var_flat: (n_valid_obs, n_lev)
        # obs_weights: (n_dst, n_valid_obs)
        # sum_obs_weights: (n_dst,)
        # result: (n_valid_GC,)
        result = obs_weights[GC_indices, :] @ var_flat
        result /= sum_obs_weights[GC_indices, None]
        return result.astype(np.float32)

    p_sat_flat = species["pressures"].reshape(-1, species["pressures"].shape[-1])[sat_mask, :]
    dryair_flat = species["dry_air_subcolumns"].reshape(-1, species["dry_air_subcolumns"].shape[-1])[sat_mask, :]
    apriori_flat = species["profile_apriori"].reshape(-1, species["profile_apriori"].shape[-1])[sat_mask, :]
    avkern_flat = species["column_AK"].reshape(-1, species["column_AK"].shape[-1])[sat_mask, :]

    p_sat_avg = weighted_avg_profiles(p_sat_flat)
    dryair_avg = weighted_avg_profiles(dryair_flat)
    apriori_avg = weighted_avg_profiles(apriori_flat)
    avkern_avg = weighted_avg_profiles(avkern_flat)

    # --- Fill dictionaries ---
    f_idx, j_idx, x_idx = np.unravel_index(GC_indices, GC_shape)
    n_lev_p = p_sat_avg.shape[1]
    n_lev_dryair = dryair_avg.shape[1]
    n_lev_apriori = apriori_avg.shape[1]
    n_lev_avkern = avkern_avg.shape[1]

    dtype = [
        ("nfi", "i4"), ("Ydimi", "i4"), ("Xdimi", "i4"),
        ("lat_sat", "f4"), ("lon_sat", "f4"), (species, "f4"),
        ("time", "U13"), ("p_sat", "f4", (n_lev_p,)),
        ("surface_pressure", "f4"),
        ("nir_albedo", "f4"),
        ("swir_albedo", "f4"),
        ("dry_air_subcolumns", "f4", (n_lev_dryair,)),
        ("apriori", "f4", (n_lev_apriori,)), ("avkern", "f4", (n_lev_avkern,)),
        ("observation_count", "f4")
    ]

    output_dicts = np.zeros(n_valid_GC, dtype=dtype)
    output_dicts["nfi"] = f_idx
    output_dicts["Ydimi"] = j_idx
    output_dicts["Xdimi"] = x_idx
    output_dicts["lat_sat"] = sat_lat_avg
    output_dicts["lon_sat"] = sat_lon_avg
    output_dicts[species] = sat_species_avg
    output_dicts["time"] = sat_time_str
    output_dicts["p_sat"] = p_sat_avg
    output_dicts["surface_pressure"] = surface_pressure_avg
    output_dicts["nir_albedo"] = nir_albedo_avg
    output_dicts["swir_albedo"] = swir_albedo_avg
    output_dicts["dry_air_subcolumns"] = dryair_avg
    output_dicts["apriori"] = apriori_avg
    output_dicts["avkern"] = avkern_avg
    output_dicts["observation_count"] = sum_obs_weights[GC_indices]

    return output_dicts


def virtual_satellite_species_and_pedge(date, gc_cache, gridcell_dict, n_elements, config):
    """
    Read species and pressure edge data from GEOS-Chem output for the grid cells of interest. 
    """
    UseGCHP = config['UseGCHP']

    # Assemble file paths to GEOS-Chem output collections for input data
    file_species = f"GEOSChem.SpeciesConc.{date}00z.nc4"
    file_pedge = f"GEOSChem.StateMetLevEdge.{date}00z.nc4"

    # Read lat, lon, species from the SpeciecConc collection
    filename = f"{gc_cache}/{file_species}"

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(filename, chunks='auto') as gc_data_all:
            if gc_data_all.sizes.get("time", 0) == 0:
                print(f"ERROR: {filename}: empty time dimension", flush=True)
            if UseGCHP:
                nfi   = xr.DataArray(gridcell_dict["nfi"],   dims="obs")
                Ydimi = xr.DataArray(gridcell_dict["Ydimi"], dims="obs")
                Xdimi = xr.DataArray(gridcell_dict["Xdimi"], dims="obs")

                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    nf=nfi,
                    Ydim=Ydimi,
                    Xdim=Xdimi,
                    drop=True
                )
            else:
                jGC   = xr.DataArray(gridcell_dict["jGC"],   dims="obs")
                iGC   = xr.DataArray(gridcell_dict["iGC"],   dims="obs")
                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    lat=jGC,
                    lon=iGC,
                    drop=True
                )
            species = gc_data[f"SpeciesConcVV_{config['Species']}"].transpose("obs","lev").values

    # Read PEDGE from the StateMetLevEdge collection
    filename = f"{gc_cache}/{file_pedge}"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(filename, chunks='auto') as gc_data_all:
            if gc_data_all.sizes.get("time", 0) == 0:
                print(f"ERROR: {filename}: empty time dimension", flush=True)
            if UseGCHP:
                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    nf=nfi,
                    Ydim=Ydimi,
                    Xdim=Xdimi,
                    drop=True
                )
            else:
                gc_data = gc_data_all.isel(
                    time=0).squeeze().isel(
                    lat=jGC,
                    lon=iGC,
                    drop=True
                )
            lev_dim = "lev" if "lev" in gc_data["Met_PEDGE"].dims else "ilev"
            PEDGE = gc_data["Met_PEDGE"].transpose("obs", lev_dim).values

    return species, PEDGE

def get_virtual_satellite(
    date, gc_cache, gridcell_dict, n_elements, config, 
):
    """
    Generate virtual satellite species observations from GEOS-Chem.

    Extracts species and pressure from GEOS-Chem, remaps to satellite layers,
    and applies averaging kernels. 

    Parameters
    ----------
    date : str
        Date of interest ("YYYYMMDD_HH").
    gc_cache : str
        Path to GEOS-Chem output files.
    gridcell_dict : dict
        Gridcell info with obs indices and satellite data.
    n_elements : int
        Number of state vector elements.
    config : dict
        Inversion configuration options.

    Returns
    -------
    If build_jacobian=False:
        ndarray (N,) of virtual satellite columns.
    If build_jacobian=True:
        (perturbation columns, base columns, final columns).
    """

    species, PEDGE = virtual_satellite_species_and_pedge(date, gc_cache, gridcell_dict, n_elements, config)

    n_superobs = len(gridcell_dict)
    virtual_satellite = np.empty([n_superobs, ], dtype=np.float32)
    virtual_satellite.fill(np.nan)

    p_sat = gridcell_dict["p_sat"]
    dry_air_subcolumns = gridcell_dict["dry_air_subcolumns"]  # mol m-2
    apriori = gridcell_dict["apriori"]  # mol m-2
    avkern = gridcell_dict["avkern"]

    # (N, S, G)  sums to 1 along G, where
    # N is the number of super observations
    # S is the number of satellite pressure edges
    # G is the number of GEOS-Chem pressure edges
    vertical_weights = remapping_weights(p_sat, PEDGE)
    sat_species = np.einsum("nsg,ng->ns", vertical_weights, species)         # (N, S)
    sat_species_molm2 = sat_species * dry_air_subcolumns                     # (N, S)
    virtual_satellite = (
        np.sum(apriori + avkern * (sat_species_molm2 - apriori), axis=1) / \
        np.sum(dry_air_subcolumns, axis=1)
    ).astype(np.float32)                      # (N,), unitless mixing ratio

    return virtual_satellite


def get_virtual_satellite_pert_and_base(
    date, gc_cache, gridcell_dict, n_elements, config, 
):
    """
    Compute virtual satellite species columns from GEOS-Chem perturbation and base simulations.

    Parameters
    ----------
    date : str
        Date of interest ("YYYYMMDD_HH").
    gc_cache : str
        Path to GEOS-Chem output files.
    gridcell_dict : dict
        Gridcell info with obs indices and satellite data.
    n_elements : int
        Number of state vector elements.
    config : dict
        Inversion configuration options.

    Returns
    -------
    A tuple of (perturbation columns, base columns)
    """
    # Read sensitivity data from GEOS-Chem perturbation simulations
    emis_elements = n_elements
    if config['OptimizeOH']:
        emis_elements -= 2 if config['isRegional'] else 1
    if config['OptimizeBCs']:
        emis_elements -= 4
    ntracers = config["NumJacobianTracers"]
    opt_OH = config["OptimizeOH"]
    opt_BC = config["OptimizeBCs"]
    is_Regional = config["isRegional"]

    num_BC = 4
    if is_Regional:
        num_OH = 1
    else:
        num_OH = 2

    n_base_runs = (
        n_elements - int(opt_OH * num_OH) - (int(opt_BC) * num_BC)
    ) / ntracers

    nruns = (
        np.ceil(n_base_runs).astype(int)
        + (int(opt_OH) * num_OH)
        + (int(opt_BC) * num_BC)
    )

    # Dictionary that stores mapping of state vector elements to
    # perturbation simulation numbers
    pert_simulations_dict = {}
    for e in range(n_elements):
        # State vector elements are numbered 1..nelements
        sv_elem = e + 1

        is_OH_element = check_is_OH_element(
            sv_elem, n_elements, opt_OH, is_Regional
        )
        is_BC_element = check_is_BC_element(
            sv_elem, n_elements, opt_OH, opt_BC, is_OH_element, is_Regional
        )
        # Determine which run directory to look in
        if is_OH_element:
            if is_Regional:
                run_number = nruns
            else:
                num_back = n_elements % sv_elem
                run_number = nruns - num_back
        elif is_BC_element:
            num_back = n_elements % sv_elem
            run_number = nruns - num_back
        else:
            run_number = np.ceil(sv_elem / ntracers).astype(int)

        run_num = str(run_number).zfill(4)

        # add the element to the dictionary for the relevant simulation number
        if run_num not in pert_simulations_dict:
            pert_simulations_dict[run_num] = [sv_elem]
        else:
            pert_simulations_dict[run_num].append(sv_elem)

    _, PEDGE = virtual_satellite_species_and_pedge(date, gc_cache, gridcell_dict, n_elements, config)
    p_sat = gridcell_dict["p_sat"]
    vertical_weights = remapping_weights(p_sat, PEDGE)

    gc_date = pd.to_datetime(date, format="%Y%m%d_%H")
    virtual_satellite_pert = [
        get_virtual_satellite_pert(gc_date, k, gridcell_dict, config, v, n_elements, vertical_weights)
        for k, v in pert_simulations_dict.items()
    ]

    virtual_satellite_pert = np.concatenate(virtual_satellite_pert, axis=1)

    virtual_satellite_base = get_virtual_satellite_pert(
        gc_date, "0001", gridcell_dict, config, [0], n_elements, vertical_weights, baserun=True
    )

    return virtual_satellite_pert, virtual_satellite_base


def get_virtual_satellite_pert(gc_date, run_id, gridcell_dict, config, sv_elems, n_elements, vertical_weights, baserun=False):
    """
    Compute virtual satellite species columns from a single GEOS-Chem Jacobian run.

    Extracts species tracer(s) from the specified perturbation (or base) simulation,
    remaps them to satellite pressure layers using vertical weights, and applies
    averaging kernels to generate virtual satellite observations.

    Parameters
    ----------
    gc_date : pd.Datetime
        Date and time of the simulation output.
    run_id : str
        ID of the Jacobian GEOS-Chem run (e.g., "0001").
    gridcell_dict : dict
        Gridcell info containing observation indices and satellite data.
    config : dict
        Inversion configuration options.
    sv_elems : list
        State vector elements included in this simulation.
    n_elements : int
        Total number of state vector elements.
    vertical_weights : np.ndarray
        Weights to remap GEOS-Chem levels to satellite layers.
    baserun : bool, optional
        If True, only process the base species (default False).

    Returns
    -------
    virtual_satellite_all : np.ndarray, shape (N, n_tracers)
        Virtual satellite columns for all super-observations (N) and
        tracer elements in this run, in unitless mixing ratio.
    """
    prefix = os.path.expandvars(
        config["OutputPath"] + "/" + config["RunName"] + "/jacobian_runs"
    )
    j_dir = f"{prefix}/{config['RunName']}_{run_id}/OutputDir"
    file_stub = gc_date.strftime("GEOSChem.SpeciesConc.%Y%m%d_0000z.nc4")
    filepath = os.path.join(j_dir, file_stub)

    # Construct the list of species vars to request
    keepvars = [f"SpeciesConcVV_{config['Species']}_{i:04}" for i in sv_elems]
    if len(keepvars) == 1:
        is_Regional = config["isRegional"]
        is_OH_element = check_is_OH_element(
            sv_elems[0], n_elements, config["OptimizeOH"], is_Regional
        )
        is_BC_element = check_is_BC_element(
            sv_elems[0],
            n_elements,
            config["OptimizeOH"],
            config["OptimizeBCs"],
            is_OH_element,
            is_Regional,
        )
        if is_OH_element or is_BC_element:
            keepvars = [f"SpeciesConcVV_{config['Species']}"]

    if baserun:
        keepvars = [f"SpeciesConcVV_{config['Species']}"]

    # It would fail if open all variables with chunks with GCHP,
    # as ncontact is duplicate for GCHP output dimensions
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(filepath, decode_cf=False) as tmp:
            other_vars = [v for v in tmp.variables if f"SpeciesConcVV_{config['Species']}" not in v]

    # Open only these variables
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="xarray")
        with xr.open_dataset(
            filepath,
            drop_variables=other_vars,
            chunks="auto"
        ) as dsmf_all:
            try:
                if dsmf_all.sizes.get("time", 0) == 0:
                    print(f"ERROR: {filepath}: empty time dimension", flush=True)
                if config['UseGCHP']:
                    nfi   = xr.DataArray(gridcell_dict["nfi"],   dims="obs")
                    Ydimi = xr.DataArray(gridcell_dict["Ydimi"], dims="obs")
                    Xdimi = xr.DataArray(gridcell_dict["Xdimi"], dims="obs")
                    dsmf = dsmf_all.isel(
                        time=gc_date.hour).squeeze().isel(
                        nf=nfi,
                        Ydim=Ydimi,
                        Xdim=Xdimi,
                        drop=True
                    )
                else:
                    jGC   = xr.DataArray(gridcell_dict["jGC"],   dims="obs")
                    iGC   = xr.DataArray(gridcell_dict["iGC"],   dims="obs")
                    dsmf = dsmf_all.isel(
                        time=gc_date.hour).squeeze().isel(
                        lat=jGC,
                        lon=iGC,
                        drop=True
                    )
            except Exception as e:
                print(f"Run id {run_id}. Failed at {gc_date} with error: {e}", flush=True)
                raise

            # ---- Batch read all species and standardize to (elem, obs, lev)
            da_list = []
            for v in keepvars:
                da = dsmf[v]
                # Ensure shape is (obs, lev)
                arr = da.transpose("obs","lev").values
                da_list.append(arr)
            species_all = np.stack(da_list, axis=0)   # (elem, obs, lev)

            dry_air_subcolumns = gridcell_dict["dry_air_subcolumns"]  # (N, S)
            apriori = gridcell_dict["apriori"]                        # (N, S)
            avkern = gridcell_dict["avkern"]                          # (N, S)
            denom = np.sum(dry_air_subcolumns, axis=1)                # (N,)

            # ---- Remap levels to satellite layers in batch:
            # vertical_weights:   (N, S, G)
            # species_all:        (E, N, G)
            # -> sat_species_all: (E, N, S)
            sat_species_all = np.einsum("nsg,eng->ens", vertical_weights, species_all)

            # Convert to column units and apply AKs, batched over E
            sat_species_molm2_all = sat_species_all * dry_air_subcolumns[None, :, :]  # (E, N, S)
            numer_all = np.sum(apriori[None, :, :] +
                               avkern[None, :, :] * (sat_species_molm2_all - apriori[None, :, :]),
                               axis=2)  # (E, N)
            virtual_satellite_all = (numer_all / denom[None, :]).T.astype(np.float32)  # (N, E)

    return virtual_satellite_all # unitless mixing ratio

def apply_averaging_kernel(
        apriori,
        avkern,
        sat_species_molm2,
        dry_air_subcolumns,
        species
):
    # Derive the column-averaged mixing ratio that the satellite would see
    # over this ground cell using eq. 46 from TROPOMI Methane ATBD,
    # Hasekamp et al. 2019
    virtual_satellite = (
        sum(apriori + avkern * (sat_species_molm2 - apriori))
        / sum(dry_air_subcolumns)
        * mixing_ratio_conv_factor(species)
    )  # volumetric mixing ratio
    return virtual_satellite
