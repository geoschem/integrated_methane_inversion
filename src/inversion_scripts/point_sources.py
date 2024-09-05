# This file contains functions related to accessing point source data
import os
import datetime
from copy import copy
import pandas as pd
import geopandas as gpd
from bs4 import BeautifulSoup
from shapely.geometry import Point
from dateutil.relativedelta import relativedelta
from src.inversion_scripts.plumes import GeoFilter, PointSources
import src.inversion_scripts.plumes as imiplumes


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

        # List of datasources
        ps_datasets = copy(config["PointSourceDatasets"])

        if len(ps_datasets) > 0:

            observers = []
            for plume_dataset in ps_datasets:
                # instantiate plume observer object
                # using its string name from config file
                observer = getattr(imiplumes, plume_dataset)
                observers.append(observer(config, usecached=True))

            gf = GeoFilter(config)
            ps = PointSources(gf, observers)
            plumes += ps.get_gridded_coords(
                emission_rate_filter=int(config["EmissionRateFilter"]),
                plume_count_filter=int(config["PlumeCountFilter"]),
            )
            msg = (
                f"Found {len(plumes)} grid cells with plumes "
                f"using {[ob.myname for ob in observers]}"
            )
            print(msg, flush=True)
            got_plumes = True

        if not got_plumes:
            print(
                "No valid external point source datasets specified. "
                'Valid values are: "SRON", "CarbonMapper", "IMEO"'
            )

    # append point sources
    coords = coords + plumes

    return coords
