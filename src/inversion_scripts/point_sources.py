# This file contains functions related to accessing point source data
import os
import datetime
import requests
import pandas as pd
import geopandas as gpd
from bs4 import BeautifulSoup
from shapely.geometry import Point
from dateutil.relativedelta import relativedelta


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
    if "PointSourceDatasets" in config.keys():
        if "SRON" in config["PointSourceDatasets"]:
            print("Fetching plumes from SRON database...")
            plumes = SRON_plumes(config)
        else:
            print(
                'No valid external point source datasets specified. Valid values are: "SRON"'
            )

    # append SRON
    return coords + plumes


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

    plumes_list = plumes[["lat", "lon"]].values.tolist()
    return plumes_list
