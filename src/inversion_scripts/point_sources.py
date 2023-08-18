# python function related to point sources

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
            coords = read_point_source_csv(path)
        # handle list of lists
        elif isinstance(coord_var, list):
            coords = coord_var
        else:
            # Variable is neither a string nor a list
            print("Warning: No ForcedNativeResolutionElements specified or invalid format.")
    
    # then we read point sources from external datasources
    plumes = []
    if "PointSourceDatasets" in config.keys():
        if "SRON" in config["PointSourceDatasets"]:
            print("Fetching plumes from SRON database...")
            plumes = SRON_plumes(config)
        else:
            print('No valid external point source datasets specified. Valid values are: "SRON"')
        
    # append SRON
    return coords + plumes