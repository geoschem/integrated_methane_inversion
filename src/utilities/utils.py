import subprocess


def download_landcover_files(config):
    """
    Download landcover files from s3 given the config file
    """
    DataPath = "/home/ubuntu/ExtData"

    if config["Met"] == "GEOSFP":
        metDir = "GEOS_FP"
        constYr = "2011"
        LandCoverFileExtension = "nc"
    elif config["Met"] == "MERRA2":
        metDir = "MERRA2"
        constYr = "2015"
        LandCoverFileExtension = "nc4"

    if config["Res"] == "4.0x5.0":
        gridDir = "4x5"
        gridFile = "4x5"
    elif config["Res"] == "2.0x2.5":
        gridDir = "2x2.5"
        gridFile = "2x25"
    elif config["Res"] == "0.5x0.625":
        gridDir = "0.5x0.625"
        gridFile = "05x0625"
    elif config["Res"] == "0.25x0.3125":
        gridDir= "0.25x0.3125"
        gridFile= "025x03125"
       
    if len(config["RegionID"]) == 2:
        LandCoverFile = f"{DataPath}/GEOS_{gridDir}_{config['RegionID']}/{metDir}/{constYr}/01/{config['Met']}.{constYr}0101.CN.{gridFile}.{config['RegionID']}.{LandCoverFileExtension}"
        s3_lc_path = f"s3://gcgrid/GEOS_{gridDir}_{config['RegionID']}/{metDir}/{constYr}/01/{config['Met']}.{constYr}0101.CN.{gridFile}.{config['RegionID']}.{LandCoverFileExtension}"
    else:
        LandCoverFile = f"{DataPath}/GEOS_{gridDir}/{metDir}/{constYr}/01/{config['Met']}.{constYr}0101.CN.{gridFile}.{LandCoverFileExtension}"
        s3_lc_path = f"s3://gcgrid/GEOS_{gridDir}/{metDir}/{constYr}/01/{config['Met']}.{constYr}0101.CN.{gridFile}.{LandCoverFileExtension}"

    # run the aws command to download the files
    command = f"aws s3 cp --no-sign-request {s3_lc_path} {LandCoverFile}"
    results = subprocess.run(command.split(), capture_output=True, text=True)

    output = (
        "Successfully downloaded landcover files"
        if results.returncode == 0
        else results
    )
    print(output)


def download_hemcodiags_files(config):
    """
    Download global hemco diagnostics files from s3 given the config file
    """
    DataPath = "/home/ubuntu/ExtData"

    if config["Res"] == "4.0x5.0":
        gridFile = "4x5"
    elif config["Res"] == "2.0x2.5":
        gridFile = "2x25"
    elif config["Res"] == "0.5x0.625":
        gridFile = "05x0625"
    elif config["Res"] == "0.25x0.3125":
        gridFile = "025x03125"

    HemcoDiagFile = f"{DataPath}/HEMCO/CH4/v2023-04/HEMCO_SA_Output/HEMCO_sa_diagnostics.{gridFile}.20190101.nc"
    s3_hd_path = f"s3://gcgrid/HEMCO/CH4/v2023-04/HEMCO_SA_Output/HEMCO_sa_diagnostics.{gridFile}.20190101.nc"

    # run the aws command to download the files
    command = f"aws s3 cp --no-sign-request {s3_hd_path} {HemcoDiagFile}"
    results = subprocess.run(command.split(), capture_output=True, text=True)

    output = (
        "Successfully downloaded hemco diags files"
        if results.returncode == 0
        else results
    )
    print(output)
