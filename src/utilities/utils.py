import os
import boto3
from botocore import UNSIGNED
from botocore.client import Config


def download_landcover_files(config):
    """
    Download landcover files from s3 given the config file
    """
    # Initialize s3 service with boto3
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    bucket = "gcgrid"

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
        gridDir = "0.25x0.3125"
        gridFile = "025x03125"

    # determine the path to the landcover file
    if len(config["RegionID"]) == 2:
        s3_lc_path = f"GEOS_{gridDir}_{config['RegionID']}/{metDir}/{constYr}/01/{config['Met']}.{constYr}0101.CN.{gridFile}.{config['RegionID']}.{LandCoverFileExtension}"
    else:
        s3_lc_path = f"GEOS_{gridDir}/{metDir}/{constYr}/01/{config['Met']}.{constYr}0101.CN.{gridFile}.{LandCoverFileExtension}"
    LandCoverFile = os.path.join(config["DataPath"], s3_lc_path)
    target_dir = os.path.dirname(LandCoverFile)

    # Check if the file already exists locally
    if not os.path.exists(LandCoverFile):
        # If not, download it
        os.makedirs(target_dir, exist_ok=True)
        s3.download_file(bucket, s3_lc_path, LandCoverFile)
        print(f"File {LandCoverFile} downloaded successfully.")
    else:
        print(f"File {LandCoverFile} already exists locally. Skipping download.")


def download_hemcodiags_files(config):
    """
    Download global hemco diagnostics files from s3 given the config file
    """
    # Initialize s3 service with boto3
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    bucket = "gcgrid"

    if config["Res"] == "4.0x5.0":
        gridFile = "4x5"
    elif config["Res"] == "2.0x2.5":
        gridFile = "2x25"
    elif config["Res"] == "0.5x0.625":
        gridFile = "05x0625"
    elif config["Res"] == "0.25x0.3125":
        gridFile = "025x03125"

    s3_hd_path = (
        f"HEMCO/CH4/v2024-07/HEMCO_SA_Output/HEMCO_sa_diagnostics.{gridFile}.2023.nc"
    )
    HemcoDiagFile = os.path.join(config["DataPath"], s3_hd_path)
    target_dir = os.path.dirname(HemcoDiagFile)

    # Check if the file already exists locally
    if not os.path.exists(HemcoDiagFile):
        # If not, download it
        os.makedirs(target_dir, exist_ok=True)
        s3.download_file(bucket, s3_hd_path, HemcoDiagFile)
        print(f"File {HemcoDiagFile} downloaded successfully.")
    else:
        print(f"File {HemcoDiagFile} already exists locally. Skipping download.")
