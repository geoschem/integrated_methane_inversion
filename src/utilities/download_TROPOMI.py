#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import boto3
import pandas as pd
import multiprocessing
import datetime
from botocore import UNSIGNED
from botocore.client import Config

# Description: Download TROPOMI data from meeo S3 bucket for desired dates.
#              Function can be called from another script or run as a
#              directly as a script.
# Example Usage as a script:
#     $ python download_TROPOMI.py 20190101 20190214 TROPOMI_data

s3 = None
VALID_TROPOMI_PROCESSOR_VERSIONS = ["020400", "020500", "020600", "020701", "020800"]


def initialize_boto3():
    """
    Initialize s3 service with boto3
    """
    global s3
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))


def download_from_s3(args):
    """
    Download s3 path to local directory if it doesn't already exist locally
    Arguments
        args             [tuple] : (s3_path, bucket, storage_dir)
    """
    s3_path, bucket, storage_dir = args
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    file = os.path.basename(s3_path)
    local_file_path = os.path.join(storage_dir, file)

    # Check if the file already exists locally
    if not os.path.exists(local_file_path):
        # If not, download it
        s3.download_file(bucket, s3_path, local_file_path)
    else:
        print(f"File {local_file_path} already exists locally. Skipping download.")


def get_s3_paths(start_date, end_date, bucket):
    """
    Gets s3 paths for download.
    Arguments
        start_date  [datetime] : start date of data download (yyyymmdd)
        end_date    [datetime] : end date of data download (yyyymmdd)
        bucket           [str] : s3 bucket name
    Returns
        s3_paths   [list] : list of s3 paths for download
    """

    # Make a list of all possible dates in our date range
    years, months, days = [], [], []
    current_date = start_date
    while current_date < end_date:
        years.append(current_date.strftime("%Y"))
        months.append(current_date.strftime("%m"))
        days.append(current_date.strftime("%d"))
        current_date += datetime.timedelta(days=1)

    # Get access to S3
    initialize_boto3()

    # Collect information about each of the files
    # Need to loop through OFFL and RPRO directories
    # We'll get a ton of duplicate files here
    names, prefixes = [], []
    for type in ["RPRO", "OFFL"]:
        for i in range(len(days)):
            Prefix = f"{type}/L2__CH4___/" f"{years[i]}/{months[i]}/{days[i]}/"
            response = s3.list_objects(Bucket=bucket, Prefix=Prefix)
            if "Contents" in response:
                for key in s3.list_objects(Bucket=bucket, Prefix=Prefix)["Contents"]:
                    prefixes.append(Prefix)
                    names.append(key["Key"].split("/")[-1])

    # Organize the information about each of the files
    df = pd.DataFrame()
    df["Name"] = names
    df["Prefix"] = prefixes
    df["ProcessorVersion"] = df["Name"].str.extract(r"_(\d{6})_")
    df["ProcessingMode"] = df["Name"].str.extract(r"_(\S{4})_")
    df["OrbitNumber"] = df["Name"].str.extract(r"_(\d{5})_")
    df["ModificationDate"] = df["Name"].str.extract(r"_\d{6}_(\d{8}T\d{6})")
    df["CollectionNumber"] = df["Name"].str.extract(r"_(\d{2})_")

    # We only want files that are v02.04.00, v02.05.00, or v02.06.00.
    # Also make sure the collection number is 03 to account for some duplicates.
    df.loc[df["ProcessorVersion"].isin(VALID_TROPOMI_PROCESSOR_VERSIONS)]
    df = df.drop_duplicates(subset=["Name", "ModificationDate"])
    df = df.loc[df["CollectionNumber"] == "03"].reset_index(drop=True)

    # Deal with duplicate orbit numbers.
    unique_orbit_numbers = df["OrbitNumber"].unique()
    for unique_orbit_number in unique_orbit_numbers:
        subset = df.loc[df["OrbitNumber"] == unique_orbit_number]
        # If this id is already a unique orbit, continue.
        if len(subset) == 1:
            continue
        # If there is both RPRO and OFFL, keep RPRO.
        elif len(subset["ProcessingMode"].unique()) > 1:
            index_to_drop = subset.loc[subset["ProcessingMode"] != "RPRO"].index
            df = df.drop(index_to_drop)
    
    # If there are multiple processor versions for the same orbit number,
    # keep the one with the latest version.
    df = df.loc[df.groupby("OrbitNumber")["ProcessorVersion"].idxmax()]
    df.reset_index(drop=True, inplace=True)

    # For the final duplicates, choose the later modification date
    df["ModificationDate"] = pd.to_datetime(df["ModificationDate"])
    df = df.loc[df.groupby("OrbitNumber")["ModificationDate"].idxmax()]
    df.reset_index(drop=True, inplace=True)
        
    assert len(df) == len(df["OrbitNumber"].unique())
    df = df.reset_index(drop=True)
    df["S3Path"] = df["Prefix"] + df["Name"]
    s3_paths = sorted(df["S3Path"].to_list())

    return s3_paths


def download_operational_TROPOMI(start_date, end_date, storage_dir):
    """
    Download TROPOMI data from AWS for desired dates.

    Arguments
        start_date  [datetime] : start date of data download
        end_date    [datetime] : end date of data download
        storage_dir      [str] : local directory to store downloaded files
    """

    bucket = "meeo-s5p"
    s3_paths = get_s3_paths(start_date, end_date, bucket)
    os.makedirs(storage_dir, exist_ok=True)

    print("=============Downloading TROPOMI Operational Data=============")
    print(f"Downloading {len(s3_paths)} files - ({start_date},{end_date}].")
    # Download the files using multiple cores
    with multiprocessing.Pool(112) as pool:
        pool.map(
            download_from_s3, [(s3_path, bucket, storage_dir) for s3_path in s3_paths]
        )
        pool.close()
        pool.join()
    print("==================Finished TROPOMI Download ==================")

    return s3_paths


if __name__ == "__main__":
    start = sys.argv[1]
    end = sys.argv[2]
    Sat_datadir = sys.argv[3]
    start_date = datetime.datetime.strptime(start, "%Y%m%d")
    end_date = datetime.datetime.strptime(end, "%Y%m%d")
    download_operational_TROPOMI(start_date, end_date, Sat_datadir)
