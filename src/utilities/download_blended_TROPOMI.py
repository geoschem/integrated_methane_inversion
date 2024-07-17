#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import boto3
import multiprocessing
from datetime import datetime, timedelta
from botocore import UNSIGNED
from botocore.client import Config

# Description: Download TROPOMI data from aws S3 bucket for desired dates.
#              Function can be called from another script or run as a
#              directly as a script.
# Example Usage as a script:
#     $ python download_blended_TROPOMI.py 20190101 20190214 TROPOMI_data


def initialize_boto3():
    """
    Initialize s3 service with boto3
    Returns
        s3   [object] : boto3 s3 client
    """
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    return s3


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


def get_file_prefixes(start_date, end_date):
    """
    Gets month prefixes and file prefixes for
    supplied date range (excluding end_date)
    Arguments
        start_date  [datetime] : start date of data download
        end_date    [datetime] : end date of data download
    Returns
        months_list   [list] : prefix dirs of all months
        days_list     [list] : acceptable file path prefixes
    """

    # default prefixes of blended dataset bucket
    subdir = "data/"
    file_prefix = "S5P_BLND_L2__CH4____"

    # Initialize set to store month strings
    # using a set prevents duplicates
    months_set = set()
    days_set = set()

    # Iterate through months between start and end dates
    # use 27 days to prevent missing a month
    current_date = start_date
    while current_date < end_date:
        month_dir = current_date.strftime("%Y-%m") + "/"
        months_set.add(subdir + month_dir)
        days_set.add(subdir + month_dir + file_prefix + current_date.strftime("%Y%m%d"))
        current_date += timedelta(days=1)

    return sorted(months_set), sorted(days_set)


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
    month_prefix_list, file_prefix_list = get_file_prefixes(start_date, end_date)
    s3 = initialize_boto3()
    s3_paths = []
    for month_prefix in month_prefix_list:
        response = s3.list_objects(Bucket=bucket, Prefix=month_prefix)
        if "Contents" in response:
            for key in response["Contents"]:
                day_match = any(
                    key["Key"].startswith(prefix) for prefix in file_prefix_list
                )
                if day_match:
                    s3_paths.append(key["Key"])
        else:
            raise Exception(
                f"s3://{bucket}/{month_prefix} does not exist. Check if blended TROPOMI+GOSAT data exists for time period."
            )
    return s3_paths


def download_blended(start_date, end_date, storage_dir):
    """
    Download blended TROPOM+GOSAT dataset from s3 to desired
    directory.
    Arguments
        start_date  [datetime] : start date of data download (yyyymmdd)
        end_date    [datetime] : end date of data download (yyyymmdd)
        storage_dir      [str] : local directory to store downloaded files
    """
    bucket = "blended-tropomi-gosat-methane"
    s3_paths = get_s3_paths(start_date, end_date, bucket)
    os.makedirs(storage_dir, exist_ok=True)

    print("=============Downloading Blended TROPOMI+GOSAT Data=============")
    print(f"Downloading {len(s3_paths)} files - ({start_date},{end_date}].")
    # Download the files using multiple cores
    with multiprocessing.Pool(112) as pool:
        pool.map(
            download_from_s3, [(s3_path, bucket, storage_dir) for s3_path in s3_paths]
        )
        pool.close()
        pool.join()

    print("==============Finished Downloading Blended Dataset==============")


if __name__ == "__main__":
    start = sys.argv[1]
    end = sys.argv[2]
    Sat_datadir = sys.argv[3]
    start_date = datetime.strptime(start, "%Y%m%d")
    end_date = datetime.strptime(end, "%Y%m%d")
    download_blended(start_date, end_date, Sat_datadir)
