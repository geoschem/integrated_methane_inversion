#!/usr/bin/env python3
import os
import sys
import boto3
import multiprocessing
from botocore import UNSIGNED
from botocore.client import Config

"""
Description: Download boundary condition data from aws S3 bucket for 
             desired dates. Function can be called from another script 
             or run as a directly as a script.
Example Usage:
    $ python download_bc.py 20200501 20200508 path/to/imi/boundary-conditions v2025-06
"""


def list_missing_files(start_date, end_date, destination):
    """
    List missing files in the specified date range.
    Arguments:
        start_date       [str] : start date of data download (yyyymmdd)
        end_date         [str] : end date of data download (yyyymmdd)
        destination      [str] : local directory to store downloaded files
    """

    missing_files = []

    start_str = str(start_date)
    start_year = start_str[:4]
    start_month = start_str[4:6]
    start_day = start_str[6:8]
    end_str = str(end_date)
    end_year = end_str[:4]
    end_month = end_str[4:6]
    end_day = end_str[6:8]

    month_days = [31, [28, 29], 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    file_prefix = "GEOSChem.BoundaryConditions."
    file_suffix = "_0000z.nc4"

    for year in range(int(start_year), int(end_year) + 1):
        # skip years with definite no data
        if year < 2018:
            print(
                "Skipping BC data download for ", str(year), ": no data from this year"
            )
            continue
        init_month = 1
        final_month = 12
        if year == int(start_year):
            # only get desired months from incomplete years
            init_month = int(start_month)
        if year == int(end_year):
            final_month = int(end_month)
        for month in range(init_month, final_month + 1):
            # skip months with definite no data
            if year == 2018 and month < 4:
                print(
                    "Skipping BC data download for ",
                    str(year),
                    "/0",
                    str(month),
                    ": no data from this month",
                )
                continue
            # add 0 to month string if necessary
            month_prefix = "0" if month < 10 else ""
            init_day = 1
            final_day = month_days[month - 1]
            # leap day
            if month == 2:
                if year % 4 == 0:
                    final_day = final_day[1]
                else:
                    final_day = final_day[0]
            if month == int(start_month) and year == int(start_year):
                # only get desired days from incomplete months
                init_day = int(start_day)
            if month == int(end_month) and year == int(end_year):
                final_day = int(end_day)
            for day in range(init_day, final_day + 1):
                # add 0 to day string if necessary
                day_prefix = "0" if day < 10 else ""

                # check if file for this day already exists
                file_name = (
                    file_prefix
                    + str(year)
                    + month_prefix
                    + str(month)
                    + day_prefix
                    + str(day)
                    + file_suffix
                )
                # add file to download list if needed
                if not os.path.exists(destination + "/" + file_name):
                    missing_files.append(file_name)

    return missing_files


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


def download_data(start_date, end_date, storage_dir, version):
    """
    Download boundary condition data from aws S3 bucket for
    desired dates. Function can be called from another script
    or run as a directly as a script.
    Arguments:
        start_date       [str] : start date of data download (yyyymmdd)
        end_date         [str] : end date of data download (yyyymmdd)
        storage_dir      [str] : local directory to store downloaded files
        version          [str] : version of boundary conditions to download
    """
    bucket = "imi-boundary-conditions"
    s3_paths = list_missing_files(start_date, end_date, storage_dir)
    s3_paths = [os.path.join(version, path) for path in s3_paths]

    os.makedirs(storage_dir, exist_ok=True)

    print(f"=============Downloading Boundary Conditions {version}=============")
    print(f"Downloading {len(s3_paths)} files - [{start_date},{end_date}].")
    # Download the files using multiple cores
    with multiprocessing.Pool(112) as pool:
        pool.map(
            download_from_s3, [(s3_path, bucket, storage_dir) for s3_path in s3_paths]
        )
        pool.close()
        pool.join()


if __name__ == "__main__":
    download_data(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
