#!/usr/bin/env python3

"""
Description:
------------
This script reads a GEOS-Chem dry-run log and downloads missing data
files from the S3 bucket `s3://gcgrid` using boto3.

This script is based on the download_data.py script in GEOS-Chem, but
it has been modified to download files from S3 using boto3 instead of
using the `aws s3 cp` command. Boto3 is much faster than the aws cli.

Requirements:
-------------
- boto3
- PyYAML
"""

import os
import sys
import yaml
import boto3
import re
from botocore.exceptions import ClientError

# Globals
GEOSCHEM_INPUT_FILE = "./geoschem_config.yml"
S3_BUCKET = "gcgrid"


def read_config_file(config_file, to_str=False):
    try:
        with open(config_file, encoding="UTF-8") as stream:
            if to_str:
                return yaml.load(stream, Loader=yaml.loader.BaseLoader)
            return yaml.load(stream, Loader=yaml.loader.SafeLoader)
    except FileNotFoundError as err:
        raise FileNotFoundError(f"Error reading {config_file}: {err}") from err


def extract_pathnames_from_log(args):
    comments = ["!" * 79, "!!! LIST OF (UNIQUE) FILES REQUIRED FOR THE SIMULATION"]
    data_found = set()
    data_missing = set()
    dryrun_log = args["dryrun_log"]

    with open(dryrun_log, "r", encoding="UTF-8") as f:
        for line in f:
            line = line.replace("CHEM_INPUTS//", "CHEM_INPUTS/")
            upcaseline = line.upper()
            if ": OPENING" in upcaseline or ": READING" in upcaseline:
                data_found.add(line.split()[-1])
            elif "FILE NOT FOUND" in upcaseline:
                data_missing.add(line.split()[-1])
            elif any(
                tag in upcaseline
                for tag in ["!!! STA", "!!! END", "!!! SIM", "!!! MET", "!!! GRI"]
            ):
                comments.append(line.rstrip())

    comments.append("!" * 79)
    found = sorted(data_found)
    missing = sorted(data_missing)

    local_prefix = ""
    for path in found + missing:
        if "ExtData" in path:
            index = path.find("ExtData")
            local_prefix = path[:index]
            break

    if not local_prefix:
        raise ValueError("Could not locate the ExtData folder in local paths.")

    return {
        "comments": comments,
        "found": found,
        "missing": missing,
        "local_prefix": local_prefix,
    }


def get_run_info():
    config = read_config_file(GEOSCHEM_INPUT_FILE, to_str=True)
    run_info = {
        "start_date": int(config["simulation"]["start_date"][0]),
        "start_time": int(config["simulation"]["start_date"][1]),
        "end_date": int(config["simulation"]["end_date"][0]),
        "end_time": int(config["simulation"]["end_date"][1]),
        "met_field": config["simulation"]["met_field"],
        "sim": config["simulation"]["name"],
        "resolution": config["grid"]["resolution"],
        "tomas15": "NK15" in config["operations"]["transport"]["transported_species"],
        "tomas40": "NK40" in config["operations"]["transport"]["transported_species"],
    }
    return run_info


def write_unique_paths(paths, unique_log):
    combined = sorted(paths["found"] + paths["missing"])
    try:
        with open(unique_log, "w", encoding="UTF-8") as f:
            for comment in paths["comments"]:
                print(comment, file=f)
            for path in combined:
                print(path, file=f)
            for comment in paths["comments"]:
                print(comment, file=f)
        print(f"Log with unique file paths written to: {unique_log}")
    except Exception as e:
        raise RuntimeError(f"Could not write to {unique_log}") from e


def download_file_from_s3(s3_client, bucket_name, s3_key, local_path):
    try:
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        s3_client.download_file(bucket_name, s3_key, local_path)
        print(f"Downloaded: s3://{bucket_name}/{s3_key} -> {local_path}")
    except ClientError as e:
        print(f"ERROR: Failed to download s3://{bucket_name}/{s3_key}: {e}")
        raise


def download_the_data(args):
    run_info = get_run_info()
    paths = extract_pathnames_from_log(args)
    write_unique_paths(paths, args["dryrun_log"] + ".unique")

    if args["skip_download"]:
        return

    print(f"Downloading data from S3 bucket: {S3_BUCKET}")
    s3 = boto3.client("s3")

    for path in paths["missing"]:
        if "ExtData" not in path:
            continue

        local_path = path.split("-->")[0].strip() if "-->" in path else path

        match = re.search(r"ExtData/(.+)", path)
        if not match:
            print(f"WARNING: Could not extract S3 key from: {path}")
            continue
        s3_key = match.group(1)

        substitutions = {
            "IPMN": "PMN",
            "NPMN": "PMN",
            "RIPA": "RIP",
            "RIPB": "RIP",
            "RIPD": "RIP",
        }
        filename = os.path.basename(s3_key)
        for wrong, correct in substitutions.items():
            if wrong in filename:
                corrected_key = s3_key.replace(wrong, correct)
                temp_path = os.path.join(
                    os.path.dirname(local_path), filename.replace(wrong, correct)
                )
                download_file_from_s3(s3, S3_BUCKET, corrected_key, temp_path)
                os.rename(temp_path, local_path)
                break
        else:
            download_file_from_s3(s3, S3_BUCKET, s3_key, local_path)

    chem_dir = os.path.join(paths["local_prefix"], "CHEM_INPUTS")
    os.makedirs(chem_dir, exist_ok=True)


def parse_args():
    dryrun_log = None
    skip_download = False

    for arg in sys.argv[1:]:
        if dryrun_log is None and not arg.startswith("-"):
            dryrun_log = arg
        elif "skip" in arg.lower():
            skip_download = True

    if dryrun_log is None:
        raise ValueError("You must provide a dryrun log file as the first argument.")

    return {
        "dryrun_log": dryrun_log,
        "skip_download": skip_download,
    }


def main():
    download_the_data(parse_args())


if __name__ == "__main__":
    main()
