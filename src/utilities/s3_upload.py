import os
import sys
import glob
import yaml
import tarfile
import boto3

"""
A simple utility script that uploads specified files to s3 as a tar archive.
Arguments
    config_path   [String]   : path to yaml config file
"""


def extract_s3_part(s3_path, s3_part):
    """
    Extracts the bucket or key from an s3 path
    Arguments
        s3_path  [str] : s3 path
        s3_part  [str] : "bucket" or "key"

    Returns
        part     [str] : processed s3 bucket or key
    """
    # Remove the s3:// prefix if present
    s3_path = s3_path.replace("s3://", "")

    # Split the path into bucket and key
    parts = s3_path.split("/", 1)

    # Return the bucket name
    if s3_part == "bucket":
        return parts[0]
    elif s3_part == "key":
        return parts[1]
    else:
        # Print error message and exit if part is not bucket or key
        sys.exit(f'Invalid s3 part specified: {s3_part}. Must be "bucket" or "key".')


def zip_and_upload_to_s3(file_paths, bucket_name, s3_key):
    """
    Zips specified files into a tar archive and uploads it to the specified s3 path.
    Arguments
        file_paths      [] : s3 paths to files to be archived and uploaded
        bucket_name  [str] : s3 bucket name
        s3_key       [str] : s3 key
    """

    # Create a temporary zip file
    tar_file_name = "temp.tar.gz"
    expanded_file_paths = []

    # use glob to expand wildcards for file path
    for file_path in file_paths:
        expanded_file_paths.extend(glob.glob(file_path))

    if len(expanded_file_paths) == 0:
        sys.exit(
            f"No files found to archive/upload at the specified path(s): {file_paths}"
        )

    print(f"Archiving specified output files into {tar_file_name}...")
    with tarfile.open(tar_file_name, "w:gz") as tar_file:
        # Add each file to the tar file
        for file_path in expanded_file_paths:
            tar_file.add(file_path, arcname=os.path.basename(file_path))

    print(f"Uploading tar file to {s3_key}...")
    # Upload the tar file to S3
    s3_client = boto3.client("s3")
    with open(tar_file_name, "rb") as file_data:
        s3_client.upload_fileobj(file_data, bucket_name, s3_key)

    # Delete the temporary zip file
    os.remove(tar_file_name)
    print(f"Upload Complete and temporary tar file removed.")


if __name__ == "__main__":
    config_path = sys.argv[1]
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    bucket = extract_s3_part(config["S3UploadPath"], "bucket")
    key = extract_s3_part(config["S3UploadPath"], "key")

    # add / if key does not already end in /
    if key[-1] != "/":
        key = key + "/" + config["RunName"] + ".tar.gz"
    else:
        key = key + config["RunName"] + ".tar.gz"

    zip_and_upload_to_s3(config["OutputFiles"], bucket, key)
