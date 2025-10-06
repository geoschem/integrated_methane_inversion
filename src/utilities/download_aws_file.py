import sys
import os
import boto3
from botocore import UNSIGNED
from botocore.client import Config
"""
A simple utility script that downloads a file from an S3 bucket to a local directory.
Usage:
    python download_aws_file.py s3://bucket-name/path/to/file local/file/path

"""

def download_from_s3(s3_path, bucket, local_file_path):
    """
    Download s3 path to local directory if it doesn't already exist locally
    Arguments
        s3_path             [str] : path excluding bucket name
        bucket              [str] : name of the bucket
        local_file_path     [str] : path to save the file locally
    """
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    # Check if the file already exists locally
    if not os.path.exists(local_file_path):
        # Create the directory if it doesn't exist
        os.makedirs(os.path.dirname(local_file_path), exist_ok=True)
        # Download the file from S3
        s3.download_file(bucket, s3_path, local_file_path)
    else:
        print(f"File {local_file_path} already exists locally. Skipping download.")

if __name__ == "__main__":
    aws_file = sys.argv[1]
    file_des = sys.argv[2]
    aws_file = aws_file.replace("s3://", "")
    bucket_name = aws_file.split("/")[0]
    aws_file = "/".join(aws_file.split("/")[1:])
    
    download_from_s3(aws_file, bucket_name, file_des)
