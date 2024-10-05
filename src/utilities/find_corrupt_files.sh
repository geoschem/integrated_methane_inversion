#!/bin/bash

# Usage: ./find_corrupt_files.sh > corrupt.log
# This shell script is used to verify the integrity of
# netcdf s3 datasets. It is meant to download the
# specified data from s3 and run an nccopy command to
# verify the file is in the correct format. All corrupted
# files are written out to a txt file

# ************** Configuration **************************
S3_BUCKET=gcgrid
RESOLUTION=GEOS_0.25x0.3125_AS
MET=GEOS_FP
YEAR=2019
DEST="/home/ubuntu/ExtData/${RESOLUTION}/${MET}/${YEAR}/"
S3_PATH="s3://${S3_BUCKET}/${RESOLUTION}/${MET}/${YEAR}/"
# ************** Configuration End **********************

# Usage: download specified s3 data recursively
# Arguments:
#   s3_path: path of s3 files to recursively download
#   dest: path files will be downloaded to
download_aws_files() {
    echo "starting aws download"
    aws s3 cp $1 $2 --no-sign-request --recursive
    echo "finished aws download"
}

# Usage: report errors from the corruption test
report() {
    corrupt=1
} >&2

mkdir -p $DEST
download_aws_files $S3_PATH $DEST

# trap will automatically send any error exit codes to the
# report function which sets the corruption boolean to 1 (True)
corrupt=0
trap report ERR

# Get a sorted list of files from the destination directory
FILE_LIST=$(find "${DEST}" -type f | sort)
echo "starting file corruption testing"

# Iterate through the list and test whether each file successfully
# performs the nccopy command
for file in $FILE_LIST; do
    echo "testing file: ${file}"
    nccopy $file testfile.nc
    if [ $corrupt -gt 0 ]; then
        echo "$file is corrupt"
        echo $file >>"corrupted_files${RESOLUTION}_${YEAR}.txt"
        corrupt=0
    fi
    rm -f testfile.nc
done
echo "Finished file corruption testing for specified data"
