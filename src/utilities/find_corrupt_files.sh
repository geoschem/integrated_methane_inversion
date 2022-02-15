#!/bin/bash

# Arguments:
#   s3_path: path of s3 files to recursively download
#   dest: path files will be downloaded to
download_aws_files() {
    aws s3 cp $1 $2 --request-payer --recursive --dryrun
}

report() {
  corrupt=1
} >&2

# Configuration
S3_BUCKET=gcgrid
RESOLUTION=GEOS_0.25x0.3125_AS
MET=GEOS_FP
YEAR=2019
DEST="/home/ubuntu/ExtData/${RESOLUTION}/${MET}/${YEAR}/"
S3_PATH="s3://${S3_BUCKET}/${RESOLUTION}/${MET}/${YEAR}/"

mkdir -p $DEST
download_aws_files $S3_PATH $DEST

corrupt=0
trap report ERR

FILE_LIST=$(find "${DEST}" -type f)
echo $DEST
for file in $FILE_LIST; do
    echo "testing file: ${file}"
    nccopy $file testfile.nc
    if [ $corrupt -gt 0 ]; then  
        echo "$file is corrupt"
        echo $file >> "corrupted_files${RESOLUTION}_${YEAR}.txt"
        corrupt=0
    fi
    rm -f testfile.nc
done



