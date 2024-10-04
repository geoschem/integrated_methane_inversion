#!/bin/bash

# This script will download 0.25x0.3125 degree global GEOS-FP meteorogical
# files in netCDF format and crop them to a user-defined region. Once files
# have been cropped, the global files may be removed to free up disk space
# using the RemoveGlobalMet option.

##############################################################################
# Custom to Harvard FAS RC cluster:
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-6:00
#SBATCH -p huce_cascade
#SBATCH --mem=2000
#SBATCH --mail-type=END

# Load modules for CDO
module load intel/17.0.4-fasrc01
module load cdo/1.9.4-fasrc02
##############################################################################

##############################################################################
# USER SETTINGS:
##############################################################################

# Create a 2-letter region ID (e.g. NA=North America, EU=Europe, AS=Asia)
RegionID="SA" # South America

# Bounds of the cropped domain
LonMin=-88.0
LonMax=-31.0
LatMin=-59.0
LatMax=16.0

# Directory containing global met fields
DownloadGlobalMet=true
RemoveGlobalMet=false
GlobalDir="./GEOS_0.25x0.3125/GEOS_FP"

# Directory for storing cropped met fields
RegionalDir="./GEOS_0.25x0.3125_${RegionID}/GEOS_FP"

# Download and crop constant meteorology fields? This file is required
# by GEOS-Che, but only needs to be processed once per region (i.e. set
# to false upon subsequent executions of this script)
ProcessConstantFields=true

# Dates to process
Year=2022
StartMonth=1
EndMonth=1

##############################################################################
# Routine crop_met begins here!
##############################################################################

# Create directories if they don't already exist
if [[ ! -d $GlobalDir/$Year ]]; then
    mkdir -p -v $GlobalDir
fi
if [[ ! -d $RegionalDir/$Year ]]; then
    mkdir -p -v $RegionalDir/$Year
fi

# Download and crop constant met fields
# For 0.25x0.3125 GEOS-FP meteorology, the timestamp of these fields is 20110101
if "$ProcessConstantFields"; then
    wget -r -nH --cut-dirs=3 -e robots=off -q -P "${GlobalDir}" "http://geoschemdata.wustl.edu/ExtData/GEOS_0.25x0.3125/GEOS_FP/2011/01/GEOSFP.20110101.CN.025x03125.nc"
    mkdir -p -v $RegionalDir/2011/01
    cdo --no_history sellonlatbox,$LonMin,$LonMax,$LatMin,$LatMax ${GlobalDir}/2011/01/GEOSFP.20110101.CN.025x03125.nc ${RegionalDir}/2011/01/GEOSFP.20110101.CN.025x03125.${RegionID}.nc
fi

# Loop over months
for ((Month = $StartMonth; Month <= $EndMonth; Month++)); do

    # Set mm string
    if (($Month < 10)); then
        mm="0${Month}"
    else
        mm="${Month}"
    fi

    printf "Processing month $mm\n"

    # Download global meteorology fields for this month
    if "$DownloadGlobalMet"; then
        printf "Downloading global meteorology files for ${Year}/${mm}\n"
        wget -r -nH --cut-dirs=3 -e robots=off -nv -o "download_met.log" -P "${GlobalDir}" -A "*.nc" "http://geoschemdata.wustl.edu/ExtData/GEOS_0.25x0.3125/GEOS_FP/${Year}/${mm}/"
        printf "Done downloading global meteorology files\n"
    fi
    # Path to global files
    InPath="${GlobalDir}/${Year}/${mm}"

    # Create directory if it doesn't exist
    if [[ ! -d $RegionalDir/$Year/$mm ]]; then
        mkdir -p -v $RegionalDir/$Year/$mm
    fi

    for file in $InPath/*.nc*; do

        # Output filename
        fOut=${file##*/}                  # Remove input path
        fOut=${fOut%.*}                   # Remove file suffix
        fOut=${fOut}.${RegionID}.nc       # Append region ID
        fOut=$RegionalDir/$Year/$mm/$fOut # Prepend output file path

        echo $fOut

        # Crop global files
        cdo --no_history sellonlatbox,$LonMin,$LonMax,$LatMin,$LatMax $file $fOut

        # Chunk and compress (deflate level=5) files
        nc_chunk.pl $fOut 5

    done

    # Download global met fields for this month
    if "$RemoveGlobalMet"; then
        rm -rf $InPath
    fi

done
