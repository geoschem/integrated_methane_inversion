#!/bin/bash

# Common shell function for the IMI
# Functions available in this file include:
#   - print_stats
#   - imi_failed
#   - ncmax
#   - ncmin

# Description:
#   Print runtime stats based on existing variables
# Usage:
#   print_stats
print_stats() {
    printf "\nRuntime statistics (s):"
    printf "\n Setup     : $([[ ! -z $setup_end ]] && echo $(($setup_end - $setup_start)) || echo 0)"
    printf "\n Spinup     : $([[ ! -z $spinup_end ]] && echo $(($spinup_end - $spinup_start)) || echo 0)"
    printf "\n Jacobian     : $([[ ! -z $jacobian_end ]] && echo $(($jacobian_end - $jacobian_start)) || echo 0)"
    printf "\n Inversion     : $([[ ! -z $inversion_end ]] && echo $(($inversion_end - $inversion_start)) || echo 0)"
    printf "\n Posterior     : $([[ ! -z $posterior_end ]] && echo $(($posterior_end - $posterior_start)) || echo 0)\n\n"
}

# Description: Print error message for if the IMI fails
#   Copy output file to output directory if it exists
# Usage:
#   imi_failed
imi_failed() {
    file=$(basename "$0")
    printf "\nFATAL ERROR on line number ${1} of ${file}: IMI exiting."
    if [ -d "${OutputPath}/${RunName}" ]; then
        cp "${InversionPath}/imi_output.log" "${OutputPath}/${RunName}/imi_output.log"
    fi
    exit 1
}

# Description: Print max value of given variable in netCDF file
#   Returns int if only trailing zeros, float otherwise
# Usage:
#   ncmax <variable> <netCDF file path>
ncmax() {
    python -c "import sys; import xarray;\
    print('%g' % xarray.open_dataset(sys.argv[2])[sys.argv[1]].max())" $1 $2
}

# Description: Print min value of given variable in netCDF file
#   Returns int if only trailing zeros, float otherwise
# Usage:
#   ncmax <variable> <netCDF file path>
ncmin() {
    python -c "import sys; import xarray; \
    print('%g' % xarray.open_dataset(sys.argv[2])[sys.argv[1]].min())" $1 $2
}

# Description: Add/Subtract half the spacing between coordinate values
#   of the given NetCDF variable to the min and max values.
#   This is useful for adjusting lat/lon bounds because GEOS-Chem
#   uses the grid cell edges, not centers, for the lat/lon bounds
# Usage:
#   calculate_geoschem_domain <variable> <netCDF file path> <min> <max>
calculate_geoschem_domain() {
    python -c "import sys; import xarray; import numpy as np; \
    sv = xarray.open_dataset(sys.argv[2])[sys.argv[1]]; \
    diff = np.unique(sv.diff(sys.argv[1])).item()/2; \
    print('%f, %f' % (float(sys.argv[3]) - diff, float(sys.argv[4]) + diff))" $1 $2 $3 $4
}

# Description: Generate gridded global OH scaling factor
# Usage:
#   gridded_optimized_OH <NH_scale> <SH_scale> <Hemis_mask_fpath> <Output_fpath> <OptimizeNorth> <OptimizeSouth> [STRETCH_GRID STRETCH_FACTOR TARGET_LAT TARGET_LON]
gridded_optimized_OH() {
    if [ "$#" -lt 6 ]; then
        echo "Usage: gridded_optimized_OH <NH_scale> <SH_scale> <Hemis_mask_fpath> <Output_fpath> <OptimizeNorth> <OptimizeSouth> [STRETCH_GRID STRETCH_FACTOR TARGET_LAT TARGET_LON]"
        return 1
    fi

    local NH_scale="$1"
    local SH_scale="$2"
    local Hemis_mask_fpath="$3"
    local Output_fpath="$4"
    local OptimizeNorth="$5"
    local OptimizeSouth="$6"
    local STRETCH_GRID="${7:-False}"
    local STRETCH_FACTOR="${8:-1.0}"
    local TARGET_LAT="${9:--90.0}"
    local TARGET_LON="${10:-170.0}"

    python - "$NH_scale" "$SH_scale" "$Hemis_mask_fpath" "$Output_fpath" "$OptimizeNorth" "$OptimizeSouth" "$STRETCH_GRID" "$STRETCH_FACTOR" "$TARGET_LAT" "$TARGET_LON" <<EOF
import sys
import xarray as xr
import numpy as np
N_HEMIS = float(sys.argv[1])
S_HEMIS = float(sys.argv[2])
maskfpath = sys.argv[3]
outfpath = sys.argv[4]
OptimizeNorth = sys.argv[5]
OptimizeSouth = sys.argv[6]
STRETCH_GRID = sys.argv[7].lower() == "true"
STRETCH_FACTOR = float(sys.argv[8])
TARGET_LAT = float(sys.argv[9])
TARGET_LON = float(sys.argv[10])
maskds = xr.open_dataset(maskfpath)
hemis = maskds['Hemisphere']
# Define masks using approximate equality (tolerance 1e-6)
n_mask = abs(hemis - 1.0) < 1e-6
s_mask = abs(hemis - 2.0) < 1e-6
# Initialize oh_scale as all ones with the same dims and coords as hemis
oh_scale = xr.ones_like(hemis)
if OptimizeNorth.lower() == 'true':
    oh_scale = oh_scale.where(~n_mask, other=N_HEMIS)
if OptimizeSouth.lower() == 'true':
    oh_scale = oh_scale.where(~s_mask, other=S_HEMIS)
# Drop Hemisphere and add oh_scale with attributes
scaleds = maskds.drop_vars('Hemisphere').assign(
    oh_scale=oh_scale
)
scaleds['oh_scale'].attrs = dict(long_name='Scaling factor for OH', units='1')
if STRETCH_GRID:
    scaleds.attrs['STRETCH_FACTOR'] = np.float32(STRETCH_FACTOR)
    scaleds.attrs['TARGET_LAT'] = np.float32(TARGET_LAT)
    scaleds.attrs['TARGET_LON'] = np.float32(TARGET_LON)
if outfpath is not None:
    print(f"Saving file {outfpath}")
    scaleds.to_netcdf(
        outfpath,
        encoding={v: {"zlib": True, "complevel": 1} for v in scaleds.data_vars},
    )
EOF
}