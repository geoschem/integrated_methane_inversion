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

# Description: Generate grid prefix for CS or SGCS grid files
# Usage:
# get_GridSpec_prefix <CS_RES> <STRETCH_GRID> [STRETCH_FACTOR TARGET_LAT TARGET_LON]
get_GridSpec_prefix() {
    local CS_RES="$1"
    local STRETCH_GRID="$2"
    local STRETCH_FACTOR="${3:-1.0}"
    local TARGET_LAT="${4:--90.0}"
    local TARGET_LON="${5:-170.0}"

    if [ "$#" -lt 2 ]; then
        echo "Usage: get_GridSpec_prefix <CS_RES> <STRETCH_GRID> [STRETCH_FACTOR TARGET_LAT TARGET_LON]" >&2
        return 1
    fi

    if "$STRETCH_GRID"; then
        # Format stretch factor: e.g., 1.50 -> 1d50
        local sf_formatted
        sf_formatted=$(printf "%.2f" "$STRETCH_FACTOR" | sed 's/\./d/')

        # Compute geohash using Python, safely quoting lat/lon
        local target_geohash
        target_geohash=$(python3 -c "import pygeohash as pgh; print(pgh.encode(float('$TARGET_LAT'), float('$TARGET_LON')))" )

        echo "c${CS_RES}_s${sf_formatted}_t${target_geohash}"
    else
        echo "c${CS_RES}"
    fi
}

# Description: Create CS or SGCS grid and save to netCDF
# Usage: generate_grid_from_GridSpec <CS_RES> <OUTPUT_FILE> <STRETCH_GRID> [STRETCH_FACTOR TARGET_LAT TARGET_LON]
generate_grid_from_GridSpec() {
    local CS_RES="$1"
    local OUTPUT_FILE="$2"
    local STRETCH_GRID="$3"
    local STRETCH_FACTOR="${4:-1.0}"
    local TARGET_LAT="${5:--90.0}"
    local TARGET_LON="${6:-170.0}"

    if [ "$#" -lt 3 ]; then
        echo "Usage: generate_grid_from_GridSpec <CS_RES> <OUTPUT_FILE> <STRETCH_GRID> [STRETCH_FACTOR TARGET_LAT TARGET_LON]"
        return 1
    fi

    local prefix
    prefix=$(get_GridSpec_prefix "$CS_RES" "$STRETCH_GRID" "$STRETCH_FACTOR" "$TARGET_LAT" "$TARGET_LON")

    python - "$CS_RES" "$OUTPUT_FILE" "$STRETCH_GRID" "$STRETCH_FACTOR" "$TARGET_LAT" "$TARGET_LON" "$prefix" <<EOF
import sys
import os
import numpy as np
import xarray as xr
import pygeohash as pgh
CS_RES = int(sys.argv[1])
OUTPUT_FILE = sys.argv[2]
STRETCH_GRID = sys.argv[3].lower() == "true"
STRETCH_FACTOR = float(sys.argv[4])
TARGET_LAT = float(sys.argv[5])
TARGET_LON = float(sys.argv[6])
prefix = sys.argv[7]
lon   = np.zeros((6, CS_RES, CS_RES))
lat   = np.zeros((6, CS_RES, CS_RES))
lon_b = np.zeros((6, CS_RES + 1, CS_RES + 1))
lat_b = np.zeros((6, CS_RES + 1, CS_RES + 1))
area  = np.zeros((6, CS_RES, CS_RES))
gridspec_file = f'{prefix}_gridspec.nc'
if not os.path.isfile(gridspec_file):
    print(f'ERROR: {gridspec_file} does not exist. Use GridSpec to generate it first!', file=sys.stderr)
    sys.exit(1)
for i in range(6):
    tile_idx = i + 1
    filepath = f"{prefix}.tile{tile_idx}.nc"
    if not os.path.isfile(filepath):
        print(f"ERROR: Tile file {filepath} does not exist!", file=sys.stderr)
        sys.exit(1)
    ds = xr.open_dataset(filepath)
    super_lats = ds['lats'].values   # shape (2*CS_RES+1, 2*CS_RES+1)
    super_lons = ds['lons'].values
    tile_area = ds['area'].values    # shape (CS_RES, CS_RES)
    lat_b[i] = super_lats[::2, ::2]
    lon_b[i] = super_lons[::2, ::2]
    lat[i]   = super_lats[1::2, 1::2]
    lon[i]   = super_lons[1::2, 1::2]
    area[i]  = tile_area
# Adjust longitudes to [0, 360)
lon[lon < 0] += 360
lon_b[lon_b < 0] += 360
# Create xarray DataArrays
area_da = xr.DataArray(
    area, 
    dims=['nf', 'Ydim', 'Xdim'], 
    coords={'lats': (['nf', 'Ydim', 'Xdim'], lat),
            'lons': (['nf', 'Ydim', 'Xdim'], lon)},
    attrs=dict(units='m2', long_name='Surface area of each grid box')
)
corner_lons_da = xr.DataArray(
    lon_b, dims=['nf', 'YCdim', 'XCdim'], 
    attrs=dict(units='degrees_east', long_name='Longitude')
)
corner_lats_da = xr.DataArray(
    lat_b, dims=['nf', 'YCdim', 'XCdim'], 
    attrs=dict(units='degrees_north', long_name='Latitude')
)
# Combine into dataset
data = xr.Dataset({
    'area': area_da,
    'corner_lons': corner_lons_da,
    'corner_lats': corner_lats_da
})
# Assign coordinate metadata
data['lats'].attrs = dict(units='degrees_north', long_name='Latitude')
data['lons'].attrs = dict(units='degrees_east', long_name='Longitude')
# Add global attributes if stretched grid
if STRETCH_GRID:
    data.attrs['STRETCH_FACTOR'] = np.float32(STRETCH_FACTOR)
    data.attrs['TARGET_LAT'] = np.float32(TARGET_LAT)
    data.attrs['TARGET_LON'] = np.float32(TARGET_LON)
# Save to netCDF
data.to_netcdf(OUTPUT_FILE)
print(f"Combined c{CS_RES} cubed-sphere grid data saved to {OUTPUT_FILE}")
EOF
}