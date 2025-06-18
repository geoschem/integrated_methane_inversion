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
    printf "\nRuntime statistics (s):\n"
    printf " Setup      : %s\n" $(( ${setup_end:-0} - ${setup_start:-0} ))
    printf " Spinup     : %s\n" $(( ${spinup_end:-0} - ${spinup_start:-0} ))
    printf " Jacobian   : %s\n" $(( ${jacobian_end:-0} - ${jacobian_start:-0} ))
    printf " Inversion  : %s\n" $(( ${inversion_end:-0} - ${inversion_start:-0} ))
    printf " Posterior  : %s\n\n" $(( ${posterior_end:-0} - ${posterior_start:-0} ))
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

# Description: create CS grids
generate_grid() {
    python - "$1" "$2" <<EOF
import sys
import gcpy
import xarray as xr

CS_RES = int(sys.argv[1])
output_file = sys.argv[2]

grid = gcpy.gen_grid(CS_RES)
area = gcpy.grid_area(cs_grid=grid)
grid = grid.rename({
    'Ydim_b': 'YCdim',
    'Xdim_b': 'XCdim',
    'lat': 'lats',
    'lon': 'lons',
    'lat_b': 'corner_lats',
    'lon_b': 'corner_lons'
})

grid_out = xr.Dataset({})
# Add area as a DataArray
grid_out['area'] = xr.DataArray(
    area,
    dims=['nf', 'Ydim', 'Xdim'],
    coords={
        'lats': (['nf', 'Ydim', 'Xdim'], grid['lats'].values),
        'lons': (['nf', 'Ydim', 'Xdim'], grid['lons'].values)
    },
    attrs={
        'long_name': 'surface area of each grid box',
        'units': 'm2'
    }
)
grid_out['lats'].attrs = dict(long_name='Latitude', units='degree_north')
grid_out['lons'].attrs = dict(long_name='Longitude', units='degree_east')
grid_out['corner_lons'] = xr.DataArray(
    grid['corner_lons'].values,
    dims=['nf', 'YCdim', 'XCdim'],
    attrs={
        'long_name': 'Longitude',
        'units': 'degree_east'
    }
)
grid_out['corner_lats'] = xr.DataArray(
    grid['corner_lats'].values,
    dims=['nf', 'YCdim', 'XCdim'],
    attrs={
        'long_name': 'Latitude',
        'units': 'degree_north'
    }
)
grid_out.to_netcdf(output_file)
EOF
}

# Description: Generate gridded global OH scaling factor
# Usage:
#   gridded_optimized_OH <NH_scale> <SH_scale> <Hemis_mask_fpath> <Output_fpath> <OptimizeNorth> <OptimizeSouth>
gridded_optimized_OH() {
    python - "$1" "$2" "$3" "$4" "$5" "$6" <<EOF
import sys
import xarray as xr
import numpy as np

N_HEMIS = float(sys.argv[1])
S_HEMIS = float(sys.argv[2])
maskfpath = sys.argv[3]
outfpath = sys.argv[4]
OptimizeNorth = sys.argv[5]
OptimizeSouth = sys.argv[6]

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

refyear = 2000
scaleds['oh_scale'].attrs = dict(long_name='Scaling factor for OH', units='1')
scaleds['time'].attrs = dict(units='days since {}-01-01 00:00:00'.format(refyear),
                            delta_t='0000-01-00 00:00:00', axis='T', standard_name='Time',
                            long_name='Time', calendar='standard')
if outfpath is not None:
    print(f"Saving file {outfpath}")
    scaleds.to_netcdf(
        outfpath,
        encoding={v: {"zlib": True, "complevel": 1} for v in scaleds.data_vars},
    )
EOF
}

# Calculate RunDuration (YYYYMMDD) between two dates (YYYYMMDD)
get_run_duration() {
  local start_date=$1
  local end_date=$2

  # Convert to seconds since epoch
  local start_sec=$(date -d "$start_date" +%s)
  local end_sec=$(date -d "$end_date" +%s)

  # Calculate total days difference
  local diff_days=$(( (end_sec - start_sec) / 86400 ))

  # Calculate years, months, days assuming:
  # 1 year = 360 days (12 * 30), 1 month = 30 days
  local years=$(( diff_days / 360 ))
  local remainder=$(( diff_days % 360 ))
  local months=$(( remainder / 30 ))
  local days=$(( remainder % 30 ))

  printf "%04d%02d%02d\n" $years $months $days
}

# Add one day to a RunDuration string YYYYMMDD (duration, not date)
add_one_day_to_duration() {
  local run_duration=$1

  local years=${run_duration:0:4}
  local months=${run_duration:4:2}
  local days=${run_duration:6:2}

  years=$((10#$years))
  months=$((10#$months))
  days=$((10#$days + 1))

  # Handle day rollover (30 days = 1 month)
  if [ $days -ge 30 ]; then
    days=0
    months=$((months + 1))
  fi

  # Handle month rollover (12 months = 1 year)
  if [ $months -ge 12 ]; then
    months=0
    years=$((years + 1))
  fi

  printf "%04d%02d%02d\n" $years $months $days
}