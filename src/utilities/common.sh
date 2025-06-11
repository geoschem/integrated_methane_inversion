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

# Description: create CS grids
generate_grid() {
    python - "$1" "$2" <<EOF
import sys
import gcpy

CS_RES = sys.argv[1]
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

# Add area as a DataArray
grid['area'] = xr.DataArray(
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

grid.to_netcdf(output_file)
EOF
}

# Description: Generate gridded global OH scaling factor
# Usage:
#   gridded_optimized_OH <NH_scale> <SH_scale> <Hemis_mask_fpath> <Output_fpath>
gridded_optimized_OH() {
    python - "$1" "$2" "$3" "$4" <<EOF
import sys
import xarray as xr
import numpy as np

N_HEMIS = float(sys.argv[1])
S_HEMIS = float(sys.argv[2])
maskfpath = sys.argv[3]
outfpath = sys.argv[4]

maskds = xr.open_dataset(maskfpath)
hemis = maskds['Hemisphere'].values

# Create masks
n_mask = np.isclose(hemis, 1.0)
s_mask = np.isclose(hemis, 2.0)

# Construct OH scaling field
oh_scale = np.ones_like(hemis)
oh_scale[n_mask] = N_HEMIS
oh_scale[s_mask] = S_HEMIS

# Add as new variable
scaleds = maskds.copy(deep=True)
scaleds['oh_scale'] = xr.DataArray(
    oh_scale, dims=maskds['Hemisphere'].dims,
    attrs=dict(long_name='Scaling factor for OH', units='1')
)

# Save to file
scaleds.to_netcdf(outfpath)
EOF
}