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
    printf "\n Setup     : $( [[ ! -z $setup_end ]] && echo $(( $setup_end - $setup_start )) || echo 0 )"
    printf "\n Spinup     : $( [[ ! -z $spinup_end ]] && echo $(( $spinup_end - $spinup_start )) || echo 0 )"
    printf "\n Jacobian     : $( [[ ! -z $jacobian_end ]] && echo $(( $jacobian_end - $jacobian_start )) || echo 0 )"
    printf "\n Inversion     : $( [[ ! -z $inversion_end ]] && echo $(( $inversion_end - $inversion_start )) || echo 0 )"
    printf "\n Posterior     : $( [[ ! -z $posterior_end ]] && echo $(( $posterior_end - $posterior_start )) || echo 0 )\n\n"
}

# Description: Print error message for if the IMI fails
#   Copy output file to output directory if it exists
# Usage:
#   imi_failed
imi_failed() {
    file=`basename "$0"`
    printf "\nFATAL ERROR on line number ${1} of ${file}: IMI exiting."
    if [ -d "${OutputPath}/${RunName}" ]; then
        cp "${InversionPath}/imi_output.log" "${OutputPath}/${RunName}/imi_output.log"
    fi
    exit 1
}

# Description: Print max value of given variable in netCDF file
#   Returns int if only trailing zeros, float otherwise
#   if (optional) 3rd argument is true, add 4 to the max value
#   this is useful for adding in elements to account for 
#   optimization of BCs
# Usage:
#   ncmax <variable> <netCDF file path> <optimize BCs>
ncmax() {
    python -c "import sys; import xarray;\
    bc_offset = 4 if len(sys.argv) > 3 and sys.argv[3] == 'true' else 0;\
    print('%g' % (xarray.open_dataset(sys.argv[2])[sys.argv[1]].max()+bc_offset))" $1 $2 $3
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
    print('%g, %g' % (float(sys.argv[3]) - diff, float(sys.argv[4]) + diff))" $1 $2 $3 $4
}