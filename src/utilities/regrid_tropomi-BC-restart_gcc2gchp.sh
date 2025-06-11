#!/bin/bash

# Regrid a GCClassic restart file (from TROPOMI boundary conditions) to a GCHP cube-sphere restart file
#
# Calling sequence:
# regrid_tropomi-BC-restart_gcc2gchp.sh {TROPOMI-BC} {FILE-PREFIX} {CS-RESOLUTION}
#
# Example:
# tropomi_bc="/n/holylfs05/LABS/jacob_lab/imi/ch4/tropomi-boundary-conditions/v2024-06/GEOSChem.BoundaryConditions.20180501_0000z.nc4"
# template="GEOSChem.Restart.20190101_0000z.c48.nc4"
# ./regrid_tropomi-BC-restart_gcc2gchp.sh $tropomi_bc $template GEOSChem.Restart.20180501_0000z 48

if [[ $# == 4 ]] ; then
    tropomi_bc=$1
    template=$2
    prefix=$3
    cs_res=$4
else
    echo "Usage: ./regrid_rst_gcc2gchp.sh {TROPOMI-BC} {Template} {FILE-PREFIX} {CS-RESOLUTION}"
    exit 1
fi

#mamba activate gcpy_env
cp $tropomi_bc ${prefix}.nc4
cp $template temp.nc4
ncrename -v SpeciesBC_CH4,SpeciesRst_CH4 ${prefix}.nc4
python -m gcpy.file_regrid --filein ${prefix}.nc4 \
       --dim_format_in classic \
       --fileout ${prefix}.c${cs_res}.nc4 \
       --cs_res_out ${cs_res} \
       --dim_format_out checkpoint

# Remove regridding files
rm -rf conservative_*.nc ${prefix}.nc4
# Add SPC_CH4 field from regridded file into template, and rename final output
ncks -A -v SPC_CH4 ${prefix}.c${cs_res}.nc4 temp.nc4
mv temp.nc4 ${prefix}.c${cs_res}.nc4
#cdo replace $template ${prefix}.c${cs_res}.nc4 ${prefix}.c${cs_res}.nc4

#mamba deactivate
exit 0