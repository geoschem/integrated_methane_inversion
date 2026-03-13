#!/bin/bash

# Common shell function for the IMI
# Functions available in this file include:
#   - submit_job
#       - submit_slurm_job
#       - submit_pbs_job
#   - print_stats
#   - imi_failed
#   - ncmax
#   - ncmin

# Description: 
#   Submit a job with default ICI settings using either SBATCH or PBS
# Usage:
#   submit_job <SchedulerType> <output_bool> <Memory> <CPUs> <Time> <Partition> <Script> <AdditionalArguments>
submit_job() {
    if [[ $1 = "slurm" || $1 = "tmux" ]]; then
        submit_slurm_job "${@:2}"
    elif [[ $1 = "PBS" ]]; then
        submit_pbs_job "${@:2}"
    else
        echo "Scheduler type $1 not recognized."
    fi

    # If output was saved, concatenate it to imi_output
    if [[ $2 = "true" ]]; then
        cat imi_output.tmp >> ${InversionPath}/imi_output.log
        rm imi_output.tmp
    fi
}

# Description: 
#   Submit a job with default ICI settings using SBATCH
# Usage:
#   submit_slurm_job $JobArguments
submit_slurm_job() {
    if [[ $1 = "true" ]]; then
        sbatch -N 1 \
            --mem $2 \
            -c $3 \
            -t $4 \
            -p $5 \
            -o imi_output.tmp \
            -W ${@:6}; wait;
    else
        sbatch -N 1 \
            --mem $2 \
            -c $3 \
            -t $4 \
            -p $5 \
            -W ${@:6}; wait;
    fi
}

# Description: 
#   Submit a job with default ICI settings using PBS
# Usage:
#   submit_pbs_job $JobArguments
submit_pbs_job() {
    # If save output
    if [[ $1 = "true" ]]; then
        qsub -lselect=1:ncpus=$3:mem=$2:model=ivy \
            -l walltime=$4 -q devel -o imi_output.tmp \
            -Wblock=true -- ${@:6}; wait;
    else
        qsub -lselect=1:ncpus=$3:mem=$2:model=ivy \
            -l walltime=$4 -q devel \
            -Wblock=true -- ${@:6}; wait;
    fi
}

convert_sbatch_to_pbs() {
    DataPaths=($OutputPath $DataPath $DataPathObs $HOME)
    declare -a SitesNeeded=()
    for DP in ${DataPaths[@]}; do
        SitesNeeded_DP=$( find $DP/ -type l -exec realpath {} \; | cut -d/ -f2 | sort -u )
        for NS in ${SitesNeeded_DP[*]}; do
            if ! [[ ${SitesNeeded[@]} =~ $NS ]]; then
                SitesNeeded+=("${NS}+")
            fi
        done
    done
    SitesNeeded=$(IFS=/ ; echo "${SitesNeeded[*]}")
    SitesNeeded="/${SitesNeeded::-1}"

    # Get files containing SBATCH
    current_dir=$(pwd)
    sbatch_files=($(grep -rl "SBATCH" . --exclude-dir={"GCClassic",".git","*utilities*"}))
    echo "Replacing SBATCH with PBS in the following files:"
    for file in ${sbatch_files[@]}; do
        f=${current_dir}${file:1}
        echo "    ${f}"

        # First, insert needed sites at the top of every file
        if ! grep -q "PBS -l site=needed" $file; then
            awk -i inplace 'FNR==NR{ if (/^##SBATCH/) p=NR; next} 1; FNR==p{ print "##PBS -l site=needed='${SitesNeeded}'" }' ${f} ${f}
            awk -i inplace 'FNR==NR{ if (/^#SBATCH/) p=NR; next} 1; FNR==p{ print "#PBS -l site=needed='${SitesNeeded}'" }' ${f} ${f}
        fi

        # Replace SBATCH options
        sed -i -e "s/SBATCH -J /PBS -N /g" \
            -e "s/SBATCH -N /PBS -l nodes=/g" \
            -e "s/SBATCH -c /PBS -l ncpus=/g" \
            -e "s/SBATCH --mem /PBS -l mem=/g" \
            -e "s/SBATCH --mem=/PBS -l mem=/g" \
            -e "s/SBATCH -t /PBS -l walltime=/g" \
            -e "s/SBATCH -n /PBS -l nodes=1:ppn=/g" \
            -e "s/SBATCH --ntasks-per-node/PBS -l nodes=1:ppn/g" \
            -e "s/SBATCH -p /PBS -q /g" \
            -e "s/SBATCH -o /PBS -o /g" \
            -e "s/SBATCH --mail-type=END/PBS -m e/g" ${f}
    done
}

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
    printf "\nFATAL ERROR on line number ${1} of ${2}: IMI exiting."
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