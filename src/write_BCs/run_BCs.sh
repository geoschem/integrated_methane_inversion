#!/bin/bash

# Read in the config file
source ../utilities/parse_yaml.sh
eval $(parse_yaml config_write_BCs.yml)

# Make directories if they don't exists (and throw an error if they are not empty so we don't accidentally overwrite anything)
function create_and_check_dirs() {
    local dirs=("$@")
    for dir in "${dirs[@]}"; do
        mkdir -p "$dir"
        if [[ $(find "$dir" -mindepth 1 -print -quit) ]]; then
            echo "Error: Directory $dir contains files."
            exit 1
        fi
    done
}

if "$RunGEOSChem"; then
    # Load modules
    source ${GEOSChemEnv}

    # Make sure runGCC1402 is empty so we don't accidentally overwrite it
    dirs=("${workdir}/runGCC1402")
    create_and_check_dirs "${dirs[@]}"

    # Remove GCClassic if it exists
    rm -rf ${workdir}/runGCC1402
    rm -rf ${imidir}/src/write_BCs/GCClassic

    # Download GEOS-Chem 14.0.2
    git clone https://github.com/geoschem/GCClassic.git
    cd GCClassic
    git checkout 14.0.2
    git submodule update --init --recursive
    cd run
    cmd="3\n2\n1\n2\n${workdir}\nrunGCC1402\nn\n" # these defaults will be replaced by the files in GC_config_files anyway
    printf ${cmd} | ./createRunDir.sh >> createRunDir.log 2>&1
    rm -f createRunDir.log

    # Overwrite the default config files
    cd ${imidir}/src/write_BCs/GC_config_files
    cp geoschem_config.yml ${workdir}/runGCC1402/
    cp HEMCO_Config.rc ${workdir}/runGCC1402/
    cp HISTORY.rc ${workdir}/runGCC1402/

    # Modify geoschem_config.yml to match our start and end date specified in the config file
    # e.g., if you want BCs for 20180401-20230331, run GC for 20180401T000000-20230420T000000 (accounts for +/- 15 day temporal averaging in calculate_bias.py)
    cd ${workdir}/runGCC1402
    gc_enddate=$(date -d "$enddate +20 days" +%Y%m%d)
    sed -i -e "s|start_date: \[[ ]*.*[ ]*\]|start_date: \[${startdate}, 000000\]|g"\
           -e "s|end_date: \[[ ]*.*[ ]*\]|end_date: \[${gc_enddate}, 000000\]|g" geoschem_config.yml

    # Compile GEOS-Chem
    cd ${workdir}/runGCC1402/build
    cmake ${imidir}/src/write_BCs/GCClassic >> build_geoschem.log 2>&1
    cmake . -DRUNDIR=..  >> build_geoschem.log 2>&1
    make -j install >> build_geoschem.log 2>&1

    # Run GEOS-Chem
    cd ${workdir}/runGCC1402/
    cp ${RestartDir}/GEOSChem.Restart.${startdate:0:8}_0000z.nc4 ./Restarts/
    sbatch -W ${imidir}/src/write_BCs/GC_config_files/geoschem.run; wait;
fi

if "$WriteBCs"; then

    # Make sure dirs are empty so we don't accidentally overwrite them
    dirs=("${workdir}/step1" "${workdir}/step2" "${workdir}/step3" "${workdir}/smoothed-boundary-conditions")
    create_and_check_dirs "${dirs[@]}"
    
    # Run python scripts
    cd ${imidir}/src/write_BCs
    sbatch -W -p ${Partition} -t 2-00:00 --mem 190000 -c 48 --wrap "source ~/.bashrc; conda activate $CondaEnv; python write_tropomi_GC_daily_avgs.py"; wait;
    sbatch -W -p ${Partition} -t 2-00:00 --mem 64000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python calculate_bias.py"; wait;
    sbatch -W -p ${Partition} -t 2-00:00 --mem 64000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python write_boundary.py"; wait;

    # Replace the days we don't have TROPOMI data with initial GC outputs
    cp ${workdir}/runGCC1402/OutputDir/GEOSChem.BoundaryConditions.201804{01..29}_0000z.nc4 ${workdir}/smoothed-boundary-conditions

fi