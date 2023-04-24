#!/bin/bash

# Read in the config file
source ../utilities/parse_yaml.sh
eval $(parse_yaml config_write_BCs.yml)

# Make directories if they don't exists
mkdir -p ${workdir}/step1 ${workdir}/step2 ${workdir}/step3 ${workdir}/smoothed-boundary-conditions

if "$RunGEOSChem"; then
    # Load modules
    source ${imidir}/envs/Harvard-Cannon/gcc.gfortran10.2_cannon.env

    # Remove directories if they exist
    rm -rf ${imidir}/src/write_BCs/GCClassic
    rm -rf ${workdir}/runGCC1402

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
    cd ${workdir}/runGCC1402
    sed -i -e "s|start_date: \[[ ]*.*[ ]*\]|start_date: \[${startdate:0:8}, ${startdate:9:15}\]|g"\
           -e "s|end_date: \[[ ]*.*[ ]*\]|end_date: \[${enddate:0:8}, ${enddate:9:15}\]|g" geoschem_config.yml

    # Modify HEMCO_Config.rc.gmao_metfields to allow us to run right up until the end of our current MET field timeseries
    sed -i -e "s|EFY|CY|g" HEMCO_Config.rc.gmao_metfields

    # Compile GEOS-Chem
    cd ${workdir}/runGCC1402/build
    cmake ${imidir}/src/write_BCs/GCClassic >> build_geoschem.log 2>&1
    cmake . -DRUNDIR=..  >> build_geoschem.log 2>&1
    make -j install >> build_geoschem.log 2>&1

    # Run GEOS-Chem
    cd ${workdir}/runGCC1402/
    cp /n/holylfs05/LABS/jacob_lab/imi/ch4/boundary-conditions/restarts/GEOSChem.Restart.${startdate:0:8}_0000z.nc4 ./Restarts/
    sbatch -W ${imidir}/src/write_BCs/GC_config_files/geoschem.run; wait;
fi

if "$WriteBCs"; then
    
    # Run python scripts
    cd ${imidir}/src/write_BCs
    sbatch -W -p seas_compute -t 2-00:00 --mem 64000 -c 48 --wrap "source ~/.bashrc; conda activate $CondaEnv; python write_tropomi_GC_daily_avgs.py"; wait;
    sbatch -W -p seas_compute -t 2-00:00 --mem 64000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python calculate_bias.py"; wait;
    sbatch -W -p seas_compute -t 2-00:00 --mem 64000 --wrap "source ~/.bashrc; conda activate $CondaEnv; python write_boundary.py"; wait;

    # Replace the days we don't have TROPOMI data with initial GC outputs
    cp ${workdir}/runGCC1402/OutputDir/GEOSChem.BoundaryConditions.201804{01..29}_0000z.nc4 ${workdir}/smoothed-boundary-conditions

fi