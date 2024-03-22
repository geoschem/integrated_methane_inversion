#!/bin/bash

# Functions available in this file include:
#   - run_prior

# Description: Run a HEMCO standalone simulation to generate prior emissions
#
# Usage:
#   run_prior
run_prior() {
    prior_start=$(date +%s)
    if [[ -d ${RunPrior} ]]; then
	printf "\nERROR: ${PriorDir} already exists. Please remove or set 'DoPriorEmis: false' in config.yml.\n"
	exit 9999
    fi

    printf "\n=== GENERATING PRIOR EMISSIONS ===\n"

    cd ${GCClassicPath}/src/HEMCO/run

    # Define the run directory name
    PriorName="${RunName}_Prior"

    # Commands to feed to createRunDir.sh
    # HEMCO standalone run directories are created for the global domain by default
    if [[ "$Met" == "MERRA2" || "$Met" == "MERRA-2" || "$Met" == "merra2" ]]; then
	metNum="1"
    elif [[ "$Met" == "GEOSFP" || "$Met" == "GEOS-FP" || "$Met" == "geosfp" ]]; then
	metNum="2"
    else
	printf "\nERROR: Meteorology field ${Met} is not supported by the IMI. "
	printf "\n Options are GEOSFP or MERRA2.\n"
	exit 1
    fi	
    if [ "$Res" = "4.0x5.0" ]; then
	resnum="1"
    elif [ "$Res" == "2.0x2.5" ]; then
	resnum="2"
    elif [ "$Res" == "0.5x0.625" ]; then
	resnum="3"
    elif [ "$Res" == "0.25x0.3125" ]; then
	resnum="4"
    else
	printf "\nERROR: Grid resolution ${Res} is not supported by the IMI. "
	printf "\n Options are 0.25x0.3125, 0.5x0.625, 2.0x2.5, or 4.0x5.0.\n"
	exit 1
    fi
    HEMCOconfig=${RunTemplate}/HEMCO_Config.rc
    
    cmd="${metNum}\n${resnum}\n${HEMCOconfig}\n${RunDirs}\n${priorDir}\nn\n"

    # Create HEMCO standalone directory
    printf ${cmd} | ./createRunDir.sh >> createHemcoDir.log 2>&1
    rm -f createHemcoDir.log
    printf "\nCreated ${RunPrior}\n" 

    cd ${RunPrior}

    # Modify HEMCO files based on settings in config.yml
    sed -i -e "s:2019-07-01:${StartDate:0:4}-${StartDate:4:2}-${StartDate:6:2}:g" \
           -e "s:2019-08-01 00:${StartDate:0:4}-${StartDate:4:2}-${StartDate:6:2} 01:g" HEMCO_sa_Time.rc

    sed -i -e "s:_NA::g" -e "s:.NA.:.:g" HEMCO_Config.rc.gmao_metfields

    sed -i -e "/DiagnFreq:           00000100 000000/d" \
	   -e "/Negative values:     0/d" \
	   -e "s/Verbose:             false/Verbose:             true/g" HEMCO_sa_Config.rc
    sed -i -e "/#SBATCH -c 8/d" runHEMCO.sh
    sed -i -e "/#SBATCH -t 0-12:00/d" runHEMCO.sh
    sed -i -e "/#SBATCH -p huce_intel/d" runHEMCO.sh
    sed -i -e "/#SBATCH --mem=15000/d" runHEMCO.sh
    mv runHEMCO.sh ${PriorName}.run

    # Compile HEMCO and store executable in template run directory
    printf "\nCompiling HEMCO...\n"
    cd build
    cmake ${InversionPath}/GCClassic/src/HEMCO >> build_hemco.log 2>&1
    cmake . -DRUNDIR=..  >> build_hemco.log 2>&1 
    make -j install >> build_hemco.log 2>&1
    cd ..
    if [[ -f hemco_standalone ]]; then
	mkdir HEMCO_build_info
	mv build/CMakeCache.txt HEMCO_build_info
        rm -rf build
    else
        printf "\nHEMCO build failed! \n\nSee ${RunTemplate}/build/build_hemco.log for details\n"
        exit 999
    fi
    printf "\nDone compiling HEMCO \n\nSee ${RunDirs}/HEMCO_build_info for details\n\n"

    printf "\nSubmitting prior emissions simulation\n\n"

    # Submit job to job scheduler
    sbatch --mem $PriorMemory \
    -c $SimulationCPUs \
    -t $RequestedTime \
    -p $SchedulerPartition \
    -W ${PriorName}.run; wait;

    # check if exited with non-zero exit code
    [ ! -f ".error_status_file.txt" ] || imi_failed $LINENO
    
    printf "\nDone prior emissions simulation\n\n"
    
    printf "\n=== DONE GENERATING PRIOR EMISSIONS ===\n"
    prior_end=$(date +%s)
}
