#!/bin/bash
#SBATCH --job-name=boundary_conditions
#SBATCH --mem=4000
#SBATCH --time=07-00:00
#SBATCH --output=debug.log

cwd="$(pwd)"

# Read in the config file and source the environment file
condaEnv=$(grep -Po 'condaEnv:\s*\K.*' config_boundary_conditions.yml)
condaFile=$(grep -Po 'condaFile:\s*\K.*' config_boundary_conditions.yml)
condaFile=$(eval echo "$condaFile")
source ${condaFile}
conda activate ${condaEnv}
eval $(python ../../src/utilities/parse_yaml.py config_boundary_conditions.yml)
source ${geosChemEnv}
echo "Environment file  --> ${geosChemEnv}" >> "${cwd}/boundary_conditions.log"

# Information needed for GEOS-Chem simulations
export GC_USER_REGISTERED=true
export GC_DATA_ROOT="${geosChemDataPath}"

# As long as it doesn't exist, make the working directory and go to it
if [[ -d "${workDir}" ]]; then
    echo "ERROR             --> Directory ${workDir} exists." >> "${cwd}/boundary_conditions.log"
    exit 1
fi
mkdir -p "${workDir}"
echo "Working directory --> ${workDir}" >> "${cwd}/boundary_conditions.log"
mkdir -p "${workDir}/tropomi-boundary-conditions"
mkdir -p "${workDir}/blended-boundary-conditions"
cd "${workDir}"

# Get GCClassic v14.4.1 and create the run directory
git clone https://github.com/geoschem/GCClassic.git
cd GCClassic
git checkout 14.4.1
git submodule update --init --recursive
cd run
runDir="gc_run"
c="9\n2\n2\n2\n${workDir}\n${runDir}\nn\n" # CH4, GEOS-FP, 2.0 x 2.5, 47L
printf ${c} | ./createRunDir.sh
cd "${workDir}/${runDir}/build"
cmake ../CodeDir -DRUNDIR=..
make -j
make install
cd "${workDir}/${runDir}"

# Modify HISTORY.rc (hourly instantaneous CH4/pressure/air and 3-hourly BCs)
sed -i -e "s|'CH4',|#'CH4',|g" \
    -e "s|'Metrics',|#'Metrics',|g" \
    -e "s|#'LevelEdgeDiags',|'LevelEdgeDiags',|g" \
    -e "s|Restart.frequency:          'End',|Restart.frequency:          '00000001 000000',|g" \
    -e "s|Restart.duration:           'End',|Restart.duration:           '00000001 000000',|g" \
    -e "s|00000100 000000|00000000 010000|g" \
    -e "s|time-averaged|instantaneous|g" \
    -e "s|Met_CMFMC|Met_PEDGE|g" \
    -e "s|#'BoundaryConditions',|'BoundaryConditions',|g" \
    -e "s|'Met_AD                        ',|'Met_AIRVOL                    ',|g" HISTORY.rc

# Remove unnecessary StateMet and then LevelEdge variables
sed -i '265,340d' HISTORY.rc
sed -i '195,200d' HISTORY.rc

# Modify HEMCO_Config.rc so that GEOS-Chem can run into 2024
sed -i '/GFED4/s/ RF/ C/g' HEMCO_Config.rc

# Modify geoschem_config.yml
# - run GC earlier than you want BCs to accomodate a 15/30 day average going back in time
# - e.g., the BCs for 15 May 2023 require 1 May 2023-15 May 2023 data (and sometimes 16 April 2023-15 May 2023)
# - run one day earlier than that because GEOS-Chem won't write SpeciesConc/LevelEdgeDiag for t = 0
if [[ ${startDate} -ge "20180501" ]]; then
    gcStartDate=$(date -d "$startDate -30 days" +%Y%m%d)
else
    gcStartDate="20180401"
fi
# - run GC to 00:00:00 the day after you want BCs
# - e.g., you want BCs for 31 Aug 2023 -> run GC to 1 Sep 2023 00:00:00
gcEndDate=$(date -d "$endDate +1 days" +%Y%m%d)
sed -i -e "s|start_date: \[20190101, 000000\]|start_date: [${gcStartDate}, 000000]|g" \
    -e "s|end_date: \[20190201, 000000\]|end_date: [${gcEndDate}, 000000]|g" geoschem_config.yml
echo "GC run times      --> ${gcStartDate} 00:00:00 until ${gcEndDate} 00:00:00" >> "${cwd}/boundary_conditions.log"

# Prepare the restart file
rm Restarts/GEOSChem.Restart.20190101_0000z.nc4
if [[ ${gcStartDate} -eq "20180401" ]]; then #  use the restart file provided with the IMI
    cp "${cwd}/GEOSChem.Restart.20180401_0000z.nc4" Restarts/
else 
    if [[ -e ${restartFilePath} ]]; then # use your own restart file
        restartFileDate=$(echo "${restartFilePath}" | grep -oP '\d{8}')
        if [[ ${restartFileDate} -ne ${gcStartDate} ]]; then
            echo "ERROR             --> your restart file date (${restartFileDate}) doesn't match the simulation start date (${gcStartDate})" >> "${cwd}/boundary_conditions.log"
            exit 1
        else
            cp "${restartFilePath}" Restarts/
        fi
    else
        echo "ERROR             --> the restart file does not exist!" >> "${cwd}/boundary_conditions.log"
        exit 1
    fi
fi

# Remove debug log file if not debug (everything written via >> debug.log 2>&1)
if ! ${debug}; then
    rm "${cwd}/debug.log"
fi

# Modify and submit the run script
cp runScriptSamples/operational_examples/harvard_cannon/geoschem.run .
sed -i -e "s|sapphire,huce_cascade,seas_compute,shared|${partition}|g" \
    -e "s|--mem=15000|--mem=128000|g" \
    -e "s|-t 0-12:00|-t 07-00:00|g"\
    -e "s|-c 8|-c 48|g" geoschem.run
sbatch -W geoschem.run; wait;

# Write the boundary conditions using write_boundary_conditions.py
cd "${cwd}"
sbatch -W -J blended -o boundary_conditions.log --open-mode=append -p ${partition} -t 7-00:00 --mem 96000 -c 40 --wrap "source $condaFile; conda activate $condaEnv; python write_boundary_conditions.py True $blendedDir $gcStartDate $gcEndDate"; wait; # run for Blended TROPOMI+GOSAT
sbatch -W -J tropomi -o boundary_conditions.log --open-mode=append -p ${partition} -t 7-00:00 --mem 96000 -c 40 --wrap "source $condaFile; conda activate $condaEnv; python write_boundary_conditions.py False $tropomiDir $gcStartDate $gcEndDate"; wait; # run for TROPOMI data
echo "" >> "${cwd}/boundary_conditions.log"
echo "Blended TROPOMI+GOSAT boundary conditions --> ${workDir}/blended-boundary-conditions" >> "${cwd}/boundary_conditions.log"
echo "TROPOMI boundary conditions               --> ${workDir}/tropomi-boundary-conditions" >> "${cwd}/boundary_conditions.log"