#!/bin/bash
#SBATCH --job-name=boundary_conditions
#SBATCH --mem=4000
#SBATCH --time=07-00:00
#SBATCH --output=debug.log

cwd="$(pwd)"

# Read in the config file and source the environment file
source ../utilities/parse_yaml.sh
eval $(parse_yaml config_boundary_conditions.yml)
source ${GEOSChemEnv}
echo "Environment file  --> ${GEOSChemEnv}" >> "${cwd}/boundary_conditions.log"

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

# Get GCClassic v14.2.1 and create the run directory
git clone https://github.com/geoschem/GCClassic.git
cd GCClassic
git checkout dac5a54 # most recent dev/14.2.1 @ 1 Sep 2023 12:44 PM (update this once 14.2.1 officially released)
git submodule update --init --recursive
cd run
runDir="gc_run"
c="3\n2\n2\n2\n${workDir}\n${runDir}\nn\n" # CH4, GEOS-FP, 2.0 x 2.5, 47L
printf ${c} | ./createRunDir.sh
cd "${workDir}/${runDir}/build"
cmake ../CodeDir -DRUNDIR=..
make -j
make install
cd "${workDir}/${runDir}"

# Modify HISTORY.rc (hourly instantaneous CH4/pressure and 3-hourly BCs)
sed -i -e "s|'CH4',|#'CH4',|g" \
    -e "s|'Metrics',|#'Metrics',|g" \
    -e "s|'StateMet',|#'StateMet',|g" \
    -e "s|#'LevelEdgeDiags',|'LevelEdgeDiags',|g" \
    -e "s|Restart.frequency:          'End',|Restart.frequency:          '00000001 000000',|g" \
    -e "s|Restart.duration:           'End',|Restart.duration:           '00000001 000000',|g" \
    -e "s|SpeciesConc.frequency:      00000100 000000|SpeciesConc.frequency:      00000000 010000|g" \
    -e "s|SpeciesConc.duration:       00000100 000000|SpeciesConc.duration:       00000000 010000|g" \
    -e "s|SpeciesConc.mode:           'time-averaged'|SpeciesConc.mode:           'instantaneous'|g" \
    -e "s|'SpeciesConcMND_?ALL?          ',|#'SpeciesConcMND_?ALL?          ',|g" \
    -e "s|LevelEdgeDiags.frequency:   00000100 000000|LevelEdgeDiags.frequency:   00000000 010000|g" \
    -e "s|LevelEdgeDiags.duration:    00000100 000000|LevelEdgeDiags.duration:    00000000 010000|g" \
    -e "s|LevelEdgeDiags.mode:        'time-averaged'|LevelEdgeDiags.mode:        'instantaneous'|g" \
    -e "s|#'BoundaryConditions',|'BoundaryConditions',|g" HISTORY.rc

# Modify geoschem_config.yml
# - run GC earlier than you want BCs to accomodate a 15 day average going back in time
# - e.g., the BCs for 15 May 2023 require 1 May 2023-15 May 2023 data
# - run one day earlier than that because GEOS-Chem won't write SpeciesConc/LevelEdgeDiag for t = 0
if [[ ${startDate} -ge "20180416" ]]; then
    gcStartDate=$(date -d "$startDate -15 days" +%Y%m%d)
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
            echo "ERROR             --> your restart file date (${restartFileDate}) \
                doesn't match the simulation start date (${gcStartDate})" >> "${cwd}/boundary_conditions.log"
            exit 1
        else
            cp "${restartFilePath}" Restarts/
        fi
    else
        echo "ERROR             --> the restart file does not exist!" >> "${cwd}/boundary_conditions.log"
        exit 1
    fi
fi

# Remove debug log file if not debug (everything written via  >> debug.log 2>&1)
if ! ${debug}; then
    rm "${cwd}/debug.log"
fi

# Modify and submit the run script
cp runScriptSamples/operational_examples/harvard_cannon/geoschem.run .
sed -i -e "s|huce_intel,seas_compute,shared|huce_cascade|g" \
    -e "s|--mem=15000|--mem=64000|g" \
    -e "s|-t 0-12:00|-t 07-00:00|g"\
    -e "s|-c 8|-c 24|g" geoschem.run
sbatch -W geoschem.run; wait;

# Write the boundary conditions using write_boundary_conditions.py
cd "${cwd}"
sbatch -W -J blended -o=boundary_conditions.log -p ${Partition} -t 7-00:00 --mem 96000 -c 48 --wrap "source ~/.bashrc; conda activate $CondaEnv; python write_boundary_conditions.py True $blendedDir"; wait; # run for Blended TROPOMI+GOSAT
sbatch -W -J tropomi -o=boundary_conditions.log -p ${Partition} -t 7-00:00 --mem 96000 -c 48 --wrap "source ~/.bashrc; conda activate $CondaEnv; python write_boundary_conditions.py False $tropomiDir"; wait; # run for TROPOMI data
echo "" >> "${cwd}/boundary_conditions.log"
echo "Blended TROPOMI+GOSAT boundary conditions --> {!$workDir}/blended-boundary-conditions" >> "${cwd}/boundary_conditions.log"
echo "TROPOMI boundary conditions               --> {!$workDir}/tropomi-boundary-conditions" >> "${cwd}/boundary_conditions.log"