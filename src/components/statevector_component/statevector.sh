#!/bin/bash

# Functions available in this file include:
#   - create_statevector 
#   - reduce_dimension

# Description: Create a native resolution state vector
# Usage:
#   create_statevector
create_statevector() {
    printf "\n=== CREATING RECTANGULAR STATE VECTOR FILE ===\n"
    
    # Use GEOS-FP or MERRA-2 CN file to determine ocean/land grid boxes
    LandCoverFile="${DataPath}/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${NestedRegion}.${LandCoverFileExtension}"

    if "$isAWS"; then
	# Download land cover file
	s3_lc_path="s3://gcgrid/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${NestedRegion}.${LandCoverFileExtension}"
	aws s3 cp --request-payer=requester ${s3_lc_path} ${LandCoverFile}
    fi

    # Output path and filename for state vector file
    StateVectorFName="StateVector.nc"

    # Create state vector file
    cd ${RunDirs}

    # Copy state vector creation script to working directory
    cp ${InversionPath}/src/utilities/make_state_vector_file.py .
    chmod 755 make_state_vector_file.py

    printf "\nCalling make_state_vector_file.py\n"
    python make_state_vector_file.py $LandCoverFile $StateVectorFName $LatMin $LatMax $LonMin $LonMax $BufferDeg $LandThreshold $nBufferClusters

    printf "\n=== DONE CREATING RECTANGULAR STATE VECTOR FILE ===\n"
}

# Description: Reduce dimension of state vector with clustering method
# Usage:
#   reduce_dimension
reduce_dimension() {
    # First run the Preview
    # run_preview
    # Make sure template run directory exists
    if [[ ! -f ${RunTemplate}/geoschem_config.yml ]]; then
        printf "\nTemplate run directory does not exist or has missing files. Please set 'SetupTemplateRundir=true' in config.yml\n" 
        exit 9999
    fi

    # Run preview script
    config_path=${InversionPath}/${ConfigFile}
    state_vector_path=${RunDirs}/StateVector.nc
    preview_dir=${RunDirs}/preview_run
    tropomi_cache=${RunDirs}/data_TROPOMI
    aggregation_file=${InversionPath}/src/components/statevector_component/aggregation.py
    # if running end to end script with sbatch then use
    # sbatch to take advantage of multiple cores 
    if "$UseSlurm"; then
        # set number of cores to run preview with
        if "$isAWS"; then
            sed -i -e "s:#SBATCH -c 8:#SBATCH -c ${cpu_count}:g" ${InversionPath}/src/inversion_scripts/imi_preview.py
        fi
        export PYTHONPATH=${PYTHONPATH}:${InversionPath}/src/
        export PYTHONPATH=${PYTHONPATH}:${InversionPath}/src/inversion_scripts
        chmod +x $aggregation_file
        sbatch -W $aggregation_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache; wait;
    else
        python $aggregation_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache
    fi
}