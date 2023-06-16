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
    HemcoDiagFile="${DataPath}/HEMCO/CH4/v2023-04/HEMCO_SA_Output/HEMCO_sa_diagnostics.${gridRes}.20190101.nc"

    if "$isAWS"; then
	# Download land cover and Hemco diagnostics files
	s3_lc_path="s3://gcgrid/GEOS_${gridDir}/${metDir}/${constYr}/01/${metUC}.${constYr}0101.CN.${gridRes}.${NestedRegion}.${LandCoverFileExtension}"
	aws s3 cp --request-payer=requester ${s3_lc_path} ${LandCoverFile}
    s3_hd_path="s3://gcgrid/HEMCO/CH4/v2023-04/HEMCO_SA_Output/HEMCO_sa_diagnostics.${gridRes}.20190101.nc"
    aws s3 cp --request-payer=requester ${s3_hd_path} ${HemcoDiagFile}
    fi

    # Output path and filename for state vector file
    StateVectorFName="StateVector.nc"

    # Create state vector file
    cd ${RunDirs}

    # Copy state vector creation script to working directory
    cp ${InversionPath}/src/utilities/make_state_vector_file.py .
    chmod 755 make_state_vector_file.py

    # Get config path
    config_path=${InversionPath}/${ConfigFile}

    printf "\nCalling make_state_vector_file.py\n"
    python make_state_vector_file.py $config_path $LandCoverFile $HemcoDiagFile $StateVectorFName

    printf "\n=== DONE CREATING RECTANGULAR STATE VECTOR FILE ===\n"
}

# Description: Reduce dimension of state vector with clustering method
# Usage:
#   reduce_dimension
reduce_dimension() {
    printf "\n=== REDUCING DIMENSION OF STATE VECTOR FILE ===\n"

    # First run the Preview if necessary to get prior emissions
    if [[ ! -d ${RunDirs}/preview_run/OutputDir ]]; then
        printf "\nPreview Dir not detected. Running the IMI Preview as a prerequisite.\n"
        run_preview
    fi

    # set input variables
    config_path=${InversionPath}/${ConfigFile}
    state_vector_path=${RunDirs}/StateVector.nc
    native_state_vector_path=${RunDirs}/NativeStateVector.nc

    preview_dir=${RunDirs}/preview_run
    tropomi_cache=${RunDirs}/data_TROPOMI
    aggregation_file=${InversionPath}/src/components/statevector_component/aggregation.py

    if [[ ! -f ${RunDirs}/NativeStateVector.nc ]]; then
        # copy the original state vector file for subsequent statevector generations
        printf "\nCopying native state vector file to NativeStateVector.nc \n"
        cp $state_vector_path $native_state_vector_path
    else
        # replace state vector file with clean, native resolution state vector
        cp $native_state_vector_path $state_vector_path
    fi

    # if running end to end script with sbatch then use
    # sbatch to take advantage of multiple cores 
    export PYTHONPATH=${PYTHONPATH}:${InversionPath}/src/
    export PYTHONPATH=${PYTHONPATH}:${InversionPath}/src/inversion_scripts
    if "$UseSlurm"; then
        chmod +x $aggregation_file
        sbatch --mem $SimulationMemory \
        -c $SimulationCPUs \
        -t $RequestedTime \
        -p $SchedulerPartition \
        -W $aggregation_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache; wait;
    else
        python $aggregation_file $InversionPath $config_path $state_vector_path $preview_dir $tropomi_cache
    fi
    nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
    printf "\nNumber of state vector elements in this inversion = ${nElements}\n\n"
    printf "\n=== DONE REDUCING DIMENSION OF STATE VECTOR FILE ===\n"
}