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
    if "$isRegional"; then
        LandCoverFile="${DataPath}/GEOS_${gridDir}/${metDir}/${constYr}/01/${Met}.${constYr}0101.CN.${gridFile}.${RegionID}.${LandCoverFileExtension}"
    else
        LandCoverFile="${DataPath}/GEOS_${gridDir}/${metDir}/${constYr}/01/${Met}.${constYr}0101.CN.${gridFile}.${LandCoverFileExtension}"
    fi

    # Use archived HEMCO standalone emissions output
    HemcoDiagFile="${DataPath}/HEMCO/CH4/v2024-07/HEMCO_SA_Output/HEMCO_sa_diagnostics.${gridFile}.2023.nc"

    if "$isAWS"; then
        # Download land cover and HEMCO diagnostics files
        s3_lc_path="s3://gcgrid/GEOS_${gridDir}/${metDir}/${constYr}/01/${Met}.${constYr}0101.CN.${gridFile}.${RegionID}.${LandCoverFileExtension}"
        aws s3 cp --no-sign-request ${s3_lc_path} ${LandCoverFile}
        s3_hd_path="s3://gcgrid/HEMCO/CH4/v2024-07/HEMCO_SA_Output/HEMCO_sa_diagnostics.${gridFile}.2023.nc"
        aws s3 cp --no-sign-request ${s3_hd_path} ${HemcoDiagFile}
    fi

    # Output path and filename for state vector file
    StateVectorFName="StateVector.nc"

    # Create state vector file
    cd ${RunDirs}

    # Copy state vector creation script to working directory
    cp ${InversionPath}/src/utilities/make_state_vector_file.py .
    chmod 755 make_state_vector_file.py

    printf "\nCalling make_state_vector_file.py\n"
    python make_state_vector_file.py $ConfigPath $LandCoverFile $HemcoDiagFile $StateVectorFName

    printf "\n=== DONE CREATING RECTANGULAR STATE VECTOR FILE ===\n"
}

# Description: Reduce dimension of state vector with clustering method
# Usage:
#   reduce_dimension
reduce_dimension() {
    printf "\n=== REDUCING DIMENSION OF STATE VECTOR FILE ===\n"

    # set input variables
    state_vector_path=${RunDirs}/StateVector.nc
    native_state_vector_path=${RunDirs}/NativeStateVector.nc
    preview_dir=${RunDirs}/preview
    tropomi_cache=${RunDirs}/satellite_data
    aggregation_file=${InversionPath}/src/components/statevector_component/aggregation.py

    if [[ ! -f ${RunDirs}/NativeStateVector.nc ]]; then
        # copy the original state vector file for subsequent statevector generations
        printf "\nCopying native state vector file to NativeStateVector.nc \n"
        cp $state_vector_path $native_state_vector_path
    else
        # replace state vector file with clean, native resolution state vector
        cp $native_state_vector_path $state_vector_path
    fi

    # conditionally add period_i to python args
    python_args=($aggregation_file $InversionPath $ConfigPath $state_vector_path $preview_dir $tropomi_cache)
    archive_sv=false
    if ("$KalmanMode" && "$DynamicKFClustering"); then
        if [ -n "$period_i" ]; then
            archive_sv=true
            python_args+=($period_i)
        fi
    fi

    # if running end to end script with sbatch then use
    # sbatch to take advantage of multiple cores
    if "$UseSlurm"; then
        rm -f .aggregation_error.txt
        chmod +x $aggregation_file
        sbatch --mem $RequestedMemory \
            -c $RequestedCPUs \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -o imi_output.tmp \
            -W "${python_args[@]}"
        wait
        cat imi_output.tmp >>${InversionPath}/imi_output.log
        rm imi_output.tmp
        # check for any errors
        [ ! -f ".aggregation_error.txt" ] || imi_failed $LINENO
    else
        python "${python_args[@]}"
    fi

    # archive state vector file if using Kalman filter
    if "$archive_sv"; then
        mkdir -p ${RunDirs}/archive_sv
        cp $state_vector_path ${RunDirs}/archive_sv/StateVector_${period_i}.nc
    fi
    nElements=$(ncmax StateVector ${RunDirs}/StateVector.nc)
    if "$OptimizeBCs"; then
        nElements=$((nElements + 4))
    fi
    if "$OptimizeOH"; then
        if "$isRegional"; then
            nElements=$((nElements + 1))
        else
            nElements=$((nElements + 2))
        fi
    fi
    printf "\nNumber of state vector elements in this inversion = ${nElements}\n\n"
    printf "\n=== DONE REDUCING DIMENSION OF STATE VECTOR FILE ===\n"
}
