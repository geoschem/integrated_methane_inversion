#!/bin/bash

# Functions available in this file include:
#   - prepare_statevector_inputs (helper function)
#   - create_statevector
#   - reduce_dimension
#   - regrid_statevector

# Description: get the input for creating or regridding state vector 
# Usage: read LandCoverFile HemcoDiagFile < <(prepare_statevector_inputs)
prepare_statevector_inputs() {
    # Use GEOS-FP or MERRA-2 CN file to determine ocean/land grid boxes
    if "$isRegional"; then
        if [ "$Res" = "0.125x0.15625" ]; then
            LandCoverSuffix="HEMCO/CH4/v2025-03/landcover/IMERG_land_sea_mask_0125x015625.nc"
        else
            if [ "$RegionID" != "" ]; then
                LandCoverSuffix="GEOS_${gridDir}/${metDir}/${constYr}/01/${Met}.${constYr}0101.CN.${gridFile}.${RegionID}.${LandCoverFileExtension}"
            else
                LandCoverSuffix="GEOS_${gridDir}/${metDir}/${constYr}/01/${Met}.${constYr}0101.CN.${gridFile}.${LandCoverFileExtension}"
            fi
        fi
    else
        LandCoverSuffix="GEOS_${gridDir}/${metDir}/${constYr}/01/${Met}.${constYr}0101.CN.${gridFile}.${LandCoverFileExtension}"
    fi
    LandCoverFile="${DataPath}/${LandCoverSuffix}"

    # Use archived HEMCO standalone emissions output
    HemcoDiagFile="${DataPath}/HEMCO/CH4/v2025-07/HEMCO_SA_Output/HEMCO_sa_diagnostics.${gridFile}.2023.nc"

    # Ensure files exist or download them
    if [ ! -f "$LandCoverFile" ]; then
        s3_lc_path="s3://gcgrid/${LandCoverSuffix}"
        python "${InversionPath}/src/utilities/download_aws_file.py" "$s3_lc_path" "$LandCoverFile"
    fi
    if [ ! -f "$HemcoDiagFile" ]; then
        s3_hd_path="s3://gcgrid/HEMCO/CH4/v2025-07/HEMCO_SA_Output/HEMCO_sa_diagnostics.${gridFile}.2023.nc"
        python "${InversionPath}/src/utilities/download_aws_file.py" "$s3_hd_path" "$HemcoDiagFile"
    fi

    # for downstream use
    echo "$LandCoverFile $HemcoDiagFile"
}

# Description: Create a native resolution state vector
# Usage:
#   create_statevector
create_statevector() {
    if "$UseGCHP"; then
        if "$STRETCH_GRID"; then
            printf "\n=== CREATING Cubed-Sphere C${CS_RES}.s${STRETCH_FACTOR}_${TARGET_LAT}N_${TARGET_LON}E STATE VECTOR FILE ===\n"
        else
            printf "\n=== CREATING Cubed-Sphere C${CS_RES} STATE VECTOR FILE ===\n"
        fi
    else
        printf "\n=== CREATING RECTANGULAR STATE VECTOR FILE ===\n"
    fi

    # get the input path for state vector
    read LandCoverFile HemcoDiagFile < <(prepare_statevector_inputs)
    # Output path and filename for state vector file
    StateVectorFName="StateVector.nc"

    # Create state vector file
    cd ${RunDirs}

    # Copy state vector creation script to working directory
    cp ${InversionPath}/src/utilities/make_state_vector_file.py .
    chmod 755 make_state_vector_file.py

    printf "\nCalling make_state_vector_file.py\n"
    python make_state_vector_file.py $ConfigPath $LandCoverFile $HemcoDiagFile $StateVectorFName

    if "$UseGCHP"; then
        if "$STRETCH_GRID"; then
            printf "\n=== DONE CREATING Cubed-Sphere C${CS_RES}.s${STRETCH_FACTOR}_${TARGET_LAT}N_${TARGET_LON}E STATE VECTOR FILE ===\n"
        else
            printf "\n=== DONE CREATING Cubed-Sphere C${CS_RES} STATE VECTOR FILE ===\n"
        fi
    else
        printf "\n=== DONE CREATING RECTANGULAR STATE VECTOR FILE ===\n"
    fi
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
    python_args=($aggregation_file $ConfigPath $native_state_vector_path $state_vector_path $preview_dir $tropomi_cache)
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
        AggCPUs="${InversionCPUs:-$RequestedCPUs}"
        sbatch --mem $RequestedMemory \
            -c $AggCPUs \
            -t $RequestedTime \
            -p $SchedulerPartition \
            -o imi_output.tmp \
            -W "${python_args[@]}"
        wait
        cat imi_output.tmp >>${RunDirs}/imi_output.log
        rm imi_output.tmp
        # check for any errors
        [ ! -f ".aggregation_error.txt" ] || imi_failed $LINENO
    else
        python "${python_args[@]}" || {
            echo "ERROR: aggregation.py failed"
            return 1
        }
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

regrid_statevector(){
    printf "\n=== REGRID STATE VECTOR at ${ReferenceStateVectorFile} to CURRENT GRID ===\n"
    
    # get the input path for state vector
    read LandCoverFile HemcoDiagFile < <(prepare_statevector_inputs)
    # Output path and filename for state vector file
    StateVectorFName="StateVector.nc"

    # Create state vector file
    cd ${RunDirs}

    # Copy state vector regriddiing script to working directory
    cp ${InversionPath}/src/utilities/regrid_state_vector_file.py .
    
    printf "\nCalling regrid_state_vector_file.py\n"
    python regrid_state_vector_file.py $ConfigPath $LandCoverFile $HemcoDiagFile $StateVectorFName

    printf "\n=== DONE REGRID STATE VECTOR at ${ReferenceStateVectorFile} to CURRENT GRID ===\n"
}