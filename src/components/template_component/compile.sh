#!/bin/bash

# Functions available in this file include:
#   - compile_geoschem

# Description: Compile GEOS-Chem, either GC-Classic or GCHP 
# Usage:
#   compile_geoschem
compile_geoschem() {

    # Compile GEOS-Chem/GCHP and store in <src>/build/bin directory
    if "$UseGCHP"; then
	printf "\n=== COMPILING GCHP ===\n"
        cd ${GCHPPath}/run
        execname="gchp"
        build_dir=${GCHPPath}/build
        KPP_dir=${GCHPPath}/src/GCHP_GridComp/GEOSChem_GridComp/geos-chem/KPP
    else
	printf "\n=== COMPILING GC-Classic ===\n"
        cd ${GCClassicPath}/run
        execname="gcclassic"
        build_dir=${GCClassicPath}/build
        KPP_dir=${GCClassicPath}/src/GEOS-Chem/KPP
    fi

    mkdir -p ${build_dir}
    cd ${build_dir}
    # first build a default exexutable without Jacobian tracers
    if [ -f "bin/${execname}.default" ]; then
        echo "Executable bin/${execname}.default already exists — skipping rebuild."
    else
        echo "Building ${execname}.default ..."

        # remove CMakeCache.txt once to initialize JACOBIAN cmake option
        rm -f CMakeCache.txt
        cd "${KPP_dir}/carbon"
        if [ ! -f carbon.eqn.default ]; then
            mv carbon.eqn carbon.eqn.default
        fi
        ln -nsf carbon.eqn.default carbon.eqn
        # initialize KPP mechanism to be the default carbon equations
        cd ${KPP_dir}
        ./build_mechanism.sh carbon >> "${build_dir}/build_geoschem.log" 2>&1 \
            || { echo "ERROR: build_mechanism.sh carbon failed." >&2; exit 1; }

        cd "${build_dir}"
        cmake .. >> build_geoschem.log 2>&1
        cmake . -DMECH=carbon -DJACOBIAN=n >> build_geoschem.log 2>&1
        make -j >> build_geoschem.log 2>&1

        mv "bin/${execname}" "bin/${execname}.default"
    fi

    # then expand a series of carbon equations for Jacobian tracers
    # sanity check on NumJacobianTracers
    if ! [[ "$NumJacobianTracers" =~ ^[0-9]+$ ]] || [ "$NumJacobianTracers" -eq 0 ]; then
        echo "ERROR: NumJacobianTracers must be a positive integer > 0." >&2
        exit 1
    fi
    # determine number of executables to be built based on NumJacobianTracers
    # get the ceiling number rounded to 10s
    upper=$(( (NumJacobianTracers + 9) / 10 * 10 ))

    cd "${build_dir}"
    cmake . -DMECH=carbon -DJACOBIAN=y >> build_geoschem.log 2>&1
    for n in $(seq 10 10 $upper); do
        if [ -f "bin/${execname}.${n}" ]; then
            echo "Executable bin/${execname}.${n} already exists — skipping rebuild."
        else
            echo "Building ${execname}.${n} ..."

            cd ${KPP_dir}/carbon
            ${InversionPath}/src/utilities/expand_carbon_eqn.py \
                carbon.eqn.default ${n} > carbon.eqn.${n}
            ln -nsf carbon.eqn.${n} carbon.eqn
            # generate KPP carbon mechanism
            cd ${KPP_dir}
            ./build_mechanism.sh carbon >> "${build_dir}/build_geoschem.log" 2>&1 \
                || { echo "ERROR: build_mechanism.sh carbon failed." >&2; exit 1; }
            # build GCC/GCHP with expanded carbon mechanism for Jacobian tracers
            cd ${build_dir}
            make -j >>build_geoschem.log 2>&1
            mv bin/${execname} bin/${execname}.${n}
        fi
    done

    printf "\nDone compiling GEOS-Chem \n\nSee ${build_dir}/build_geoschem.log for details\n\n"

    # Navigate back to working directory
    cd ${RunDirs}

    printf "\n=== DONE BUILDING GEOS-Chem ===\n"
}
