#!/bin/bash

# name: InstallSpackEnvironment

# IMPORTANT: run 'spack load gcc@12.2.0' on the command line prior to running this script.

SpackEnvironmentName="compute_env"
SpackCompiler="gcc@12.2.0"
SpackEnvFile="geoschem_deps-gnu-openmpi-122.yml"
InstallScriptPath="/home/ubuntu/integrated_methane_inversion/resources/aws/gchp/ubuntu/install-scripts"

# Load spack
source /home/ubuntu/.bashrc

# Install spack environment
umask 022
spack env create $SpackEnvironmentName $InstallScriptPath/$SpackEnvFile
spack env activate $SpackEnvName
spack install -j32 --fail-fast --show-log-on-error
spack clean --all

# garbage collector to reduce the size of spack install
spack gc -Ey

# AutomaticallyLoadEnvironmentOnLogin
echo spack env activate $SpackEnvName >> ~/.bashrc
