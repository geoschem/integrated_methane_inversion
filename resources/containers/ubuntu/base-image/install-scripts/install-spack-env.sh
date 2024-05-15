#!/bin/bash

# name: InstallSpackEnvironment

set -e # exit 1 if error
set -x
# Name of spack env
SpackEnvironmentName="compute_env"
# Spack environment manifest file
SpackEnvironmentFile="geoschem_deps-gnu-openmpi-122.yml"

# Load spack
. ~/spack/share/spack/setup-env.sh

# name: InstallSpackEnvironment
echo "Downloading spack environment file"
umask 022

spack compiler list

spack env create $SpackEnvironmentName /home/ubuntu/install-scripts/$SpackEnvironmentFile
spack env activate $SpackEnvironmentName
spack install --fail-fast --show-log-on-error
spack clean --all

# garbage collector to reduce the size of spack install
spack gc -Ey
      
# AutomaticallyLoadEnvironmentOnLogin
echo spack env activate $SpackEnvironmentName >> ~/.bashrc
