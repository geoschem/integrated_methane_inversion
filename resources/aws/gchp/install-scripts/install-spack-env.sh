#!/bin/bash

# name: InstallSpackEnvironment

# Name of spack env
SpackEnvironmentName="compute_env"

# Spack install spec for desired compiler
SpackCompiler="gcc@12.2.0"

# Spack environment manifest file
SpackEnvironmentFile="geoschem_deps-gnu-openmpi-122.yml"

# Load spack
source /home/ubuntu/.bashrc

# Clean up, in case debugging
spack env rm $SpackEnvironmentName

# Install spack environment
umask 022
set -x
set -e # exit 1 if error
cat /home/ubuntu/install-scripts/$SpackEnvironmentFile
spack compiler list


spack env create $SpackEnvironmentName /home/ubuntu/install-scripts/$SpackEnvironmentFile
spack load $SpackCompiler
spack env activate $SpackEnvironmentName
spack install -j32 --fail-fast --show-log-on-error target=x86_64 platform=linux os=amzn2
spack clean --all
##      
## AutomaticallyLoadEnvironmentOnLogin
##echo spack load $SpackCompiler >> /etc/bashrc
##echo spack env activate $SpackEnvironmentName >> /etc/bashrc
##chmod -R +r /opt/spack
#

#spack env create $SpackEnvironmentName /home/ubuntu/install-scripts/$SpackEnvironmentFile
#spack env activate $SpackEnvironmentName
#spack install --fail-fast --show-log-on-error
#spack clean --all
#
## garbage collector to reduce the size of spack install
#spack gc -Ey
#      
## AutomaticallyLoadEnvironmentOnLogin
#echo spack env activate $SpackEnvironmentName >> ~/.bashrc
