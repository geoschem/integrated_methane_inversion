#!/bin/bash

# name: InstallSpackEnvironment

set -e # exit 1 if error
set -x
# Spack install spec for desired compiler
SpackCompiler="gcc@10.2.0"
# Name of spack env
SpackEnvironmentName="compute_env"
# Spack environment manifest file
SpackEnvironmentFile="geoschem_deps-gnu-openmpi-102.yaml"
# InstallSpack
source /etc/bashrc

# name: InstallSpackEnvironment
echo "Downloading spack environment file"
umask 022
# curl $SpackEnvironmentFile -o /spack-manifest.yaml
# . /opt/spack/share/spack/setup-env.sh
cat /home/al2/install-scripts/$SpackEnvironmentFile
spack compiler list

spack env create $SpackEnvironmentName /home/al2/install-scripts/$SpackEnvironmentFile
spack load $SpackCompiler
spack env activate $SpackEnvironmentName
spack install
spack clean --all
      
# AutomaticallyLoadEnvironmentOnLogin
echo spack load $SpackCompiler >> /etc/bashrc
echo spack env activate $SpackEnvironmentName >> /etc/bashrc
chmod -R +r /opt/spack
