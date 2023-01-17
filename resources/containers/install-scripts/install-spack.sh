#!/bin/bash

# name: InstallSpackEnvironment

# Spack install spec for desired compiler
SpackCompiler="gcc@10.2.0"
# Name of spack env
SpackEnvironmentName="compute_env"
# Spack environment manifest file
SpackEnvironmentFile="geoschem_deps-gnu-openmpi-102.yaml"
# InstallSpack
source /etc/bashrc
conda activate py38
yum install -y git curl texinfo
yum groupinstall -y "Development tools"
set -x
set -e # exit 1 if error
umask 022
git clone https://github.com/spack/spack /opt/spack
SPACK_ROOT=/opt/spack
. /opt/spack/share/spack/setup-env.sh
spack compiler find --scope system
spack external find --scope system
spack install $SpackCompiler
spack load $SpackCompiler
spack compiler find --scope system

# Add spack setup to bashrc
echo "export SPACK_ROOT=/opt/spack" >> /etc/bashrc
echo . /opt/spack/share/spack/setup-env.sh >> /etc/bashrc