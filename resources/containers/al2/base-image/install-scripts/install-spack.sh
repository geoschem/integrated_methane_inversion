#!/bin/bash

# name: InstallSpackEnvironment
set -x
set -e # exit 1 if error

# Spack version to install
# TODO: replace with SpackVersion="v0.22.0" when released
SpackVersion="561da58cea1475927b22b95aa7a6b567ef1105f3"
# Spack install spec for desired compiler
SpackCompiler="gcc@10.2.0"
# Name of spack env
SpackEnvironmentName="compute_env"
# Spack environment manifest file
SpackEnvironmentFile="geoschem_deps-gnu-openmpi-102.yaml"
# InstallSpack
source /etc/bashrc

# Spack needs python 3.6 or higher
micromamba create -n py39 python=3.9.1 -c anaconda
micromamba activate py39

yum install -y git curl curl-devel texinfo
yum groupinstall -y "Development tools"
umask 022
# Use spack v0.21.0
git clone https://github.com/spack/spack /opt/spack
(cd /opt/spack && git checkout $SpackVersion)

SPACK_ROOT=/opt/spack
. /opt/spack/share/spack/setup-env.sh
spack compiler find --scope system
spack external find --scope system
spack install --fail-fast $SpackCompiler
spack load $SpackCompiler
spack compiler find --scope system

# Add spack setup to bashrc ensuring we activate py39 first
echo "micromamba activate py39" >> /etc/bashrc
echo "export SPACK_ROOT=/opt/spack" >> /etc/bashrc
echo . /opt/spack/share/spack/setup-env.sh >> /etc/bashrc