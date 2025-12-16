#!/bin/bash

# From gc-cloud-benchmark
# name: InstallSpackEnvironment

# Name of spack env
SpackEnvironmentName="compute_env"

# Spack install spec for desired compiler
SpackCompiler="gcc@12.2.0"

# Spack version to install
SpackVersion="releases/v0.23"

# InstallSpack
source ~/.bashrc
set -x
set -e # exit 1 if error

# Spack needs python 3.6 or higher
#micromamba activate imi_env
#conda create -n pysetup python=3.12.0 -c anaconda
#conda activate pysetup

# install basic dependencies
##sudo apt install -y git curl curl-devel texinfo
#sudo apt install -y git curl
#sudo apt install -y git texinfo
##sudo apt groupinstall -y "Development tools"
umask 022

# install spack
git clone -c feature.manyFiles=true --depth=2 -b $SpackVersion https://github.com/spack/spack /home/ubuntu/spack
cd /home/ubuntu/spack
SPACK_ROOT=/home/ubuntu/spack
. /home/ubuntu/spack/share/spack/setup-env.sh

#spack compiler find --scope system
#spack external find --scope system
#spack install -j32 --fail-fast $SpackCompiler --show-log-on-error target=x86_64 platform=linux os=amzn2
#spack load $SpackCompiler
#spack compiler find --scope system
#
## Add spack setup to bashrc
#echo "export SPACK_ROOT=/opt/spack" >> ~/.bashrc
#echo . /home/ubuntu/share/spack/setup-env.sh >> ~/.bashrc

# From imi:
##git clone --depth=100 --branch=v0.22.1 https://github.com/spack/spack.git ~/spack
##cd ~/spack
##
### make etc/spack directory
##sudo mkdir /etc/spack
##sudo chown -R ubuntu:ubuntu /etc/spack
##
### start spack
##. share/spack/setup-env.sh
##
### find external packages/ compiler
##spack compiler find --scope system
##spack external find --scope system
##
##echo . ~/spack/share/spack/setup-env.sh >> ~/.bashrc
