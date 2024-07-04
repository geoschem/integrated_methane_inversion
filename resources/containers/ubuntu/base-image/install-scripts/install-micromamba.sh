#!/bin/bash

set -e

cd ~
# This component installs micromamba
# which is smaller and faster than conda
# uncomment below for intel processors
MicromambaInstallURL="https://micro.mamba.pm/api/micromamba/linux-64/latest"
# uncomment below for apple silicon
# MicromambaInstallURL="https://micro.mamba.pm/api/micromamba/linux-aarch64/latest"
curl -Ls $MicromambaInstallURL | tar -xvj bin/micromamba

# add micromamba to .bashrc
bin/micromamba shell init -s bash -p ~/micromamba

# set alias to use conda and micromamba interchangeably
echo $'conda() {\n    micromamba "$@"\n}' >> ~/.bashrc

# remove interactive check from .bashrc
sed -i '/# If not running interactively, don.t do anything/,/esac/d' ~/.bashrc

# install python dependencies for imi into conda environment
bin/micromamba env create -f /home/ubuntu/install-scripts/imi_env.yml
bin/micromamba clean -ya