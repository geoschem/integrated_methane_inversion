#!/bin/bash

# install dependencies for gcpy into conda environment

# Activate conda env
source /etc/bashrc
spack env activate compute_env

# install the environment
conda config --add channels conda-forge
conda env create -f install-scripts/gcpy_environment.yaml
# esmpy 8.4.0 is not compatible with xesmf yet
conda activate py38
conda install -y -c conda-forge esmpy">=8.0.0, <8.4.0"
conda install -c conda-forge sparselt==0.1.3