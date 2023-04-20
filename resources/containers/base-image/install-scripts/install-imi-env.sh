#!/bin/bash

set -e
# install python dependencies for imi into conda environment
# Activate conda env
source /etc/bashrc
spack env activate compute_env

# install the environment
conda env create -f install-scripts/imi_env.yml
conda activate imi_env