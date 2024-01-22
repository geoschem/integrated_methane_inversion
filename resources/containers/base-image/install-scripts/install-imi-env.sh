#!/bin/bash

set -e
# install python dependencies for imi into conda environment
# Activate micromamba and spack
source /etc/bashrc
spack env activate compute_env

# install the environment
micromamba env create -f install-scripts/imi_env.yml
micromamba clean -a
micromamba activate imi_env

# Also install yq into this environment
wget https://github.com/mikefarah/yq/releases/download/v4.32.2/yq_linux_amd64 -O /opt/micromamba/envs/imi_env/bin/yq  \
    && chmod +x /opt/micromamba/envs/imi_env/bin/yq