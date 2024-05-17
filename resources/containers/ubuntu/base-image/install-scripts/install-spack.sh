#!/bin/bash

git clone --depth=100 --branch=releases/v0.22 https://github.com/spack/spack.git ~/spack
cd ~/spack

# make etc/spack directory
sudo mkdir /etc/spack
sudo chown -R ubuntu:ubuntu /etc/spack

# start spack
. share/spack/setup-env.sh

# find external packages/ compiler
spack compiler find --scope system
spack external find --scope system

echo . ~/spack/share/spack/setup-env.sh >> ~/.bashrc