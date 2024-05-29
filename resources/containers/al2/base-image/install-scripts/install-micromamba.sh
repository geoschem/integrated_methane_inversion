#!/bin/bash

set -e

# This component installs micromamba
# which is smaller and faster than conda
MicromambaInstallURL="https://micro.mamba.pm/api/micromamba/linux-64/latest"
yum install -y bzip2 tar
mkdir -p /opt/micromamba/bin/
curl -Ls $MicromambaInstallURL | tar -xvj bin/micromamba
mv bin/micromamba /opt/micromamba/bin/micromamba
rm -rf bin

# Configuration in bashrc file
echo "export MAMBA_ROOT_PREFIX=/opt/micromamba" >> /etc/bashrc
echo 'eval "$(/opt/micromamba/bin/micromamba shell hook -s posix)"' >> /etc/bashrc

# set alias to use conda and micromamba interchangeably
echo $'conda() {\n    micromamba "$@"\n}' >> /etc/bashrc
