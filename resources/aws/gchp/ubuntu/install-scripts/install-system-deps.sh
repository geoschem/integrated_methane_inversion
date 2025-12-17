#!/bin/bash

set -e # exit 1 if error

# Retrieve system's local package index
sudo apt update

# Install core utilities
sudo apt install emacs unzip bzip2 flex bison -y
sudo apt install build-essential gcc g++ gfortran -y

# Install aws-cli
cd /home/ubuntu
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
rm awscliv2.zip
sudo ./aws/install

#yum install -y emacs wget time jq less glibc which
#
## install yq
#mkdir -p /opt/conda/envs/py39/bin/
#wget https://github.com/mikefarah/yq/releases/download/v4.32.2/yq_linux_amd64 -O /opt/conda/envs/py39/bin/yq
#chmod +x /opt/conda/envs/py39/bin/yq
#



