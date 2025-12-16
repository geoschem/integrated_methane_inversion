#!/bin/bash

set -e # exit 1 if error

# Retrieve system's local package index
sudo apt update

# Install core utilities
sudo apt install -y emacs
sudo apt install -y unzip
sudo apt install -y bzip2

#yum install -y emacs wget time jq less glibc which
#
## install yq
#mkdir -p /opt/conda/envs/py39/bin/
#wget https://github.com/mikefarah/yq/releases/download/v4.32.2/yq_linux_amd64 -O /opt/conda/envs/py39/bin/yq
#chmod +x /opt/conda/envs/py39/bin/yq
#

# Install aws-cli
cd /home/ubuntu
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
rm awscliv2.zip

