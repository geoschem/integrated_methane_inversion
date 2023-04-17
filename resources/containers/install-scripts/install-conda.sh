#!/bin/bash

set -e
# InstallConda This component install miniconda
MinicondaInstallURL="https://repo.anaconda.com/miniconda/Miniconda3-py39_23.1.0-1-Linux-x86_64.sh"
# URL of miniconda installer
curl --silent $MinicondaInstallURL --output Miniconda3-Linux-x86_64.sh
bash Miniconda3-Linux-x86_64.sh -b -p /opt/conda
echo "PATH=/opt/conda/bin:$PATH" >> /etc/bashrc
echo ". /opt/conda/bin/activate" >> /etc/bashrc
PATH=/opt/conda/bin:$PATH
