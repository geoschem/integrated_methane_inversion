#!/bin/bash

# InstallConda This component install miniconda

# URL of miniconda installer
MinicondaInstallURL="https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh"

curl --silent $MinicondaInstallURL --output Miniconda3-Linux-x86_64.sh
bash Miniconda3-Linux-x86_64.sh -b -p /opt/conda
echo "PATH=/opt/conda/bin:$PATH" >> /etc/bashrc
echo ". /opt/conda/bin/activate" >> /etc/bashrc
PATH=/opt/conda/bin:$PATH
