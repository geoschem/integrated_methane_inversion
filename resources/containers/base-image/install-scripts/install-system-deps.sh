#!/bin/bash

set -e # exit 1 if error

yum install -y emacs wget time jq less glibc which sudo

# install aws-cli
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
wget https://github.com/mikefarah/yq/releases/download/v4.32.2/yq_linux_amd64 -O /opt/conda/envs/py39/bin/yq  \
    && chmod +x /opt/conda/envs/py39/bin/yq
./aws/install

# cleanup
rm awscliv2.zip
