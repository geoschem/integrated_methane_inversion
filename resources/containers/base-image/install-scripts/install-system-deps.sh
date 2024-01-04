#!/bin/bash

set -e # exit 1 if error

yum -y update
yum install -y vim wget time jq less glibc which sudo curl

# install aws-cli
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
./aws/install

# cleanup
rm awscliv2.zip
