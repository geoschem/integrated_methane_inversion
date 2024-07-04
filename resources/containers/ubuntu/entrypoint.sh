#!/bin/bash

set -e 

# configure slurm
cd /home/ubuntu/install-scripts
cp /home/ubuntu/install-scripts/base_slurm.conf /home/ubuntu/install-scripts/new_slurm.conf
python3 /home/ubuntu/install-scripts/configure_slurm.py

# put newly generated slurm.conf in correct locations
cp /home/ubuntu/install-scripts/new_slurm.conf /etc/slurm/slurm.conf
rm /home/ubuntu/install-scripts/new_slurm.conf

# also add a cgroup.conf file
cp /home/ubuntu/install-scripts/cgroup.conf /etc/slurm/cgroup.conf

# start munge and slurm services
sudo -u munge service munge start
service slurmctld start
service slurmd start

# run test script
cd /home/ubuntu/integrated_methane_inversion
mkdir -p ~/ExtData
mkdir -p ~/imi_output_dir

# load micromamba
source ~/.bashrc
micromamba activate imi_env

# Test if using a non-default config file
if [[ "x${IMI_CONFIG_PATH}" != "x" ]]; then
	config_file=$IMI_CONFIG_PATH
else
	config_file="config.yml"
fi

# override specific config file vars with env variables of 
# the syntax IMI_<config-variable>
chmod +x src/utilities/override_config_variables.py
python src/utilities/override_config_variables.py $config_file $config_file

# runs an imi preview
./run_imi.sh | tee imi_output.log