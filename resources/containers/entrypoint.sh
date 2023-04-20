#!/usr/bin/env bash
source /etc/bashrc
spack unload -a

# start munge service
munged

# initialize slurm by dynamically generating slurm.conf
cd /home/al2/install-scripts
cp /home/al2/install-scripts/base_slurm.conf /home/al2/install-scripts/new_slurm.conf
python /home/al2/install-scripts/configure_slurm.py
# put newly generated slurm.conf in correct locations
cp /home/al2/install-scripts/new_slurm.conf /usr/local/etc/slurm.conf
rm /home/al2/install-scripts/new_slurm.conf
# also add a cgroup.conf file
cp /home/al2/install-scripts/cgroup.conf /usr/local/etc/cgroup.conf
# start slurm
slurmctld
slurmd

# Fix issue when switching instance types where node claims to be drained
# scontrol update nodename=$HOSTNAME state=idle
# scontrol update nodename=$HOSTNAME state=DOWN reason="undraining"
# scontrol update nodename=$HOSTNAME state=RESUME

while :
do
	echo "Running Forever. Shut down manually to stop"
	sleep 10
done