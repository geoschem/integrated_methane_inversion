#!/bin/bash

set -e 

# modify permissions for various slurm related dependencies
sudo chmod 777 /home/ubuntu/install-scripts
sudo chmod 777 /etc/slurm
sudo mkdir -p /run/munge
sudo chown munge:munge /run/munge

# create munge key
sudo -u munge rm /etc/munge/munge.key
sudo -u munge /bin/bash -c 'umask 077 && dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key'

# make necessary directories
sudo mkdir -p /var/spool/slurmctld
sudo mkdir -p /var/spool/slurmd
sudo mkdir -p /bin/mail

# assign ownership to ubuntu user
sudo chown -R ubuntu:ubuntu /var/spool/slurmctld
sudo chown -R ubuntu:ubuntu /var/spool/slurmd
sudo chown -R ubuntu:ubuntu /bin/mail

# create some default files
sudo touch /var/log/slurmctld.log
sudo chown ubuntu:ubuntu /var/log/slurmctld.log
sudo touch /var/log/slurmd.log
sudo chown ubuntu:ubuntu /var/log/slurmd.log
sudo touch /var/run/slurmd.pid
sudo chown ubuntu:ubuntu /var/run/slurmd.pid
sudo -u ubuntu touch /var/spool/slurmctld/node_state
sudo -u ubuntu touch /var/spool/slurmctld/job_state
sudo -u ubuntu touch /var/spool/slurmctld/job_state.old
sudo -u ubuntu touch /var/spool/slurmctld/resv_state
sudo -u ubuntu touch /var/spool/slurmctld/resv_state.old
sudo -u ubuntu touch /var/spool/slurmctld/trigger_state
sudo -u ubuntu touch /var/spool/slurmctld/trigger_state.old