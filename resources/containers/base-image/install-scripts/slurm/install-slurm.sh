#!/bin/bash
# install openssl *dependency for munge* 

## install munge *may not need for single container*
# set -e
# git clone https://github.com/dun/munge.git
# cd munge
# ./bootstrap
# ./configure --prefix=/usr --sysconfdir=/etc --localstatedir=/var --with-runstatedir=/run
# make
# make check
# sudo make install

# create slurm user and munge user
# export MUNGEUSER=991
# groupadd -g $MUNGEUSER munge
# useradd  -m -c "MUNGE Uid 'N' Gid Emporium" -d /var/lib/munge -u $MUNGEUSER -g munge  -s /sbin/nologin munge
# export SLURMUSER=992
# groupadd -g $SLURMUSER slurm
# useradd  -m -c "SLURM workload manager" -d /var/lib/slurm -u $SLURMUSER -g slurm  -s /bin/bash slurm

# install munge
set -e 
amazon-linux-extras install epel -y
yum install munge munge-libs munge-devel -y
yum install rng-tools -y
rngd -r /dev/urandom
/usr/sbin/create-munge-key -r
dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key
chown -R root:root /var/log/munge
chown root:root /run/munge
chown root:root /var/lib/munge 
chown root:root /etc/munge
chown root:root /etc/munge/munge.key
# run with $ munged

# install slurm
curl https://download.schedmd.com/slurm/slurm-23.02.1.tar.bz2 --output slurm-23.02.1.tar.bz2
tar --bzip -x -f slurm-23.02.1.tar.bz2
cd slurm-23.02.1
source /etc/bashrc # need python3 loaded
spack unload -a
./configure --sysconfdir=/usr/local/etc --with-hdf5=no
make
make install
ldconfig -n /usr/local/lib
mkdir -p /var/spool/slurmd
# slurmctld # start slurm
# slurmd
# scontrol update nodename=$HOSTNAME state=DOWN reason="undraining"
# scontrol update nodename=$HOSTNAME state=RESUME