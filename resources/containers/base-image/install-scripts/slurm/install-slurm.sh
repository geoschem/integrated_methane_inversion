#!/bin/bash

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