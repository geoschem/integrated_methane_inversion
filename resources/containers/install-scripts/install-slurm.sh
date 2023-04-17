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
export MUNGEUSER=991
groupadd -g $MUNGEUSER munge
useradd  -m -c "MUNGE Uid 'N' Gid Emporium" -d /var/lib/munge -u $MUNGEUSER -g munge  -s /sbin/nologin munge
export SLURMUSER=992
groupadd -g $SLURMUSER slurm
useradd  -m -c "SLURM workload manager" -d /var/lib/slurm -u $SLURMUSER -g slurm  -s /bin/bash slurm

# install munge
amazon-linux-extras install epel -y
yum install munge munge-libs munge-devel -y
yum install rng-tools -y
rngd -r /dev/urandom
/usr/sbin/create-munge-key -r
dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key
chown munge: /etc/munge/munge.key
chmod 400 /etc/munge/munge.key
sudo chown -R munge:munge /var/log/munge
sudo chown root:root /var/log/munge
sudo chown root:root /var/log/munge/munged.log
munged -f


# install slurm
# set -e 
# curl https://download.schedmd.com/slurm/slurm-23.02.1.tar.bz2 --output slurm-23.02.1.tar.bz2
# tar --bzip -x -f slurm-23.02.1.tar.bz2
# cd slurm-23.02.1
# ./configure --sysconfdir=/usr/local/etc --with-hdf5=$(spack location -i netcdf-c)
# make
# make install
# ldconfig -n /usr/local/lib
# slurmctld -Dvvv # start slurm
sudo mkdir /run/slurm
sudo chown slurm:slurm /run/slurm
sudo mkdir /var/spool/slurmd
sudo chown slurm:slurm /var/spool/slurmd