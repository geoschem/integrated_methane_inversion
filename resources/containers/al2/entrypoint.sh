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

spack env activate compute_env
micromamba activate imi_env

# set up environment variables for GEOS-Chem
export NETCDF_HOME=$(spack location -i netcdf-c)
export NETCDF_FORTRAN_HOME=$(spack location -i netcdf-fortran)
export LD_LIBRARY_PATH=${NETCDF_HOME}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NETCDF_FORTRAN_HOME}/lib

export FC=gfortran
export CC=gcc
export CXX=g++

# Tell GEOS-Chem where to find netCDF library files
export GC_BIN=$NETCDF_HOME/bin
export GC_INCLUDE=$NETCDF_HOME/include
export GC_LIB=$NETCDF_HOME/lib

# NOTE: If netCDF-Fortran was loaded as a separate module, then
# also define these variables.  (Otherwise comment these out.)
export GC_F_BIN=$NETCDF_FORTRAN_HOME/bin
export GC_F_INCLUDE=$NETCDF_FORTRAN_HOME/include
export GC_F_LIB=$NETCDF_FORTRAN_HOME/lib

# Max out the stack memory for OpenMP
# http://wiki.seas.harvard.edu/geos-chem/index.php/GNU_Fortran_compiler#Requesting_sufficient_stack_memory_for_GEOS-Chem
# operation not permitted
# ulimit -s unlimited
export OMP_STACKSIZE=500m

cd /home/al2/integrated_methane_inversion

if [[ "x${IMI_CONFIG_PATH}" != "x" ]]; then
	config_file=$IMI_CONFIG_PATH
else
	config_file="config.yml"
	# remove default config file and replace with docker config file
	cp resources/containers/container_config.yml config.yml
fi

# override specific config file vars with env variables of 
# the syntax IMI_<config-variable>
chmod +x src/utilities/override_config_variables.py
python src/utilities/override_config_variables.py $config_file $config_file

# sbatch -W --mem 2000 -c 1 run_imi.sh resources/containers/container_config.yml; wait;
./run_imi.sh $config_file | tee imi_output.log

# Fix issue when switching instance types where node claims to be drained
# scontrol update nodename=$HOSTNAME state=idle
# scontrol update nodename=$HOSTNAME state=DOWN reason="undraining"
# scontrol update nodename=$HOSTNAME state=RESUME

# keep container running for testing purposes
while :
do
	echo "IMI Finished. Container running forever. Shut down manually to stop"
	sleep 10
done