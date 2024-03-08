# ************************* Base Dockerfile *************************
# Installs the micromamba/spack environments and system dependencies
# Note: This will take a long time to build due to the spack 
# installations
FROM amazonlinux:2

RUN mkdir /home/al2/
RUN mkdir /home/al2/install-scripts 
WORKDIR /home/al2/

# install micromamba
ADD install-scripts/install-micromamba.sh /home/al2/install-scripts/install-micromamba.sh
RUN source install-scripts/install-micromamba.sh

# install spack
ADD install-scripts/install-spack.sh /home/al2/install-scripts/install-spack.sh
RUN source install-scripts/install-spack.sh

# install lustre
RUN uname -r && amazon-linux-extras install -y lustre epel

# install geoschem deps with spack
ADD install-scripts/install-spack-env.sh /home/al2/install-scripts/install-spack-env.sh
ADD install-scripts/geoschem_deps-gnu-openmpi-102.yaml /home/al2/install-scripts/geoschem_deps-gnu-openmpi-102.yaml
RUN source install-scripts/install-spack-env.sh

# install system deps with yum
ADD install-scripts/install-system-deps.sh /home/al2/install-scripts/install-system-deps.sh
RUN source install-scripts/install-system-deps.sh

# install python dependencies
ADD install-scripts/install-imi-env.sh /home/al2/install-scripts/install-imi-env.sh
ADD install-scripts/imi_env.yml /home/al2/install-scripts/imi_env.yml
RUN source install-scripts/install-imi-env.sh

# create user and set to default user
RUN useradd -ms /bin/bash al2
# USER al2
ADD install-scripts/slurm/install-slurm.sh /home/al2/install-scripts/install-slurm.sh
RUN source install-scripts/install-slurm.sh

# slurm installation
ADD install-scripts/slurm/base_slurm.conf /home/al2/install-scripts/base_slurm.conf
ADD install-scripts/slurm/cgroup.conf /home/al2/install-scripts/cgroup.conf
ADD install-scripts/slurm/configure_slurm.py /home/al2/install-scripts/configure_slurm.py
