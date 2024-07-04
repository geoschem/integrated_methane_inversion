# ************************* Base Dockerfile *************************
# Installs the micromamba/spack environments and system dependencies
# Note: This will take a long time to build due to the spack 
# installations
# Use Ubuntu 24.04 as base image
FROM ubuntu:24.04

# Update package list and install dependencies
RUN apt-get update && \
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
munge \
slurm-wlm \
sudo \
python3 \
vim \
libpmix2 \
curl \
tar \
bzip2 \
zip \
ca-certificates \
git \
build-essential \
coreutils \ 
environment-modules \
gfortran \
gpg \
lsb-release \
python3-venv \
unzip \
libssl-dev \
libcurl4-openssl-dev \
&& rm -rf /var/lib/apt/lists/*

# Copy sudoers file to grant passwordless sudo
COPY install-scripts/slurm/sudoers /etc/sudoers.d/

# Switch to ubuntu user
USER ubuntu

# Add slurm install scripts
ADD install-scripts/slurm/base_slurm.conf /home/ubuntu/install-scripts/base_slurm.conf
ADD install-scripts/slurm/cgroup.conf /home/ubuntu/install-scripts/cgroup.conf
ADD install-scripts/slurm/configure_slurm.py /home/ubuntu/install-scripts/configure_slurm.py
ADD install-scripts/slurm/test_slurm.sh /home/ubuntu/install-scripts/test_slurm.sh
ADD install-scripts/slurm/install-slurm.sh /home/ubuntu/install-scripts/install-slurm.sh

# Install slurm
RUN sudo chmod +x /home/ubuntu/install-scripts/install-slurm.sh
RUN /home/ubuntu/install-scripts/install-slurm.sh

# Install Spack
ADD install-scripts/install-spack.sh /home/ubuntu/install-scripts/install-spack.sh
RUN sudo chmod +x /home/ubuntu/install-scripts/install-spack.sh
RUN /home/ubuntu/install-scripts/install-spack.sh

# Install Spack env
ADD install-scripts/install-spack-env.sh /home/ubuntu/install-scripts/install-spack-env.sh
ADD install-scripts/geoschem_deps-gnu-openmpi-122.yml /home/ubuntu/install-scripts/geoschem_deps-gnu-openmpi-122.yml
RUN sudo chmod +x /home/ubuntu/install-scripts/install-spack-env.sh
RUN /home/ubuntu/install-scripts/install-spack-env.sh

# Install micromamba
ADD install-scripts/install-micromamba.sh /home/ubuntu/install-scripts/install-micromamba.sh
ADD install-scripts/imi_env.yml /home/ubuntu/install-scripts/imi_env.yml
RUN sudo chmod +x /home/ubuntu/install-scripts/install-micromamba.sh
RUN /home/ubuntu/install-scripts/install-micromamba.sh

# Install aws cli
RUN cd ~ \
&& curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
&& unzip awscliv2.zip \
&& sudo ./aws/install
