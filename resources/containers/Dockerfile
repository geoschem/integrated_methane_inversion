# ************************* Base Dockerfile *************************
# Installs the conda/spack environments and system dependencies
# Note: This will take a long time to build due to the spack 
# installations
FROM amazonlinux:2

RUN mkdir install-scripts
RUN mkdir /home/al2/

# install conda
ADD install-scripts/install-conda.sh /install-scripts/install-conda.sh
RUN source install-scripts/install-conda.sh

# install spack
ADD install-scripts/install-spack.sh /install-scripts/install-spack.sh
RUN source install-scripts/install-spack.sh

# install lustre
RUN uname -r && amazon-linux-extras install -y lustre

# install geoschem deps with spack
ADD install-scripts/install-spack-env.sh /install-scripts/install-spack-env.sh
ADD install-scripts/geoschem_deps-gnu-openmpi-102.yaml /install-scripts/geoschem_deps-gnu-openmpi-102.yaml
RUN source install-scripts/install-spack-env.sh

# install system deps with yum
ADD install-scripts/install-system-deps.sh /install-scripts/install-system-deps.sh
RUN source install-scripts/install-system-deps.sh

# install gcpy dependencies
ADD install-scripts/install-gcpy.sh /install-scripts/install-gcpy.sh
ADD install-scripts/gcpy_environment.yaml /install-scripts/gcpy_environment.yaml
RUN source install-scripts/install-gcpy.sh

WORKDIR /home/al2/

ADD entrypoint.sh /home/al2/entrypoint.sh
RUN chmod +x entrypoint.sh

CMD [ "entrypoint.sh" ]
