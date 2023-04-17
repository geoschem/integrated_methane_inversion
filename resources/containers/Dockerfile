# ************************* Base Dockerfile *************************
# Installs the conda/spack environments and system dependencies
# Note: This will take a long time to build due to the spack 
# installations
FROM amazonlinux:2

RUN mkdir /home/al2/
RUN mkdir /home/al2/install-scripts 
WORKDIR /home/al2/

# install conda
ADD install-scripts/install-conda.sh /home/al2/install-scripts/install-conda.sh
RUN source install-scripts/install-conda.sh

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

ADD entrypoint.sh /home/al2/entrypoint.sh
RUN chmod +x entrypoint.sh

USER imi_user

CMD [ "entrypoint.sh" ]
