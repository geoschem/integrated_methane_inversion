# IMI Base Environment Image
This repository contains the Dockerfile for the IMI Base Ubuntu Docker image. It is intended to be a template for testing different versions of the IMI on ubuntu. This will allow reproducible builds of AMI's on AWS. The Dockerfile is based on the official Ubuntu image and installs the necessary dependencies for the IMI. The build time for this image takes a long time to build the spack dependencies, so it is meant to contain the software dependencies that take significant time to build and are not regularly updated (gcc, esmf, conda, etc). The image is built and pushed to the public ECR repository for use in the IMI. This is then used as a base-image for the imi Dockerfile, which adds the relevant imi source code to the image.

# Dockerfile
The Dockerfiles contain a list of instructions to build an image with the necessary dependencies preinstalled. The Dockerfile starts with the ubuntu:24.04 image available on dockerhub. The Dockerfile then runs several installation scripts (located in the install-scripts/ directory) to install micromamba (conda), spack, slurm, and any other necessary dependencies

## Prerequisites
[Docker](https://www.docker.com/) must be installed to run or build IMI docker containers.

## Building the Docker image 
To build the docker image, cd into the resources/containers/base-image directory and run
`$ docker build -t imi-ubuntu-base-image . --platform=linux/amd64`
The `--platform=linux/amd64` flag specifies the architecture to mimic while building and ensures consistent builds regardless of whether the build is done on an arm system (apple silicon) or amd system (intel).

Building all of the spack dependencies can take hours. Grab a cup of coffee and let it run in the background.
When finished you can tag it and push it to the relevant remote repository (example is to aws):
```
# login to ecr repository
$ aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/w1q7j9l2
$ docker tag imi-ubuntu-base-image:latest public.ecr.aws/w1q7j9l2/imi-ubuntu-base-image:latest
$ docker push public.ecr.aws/w1q7j9l2/imi-ubuntu-base-image:latest
```
## Running the image
This runs the built image using the default `entrypoint.sh` command:

`docker run --platform=linux/amd64 imi-ubuntu-base-repository:latest`

## Pulling an existing image
A fully built base image is stored on aws and can be accessed by running:
```
$ docker pull public.ecr.aws/w1q7j9l2/imi-ubuntu-base-image:latest
```
# Conda environments
The python dependencies are built using the default imi conda imi_env.yml file.

# Spack environments
The spack dependencies are built using spack environments, which package a set of dependencies together similar to conda. Spack can be very finicky -- edit at your own risk.

# Slurm
The IMI uses the slurm scheduler to manage jobs and is specifically useful for allocating resources amongst the jacobian simulations. However, slurm and docker compatability is complicated. Slurm requires the processes of munge, slurmd, and slurmctld to be run in the background to function. Best practices with docker is to setup a separated docker container for each process and share filesystems and pass messages as needed between containers. For simplicity, and to maximize use of available hardware resources, we run everything in a single container. The current slurm installation ~works~ but it is possible that the configuration could be improved.

## Modifying version numbers for spack dependencies
The spack dependency list for installation are contained in the install-scripts/geoschem_deps-gnu-openmpi-122.yml file. You can modify the version numbers in this configuration file and rebuild the docker container to install the modified versions.


