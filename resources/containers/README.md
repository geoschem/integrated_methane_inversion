# Environments 
This directory contains the files necessary to build an image with all the necessary dependencies for running the IMI. The build time for this image takes too long to run on github actions, so it is meant to contain the software dependencies that take significant time to build and are not regularly updated (gcc, esmf, conda, etc). This is then used as a base-image for the benchmarks/Dockerfile, which adds all of the scripts automating the GCClassic and GCHP runs.

# Dockerfile
Dockerfiles contain a list of instructions to build an image with the necessary dependencies preinstalled. The Dockerfile starts with the amazonlinux:2 image available on dockerhub -- this image already contains the aws-cli. The Dockerfile then runs several installation scripts (located in the install-scripts/ directory) to install conda, spack, and any other necessary dependencies

## Building the Docker image 
To build the docker image, cd into the benchmarks/environments/ directory and run
`$ aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 753979222379.dkr.ecr.us-east-1.amazonaws.com` # login to ecr
`$ docker build -t imi-base-repository . --platform=linux/amd64`
If you are running docker on a machine with an arm64 processor (M1 mac) you may need to append the following flag:
`--platform=linux/amd64`

Building all of the spack dependencies can take hours. Get a coffee and let it run in the background.
Now tag it and push:
```
$ docker tag geoschem-base-repository:latest 753979222379.dkr.ecr.us-east-1.amazonaws.com/geoschem-base-repository:latest
$ docker push 753979222379.dkr.ecr.us-east-1.amazonaws.com/geoschem-base-repository:latest
```

# Conda environments
The gcpy dependencies are built using a conda environment.yaml file. I originally attempted to install just using `conda install geoschem-gcpy --only-deps`, but xesmf doesn't work because it is not compatible with v8.4.0 of esmpy. This works for now. And we should at some point be able to install using only the above command.

# Spack environments
The spack dependencies are built using spack environments, which package a set of dependencies together similar to conda.

## Modifying version numbers for spack dependencies
The spack dependency list for installation are contained in the install-scripts/geoschem_deps-gnu-openmpi-102.yaml file. You can modify the version numbers in this configuration file and rebuild the docker container to install the modified versions.


