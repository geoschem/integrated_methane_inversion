# Running the IMI in a docker container 
This directory contains the relevant files for containerizing the IMI using docker. The IMI docker container is a prepackaged software environment that has all the dependencies and imi source code needed to run the IMI. This can be desirable to avoid the need to setup the software environment to run the IMI and more easily facilitates automation. The container should be portable to *all* (hopefully) systems regardless of architecture or operating system. 

## Prerequisites
### software
[Docker](https://www.docker.com/) must be installed to run or build IMI docker containers. This is the only dependency. That's the point of docker containers.

### Hardware
You must give the docker container enough resources to run GEOS-Chem simulations at the relevant resolution. See the [GEOS-Chem harware requirements](https://geos-chem.readthedocs.io/en/latest/getting-started/system-req-hard.html) for more information.

## Building and running the image
Important: Make sure you are in the top-level directory of the IMI source code.
```
$ docker build -f resources/containers/Dockerfile -t imi-docker-image . --platform=linux/amd64
```
## Running the container
In order to access the files from the inversion it is best to mount a volume from your local system onto the docker container. This allows the results of the inversion to persist after the container exits.

To do so create a volume with:
`docker volume create imi_output_dir`
Then run with the mounted volume:
```
$ docker run --platform=linux/amd64 --mount source=imi_output_dir,target=/home/al2/imi_output_dir imi-docker-image:latest
```
## pushing the image to remote repository
```
$ aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 753979222379.dkr.ecr.us-east-1.amazonaws.com
$ docker tag imi-docker-image:latest 753979222379.dkr.ecr.us-east-1.amazonaws.com/imi-docker-repository:latest
$ docker push 753979222379.dkr.ecr.us-east-1.amazonaws.com/imi-docker-repository:latest
```
## pulling the image
Note: this image is currently not publically available
`$ docker pull 753979222379.dkr.ecr.us-east-1.amazonaws.com/imi-docker-repository:latest`