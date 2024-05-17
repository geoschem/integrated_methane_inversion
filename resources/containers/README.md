# IMI containers
This directory contains the files and resources for building the IMI containers. The containers are built using the Dockerfile in the subdirectories of this directory. We
build two types of containers: one using an ubuntu base image and a second using an
amazon linux 2 base image. The containers are built on each commit to the dev branch 
and with tagged versions. The containers are pushed to the IMI ECR public repository.

## Pulling the containers
The containers are available in the IMI ECR public repository. To pull the containers, you will need to have docker installed. You can pull the containers using the following commands:

For ubuntu operational image:
```bash
docker pull public.ecr.aws/w1q7j9l2/imi-ubuntu-docker-image:latest
```

For amazon linux 2 operational image:
```bash 
docker pull public.ecr.aws/w1q7j9l2/imi-al2-docker-image:latest
```
