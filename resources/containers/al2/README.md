# Running the IMI in a docker container 
This directory contains the relevant files for containerizing the IMI using docker. The IMI docker container is a prepackaged software environment that has all the dependencies and imi source code needed to run the IMI. This can be desirable to avoid the need to setup the software environment to run the IMI and more easily facilitates automation. The container should be portable to *all* (hopefully) systems regardless of architecture or operating system. 

## Prerequisites
### software
[Docker](https://www.docker.com/) must be installed to run or build IMI docker containers. This is the only dependency. That's the point of docker containers.
[Docker Compose Plugin](https://docs.docker.com/compose/install/) This allows use of the docker compose.yml file for configuring IMI runs, but is not strictly necessary if you choose to enter in all env variables and bind mounts manually via the [docker run](https://docs.docker.com/engine/reference/commandline/run/) command.

### Hardware
You must give the docker container enough resources to run GEOS-Chem simulations at the relevant resolution. See the [GEOS-Chem harware requirements](https://geos-chem.readthedocs.io/en/latest/getting-started/system-req-hard.html) for more information.

## pulling the image
To run the image you will first need to pull the image from our cloud repository
`$ docker pull public.ecr.aws/w1q7j9l2/imi-al2-docker-image:latest`

## Setting up the compose.yml file
The IMI needs access to both input data and personalized configuration variables for running the inversion for your desired region and period of interest. In order to supply these settings we use a docker [compose.yml](https://docs.docker.com/compose/compose-file/03-compose-file/) file. The compose file allows you to input environment variables and mount files/directories from your local system into the container.


### IMI input data
The IMI needs input data in order to run the inversion. If you do not have the necessary input data available locally then you will need to give the IMI container access to S3 on AWS, where the input data is available. This can be done by specifying your [aws credentials](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-envvars.html#envvars-set) in the `environment` section of the compose.yml file. Eg:

```
environment:
- AWS_ACCESS_KEY_ID=your_access_key_id
- AWS_SECRET_ACCESS_KEY=your_secret_access_key
- AWS_DEFAULT_REGION=us-east-1
```

Note: these credentials are sensitive, so do not post them publicly in any repository.

If you already have the necessary input data available locally, then you can mount it to the IMI container in the `volumes` section of the compose.yml file without setting your aws credentials. Eg:

```
volumes:
  - /local/input/data:/home/al2/ExtData # mount input data directory
```

### Storing the output data
In order to access the files from the inversion it is best to mount a volume from your local system onto the docker container. This allows the results of the inversion to persist after the container exits. We recommend making a dedicated IMI output directory using `mkdir`.

```
volumes:
  - /local/output/dir/imi_output:/home/al2/imi_output_dir # mount output directory
```
### Updating the config.yml file

There are two mechanisms to update the config.yml;

1. If you would only like to update specific variables you can pass them in as environment variables:

All environment variables matching the pattern `IMI_<config-variable-name>` will update their corresponding config.yml variable. For example:

```
environment:
  - IMI_StartDate=20200501 
  - IMI_EndDate=20200601
```
will replace the `StartDate` and `EndDate` in the IMI config.yml file.

2. Replace the entire config.yml file with one from the host system:

To apply a config.yml file from your local system to the docker container, specify it in your compose.yml file as a volume. Then set the `IMI_CONFIG_PATH` environment variable to point to that path. Eg:

```
volumes:
  - /local/path/to/config.yml:/home/al2/integrated_methane_inversion/config.yml # mount desired config file
environment:
  - IMI_CONFIG_PATH=/home/al2/integrated_methane_inversion/config.yml # should point to the path in the container
```

Note: any env variables matching the pattern specified in option 1 will overwrite the corresponding config vars in `IMI_CONFIG_PATH`.

## Running the IMI
Once you have configured the compose.yml file, you can run the IMI by running:
```
$ docker compose up
```
while in the same directory as the compose.yml file.

Alternatively, if you chose not to install `docker compose` you should be able to run the IMI using the [docker run](https://docs.docker.com/engine/reference/commandline/run/) command, but this requires specifying all env variables and volumes via flags.

## Developers info
Some users may wish to modify the IMI source code and build their own version of the IMI docker container. This section tells those brave souls how to rebuild the container once they have made their desired code changes. Note: If you introduce new software dependencies, you may need to update the IMI base-image container, which contains all the software dependencies for the IMI (see the base-image directory for more details).
### Building and running the image
Important: Make sure you are in the top-level directory of the IMI source code.
```
$ docker build -f resources/containers/al2/Dockerfile -t imi-al2-docker-image . --platform=linux/amd64
```
### Pushing the image to remote repository
Update these as you see fit for your desired aws or docker repository
```
$ aws ecr-public get-login-password --region us-east-1 | docker login --username AWS --password-stdin public.ecr.aws/w1q7j9l2
$ docker tag imi-docker-image:latest public.ecr.aws/w1q7j9l2/imi-al2-docker-image:latest
$ docker push public.ecr.aws/w1q7j9l2/imi-al2-docker-image:latest
```
