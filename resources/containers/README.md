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

### IMI input data
The IMI needs input data in order to run the inversion. If you do not have the necessary input data available locally then you will need to give the IMI container access to S3 on AWS, where the input data is available. This can be done by passing your [aws credentials](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-envvars.html#envvars-set) as environment variables or in a .env file at runtime using `-e AWS_ACCESS_KEY_ID=<your-access-key> -e AWS_SECRET_ACCESS_KEY=<your-secret-key>`. Note that your AWS account will be charged for data download outside of the us-east-1 region (eg. to your local cluster).

If you already have the necessary input data available locally, then you can mount it to the IMI container via a docker volume using the `-v` flag. Eg. `-v /path/to/local/data:/home/al2/target/path`. Just make sure to set the corresponding config.yml variable to point at the target path as well.

### Storing the output data
In order to access the files from the inversion it is best to mount a volume from your local system onto the docker container using the `-v` flag. This allows the results of the inversion to persist after the container exits.

Additionally, in order to access the input data needed to run the IMI, you will need to give it access

To do so create a volume with:
`docker volume create imi_output_dir`
Then run with the mounted volume:
```
$ docker run --platform=linux/amd64 --mount source=imi_output_dir,target=/home/al2/imi_output_dir imi-docker-image:latest
```

To view the run directory you can find the path to your volume on the host system by running:
`docker volume inspect imi_output_dir` and looking in the corresponding `MountPoint` directory.
### Updating the config.yml file
There are two mechanisms to update the config.yml;

1. If you would only like to update specific variables you can pass them in as environment variables:

All environment variables matching the pattern `IMI_<config-variable-name>` will update their corresponding config.yml variable. For example
```
$ docker run --platform=linux/amd64 \
-e IMI_StartDate="20200501" \
-e IMI_EndDate="20200601" \
--mount source=imi_output_dir,target=/home/al2/imi_output_dir imi-docker-image:latest
```
```

$ docker run --platform=linux/amd64 \
--env-file .env \
-e IMI_CONFIG_CONTENTS="$(cat /Users/lucasestrada/Projects/IMI/integrated_methane_inversion/resources/containers/container_config.yml)" \
-v /Users/lucasestrada/Projects/IMI/imi_output:/home/al2/imi_output_dir imi-docker-image:latest
```
will replace the `StartDate` and `EndDate` in the IMI config.yml file.

2. Replace the entire config.yml file with one from the host system:

The environment variable `IMI_CONFIG_CONTENTS`, if present, will be interpreted as the config.yml file for your inversion. Eg:

```
$ docker run --platform=linux/amd64 \
-e IMI_CONFIG_CONTENTS="$(cat /path/to/config.yml)" \
--mount source=imi_output_dir,target=/home/al2/imi_output_dir imi-docker-image:latest
```
Note: any env variables matching the pattern specified in option 1 will overwrite the corresponding config vars in `IMI_CONFIG_CONTENTS`.

## pushing the image to remote repository
```
$ aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 753979222379.dkr.ecr.us-east-1.amazonaws.com
$ docker tag imi-docker-image:latest 753979222379.dkr.ecr.us-east-1.amazonaws.com/imi-docker-repository:latest
$ docker push 753979222379.dkr.ecr.us-east-1.amazonaws.com/imi-docker-repository:latest
```
## pulling the image
Note: this image is currently not publically available
`$ docker pull 753979222379.dkr.ecr.us-east-1.amazonaws.com/imi-docker-repository:latest`