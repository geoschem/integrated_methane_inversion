==============================
Using the IMI Docker container
==============================

What is a container?
====================
A Docker container is a lightweight, standalone, and executable software package that encapsulates an 
application and all its dependencies, including libraries, frameworks, and system tools. Docker containers 
provide a consistent and reproducible environment, ensuring that an application can run consistently 
across different systems, such as local clusters, cloud servers, and even local computers.


Why use the IMI Docker container?
=================================
Aside from providing a consistent environment, using the IMI container can significantly ease installation 
of the IMI on a new system. This is because the container has all the necessary dependencies and source code 
for running the IMI preinstalled and preloaded. This equals easier setup for you.

Additionally, Docker containers lend themselves very well to automated workflows, so using a docker 
container version of the IMI can make it easier to set up scheduled inversions of the IMI.


How to use the IMI Docker container
===================================

-------------
Prerequisites
-------------
To use the IMI Docker container, you must have Docker installed on your system. Docker can be installed on 
Windows, Mac, and Linux systems. For instructions on how to install Docker on your system, see the 
`Docker documentation <https://docs.docker.com/get-docker/>`__.

Additionally, configuring the container for your particular application is much easier if you have the 
docker compose plugin installed as well. For instructions on how to install docker compose, see the 
`Docker documentation <https://docs.docker.com/compose/install/>`__.

Note: if your cluster does not support Docker, you can also use Singularity as an alternative to Docker. See 
the section on `Using Singularity instead of Docker <#using-singularity-instead-of-docker>`__ for more information.

-----------------
Pulling the image
-----------------
To run the container you will first need to pull the image from our cloud repository. There are two flavors of the
IMI docker container using a base operating system of ubuntu or amazon linux 2

For Amazon Linux 2::

    $ docker pull public.ecr.aws/w1q7j9l2/imi-al2-docker-image:latest


For Ubuntu::

    $ docker pull public.ecr.aws/w1q7j9l2/imi-ubuntu-docker-image:latest

-------------------------------
Setting up the compose.yml file
-------------------------------

The IMI needs access to both input data and personalized configuration variables for running the inversion for 
your desired region and period of interest. In order to supply these settings we use a docker 
`compose.yml <https://docs.docker.com/compose/compose-file/03-compose-file/>`__ file. The compose file allows 
you to input environment variables and mount files/directories from your local system into the container. This 
allows you to more easily configure the IMI and save the output directory to your local system.

IMI input data
--------------
The IMI needs input data in order to run the inversion. If you do not have the necessary input data available 
locally then it will be automatically downloaded from AWS, where the input data is publicly available.

If you already have the necessary input data available locally, then you can mount it to the IMI container in the 
`volumes` section of the compose.yml file without setting your aws credentials. Eg:::

    volumes:
        - /local/input/data:/home/al2/ExtData # mount input data directory


Storing the output data
-----------------------
In order to access the files from the inversion it is best to mount a volume from your local system onto the docker 
container. This allows the results of the inversion to persist after the container exits. We recommend making a 
dedicated IMI output directory using `mkdir`.::

    # Note replace /home/al2 with /home/ubuntu if using ubuntu flavor
    volumes:
        - /local/output/dir/imi_output:/home/al2/imi_output_dir # mount output directory
        - /local/container/config.yml:/home/al2/integrated_methane_inversion/config.yml # mount desired config file

Updating the config.yml file
----------------------------

The config.yml file configures the IMI to run according to your specific inversion requirements. There are two 
mechanisms to update the config.yml file:

1. If you would only like to update specific variables you can pass them in as environment variables:

All environment variables matching the pattern ``IMI_<config-variable-name>`` will update their corresponding config.yml 
variable. For example:::

    environment:
        - IMI_StartDate=20200501 
        - IMI_EndDate=20200601

will replace the ``StartDate`` and ``EndDate`` in the IMI config.yml file.

2. Replace the entire config.yml file with one from the host system:

To apply a config.yml file from your local system to the docker container, specify it in your compose.yml file as a 
volume. Then set the ``IMI_CONFIG_PATH`` environment variable to point to that path. Eg:::

    # Note replace /home/al2 with /home/ubuntu if using ubuntu flavor
    volumes:
        - /local/path/to/config.yml:/home/al2/integrated_methane_inversion/config.yml # mount desired config file
    environment:
        - IMI_CONFIG_PATH=/home/al2/integrated_methane_inversion/config.yml # should point to the path in the container


Note: any env variables matching the pattern specified in option 1 will overwrite the corresponding config vars in 
`IMI_CONFIG_PATH`.

Example compose.yml file
------------------------
This is an example of what a fully filled out compose.yml file looks like:::

    # IMI Docker Compose File
    # This file is used to run the IMI Docker image
    # and define important parameters for the container
    services:
      imi:
        image: public.ecr.aws/w1q7j9l2/imi-docker-image:latest
        volumes:
        # comment out any volume mounts you do not need for your system
          - /local/container/config.yml:/home/al2/integrated_methane_inversion/config.yml # mount desired config file
          - /local/input/data:/home/al2/ExtData # mount input data directory
          - /local/output/dir/imi_output:/home/al2/imi_output_dir # mount output directory
        environment:
        # comment out any environment vars you do not need for your system
          - IMI_CONFIG_PATH=config.yml # path starts from /home/al2/integrated_methane_inversions


Running the IMI
---------------
Once you have configured the compose.yml file, you can run the IMI by running:::

    $ docker compose up


from the same directory as your ``compose.yml`` file. This will start the IMI container and run the inversion. 
The output will be saved to the directory you specified in the compose.yml file. 

Alternatively, if you chose not to install ``docker compose`` you should be able to run the IMI using the 
`docker run <https://docs.docker.com/engine/reference/commandline/run/>`__ command, but this requires specifying 
all env variables and volumes via flags.

Using Singularity instead of Docker
===================================
We use Docker `Docker <https://docs.docker.com/get-started/overview/>`__ to containerize the IMI, but it may also 
be possible to use other container engines, like `Singularity <https://docs.sylabs.io/guides/3.5/user-guide/introduction.html>`__. 
Singularity is a container engine designed to run on HPC systems and local clusters, as some clusters do not allow 
Docker to be installed.
Note: using Singularity to run the IMI is untested and may not work as expected.

First pull the image:::

    $ singularity pull public.ecr.aws/w1q7j9l2/imi-ubuntu-docker-image:latest

Then run the image:::

    $ singularity run imi-ubuntu-docker-image_latest.sif
