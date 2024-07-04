# ************************* Operational Dockerfile *************************
# Starts from base image and adds automation scripts to it
# Dockerfile for base image is located in base-image directory
FROM public.ecr.aws/w1q7j9l2/imi-ubuntu-base-image:latest

USER ubuntu

# make relevant directories for the imi
RUN mkdir /home/ubuntu/imi_output_dir && mkdir /home/ubuntu/ExtData

# add imi source code to base image
ADD ./ /home/ubuntu/integrated_methane_inversion

# change ownership of the imi source code to ubuntu user
RUN sudo chown ubuntu:ubuntu -R /home/ubuntu/integrated_methane_inversion/

# make imi run scripts executable
RUN chmod -R +x /home/ubuntu/integrated_methane_inversion/

# update default config.yml file TODO remove
RUN sed -i "s/v2024-03/v2023-06/g" ~/integrated_methane_inversion/config.yml
RUN sed -i "s|CondaFile:.*|CondaFile: /home/ubuntu/.bashrc|g" ~/integrated_methane_inversion/config.yml
RUN sed -i "s|CondaEnv:.*|CondaEnv: imi_env|g" ~/integrated_methane_inversion/config.yml

CMD [ "/home/ubuntu/integrated_methane_inversion/resources/containers/ubuntu/entrypoint.sh" ]
