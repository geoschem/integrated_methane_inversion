# ************************* Operational Dockerfile *************************
# Starts from base image and adds automation scripts to it
# Dockerfile for base image is located in base-image directory
FROM public.ecr.aws/w1q7j9l2/imi-al2-base-image:latest

# # add spack environment
# ENV SPACK_ROOT="/opt/spack"

# make relevant directories for the imi
RUN mkdir /home/al2/imi_output_dir && mkdir /home/al2/ExtData

# add imi source code to base image
ADD ./ /home/al2/integrated_methane_inversion

# make imi run scripts executable
RUN chmod -R +x /home/al2/integrated_methane_inversion/

CMD [ "/home/al2/integrated_methane_inversion/resources/containers/entrypoint.sh" ]
