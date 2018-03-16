# Distro
FROM jupyter/scipy-notebook

LABEL maintainer="Hazen Babcock <hbabcock@mac.com>"

USER root

RUN apt update
RUN apt-get --yes install gcc
RUN apt-get --yes install git
RUN apt-get --yes install libfftw3-dev
RUN apt-get --yes install libgeos-dev
RUN apt-get --yes install liblapack-dev
RUN apt-get --yes install scons

USER $NB_UID

RUN pip install pillow
RUN pip install tifffile
RUN pip install shapely
RUN pip install pytest
RUN pip install astropy
RUN pip install randomcolor

# Get current storm-analysis
RUN mkdir code
WORKDIR code
RUN git clone https://github.com/ZhuangLab/storm-analysis.git

# Set working directory
WORKDIR storm-analysis

# run scons to build the C libraries.
RUN scons

# install storm-analysis.
RUN /bin/bash -c "python setup.py install;"

# Record when this image was made
RUN date > image_date.txt
RUN git log -1 --stat > sa_version.txt

RUN fix-permissions $CONDA_DIR
RUN fix-permissions /home/$NB_USER
    
USER $NB_UID


