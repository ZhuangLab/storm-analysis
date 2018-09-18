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

RUN pip install --upgrade pip
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

# copy storm-analysis jupyter notebooks to the work directory.
RUN mkdir ../../work/sa_notebooks
RUN cp jupyter_notebooks/* ../../work/sa_notebooks/.

# create local mount point
RUN mkdir ../../work/share

# record when this image was made
RUN date > ../../work/image_date.txt
RUN git rev-parse HEAD > ../../work/sa_version.txt

RUN fix-permissions $CONDA_DIR
RUN fix-permissions /home/$NB_USER

# start in home directory.
WORKDIR /home/$NB_USER/work

USER $NB_UID


