# Distro
FROM ubuntu:bionic
MAINTAINER Hazen Babcock <hbabcock@mac.com>

# Update sources.
RUN apt update

# Get tzdate non-interactive.
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata

# Get dependencies.
RUN apt-get --yes install gcc
RUN apt-get --yes install git
RUN apt-get --yes install libfftw3-dev
RUN apt-get --yes install libgeos-dev
RUN apt-get --yes install liblapack-dev
RUN apt-get --yes install scons
RUN apt-get --yes install python
RUN apt-get --yes install python3
RUN apt-get --yes install python3-dev
RUN apt-get --yes install python3-pip
RUN apt-get --yes install python3-setuptools
RUN apt-get --yes install python3-tk
RUN apt-get --yes install zlib1g

# Configure python3 virtual environment with 3D-DAOSTORM dependencies.
RUN pip3 install virtualenv
RUN virtualenv -p /usr/bin/python3 venv
RUN /bin/bash -c "source venv/bin/activate; pip3 install numpy;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install scipy;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install matplotlib;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install pillow;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install tifffile;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install shapely;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install randomcolor;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install pywavelets;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install h5py;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install pytest;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install astropy;"

# Get current storm-analysis
RUN git clone https://github.com/ZhuangLab/storm-analysis.git

# Set working directory
WORKDIR /storm-analysis

# run scons to build the C libraries.
RUN scons

# install storm-analysis.
RUN /bin/bash -c "source /venv/bin/activate; python setup.py install;"

# Record when this image was made
RUN date > image_date.txt
RUN git log -1 --stat > sa_version.txt

# Volume for accessing STORM data
VOLUME ["/data"]

