# Distro
FROM debian:latest
MAINTAINER Hazen Babcock <hbabcock@mac.com>

# Update sources
RUN apt update

# Get dependencies
RUN apt --yes install gcc
RUN apt --yes install git
RUN apt --yes install libfftw3-dev
RUN apt --yes install liblapack-dev
RUN apt --yes install scons
RUN apt --yes install python
RUN apt --yes install python3.4
RUN apt --yes install python3-dev
RUN apt --yes install python3-pip
RUN apt --yes install python3-setuptools
RUN apt --yes install zlib1g

# Configure python3 virtual environment with 3D-DAOSTORM dependencies.
RUN pip3 install virtualenv
RUN virtualenv -p /usr/bin/python3 venv
RUN /bin/bash -c "source venv/bin/activate; pip3 install numpy;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install pillow;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install scipy;"
RUN /bin/bash -c "source venv/bin/activate; pip3 install tifffile;"

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

# Volume for accessing STORM data
VOLUME ["/data"]

# Copy in the analysis shell script
COPY ./daostorm_3d.sh ./

# Run 3D-DAOSTORM analysis.
CMD ["/bin/bash", "daostorm_3d.sh"]
