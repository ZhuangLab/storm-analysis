#!/bin/bash
#
# This file is based on the following example at conda.io:
#
# http://conda.pydata.org/docs/travis.html#using-conda-with-travis-ci
#

# Update apt and install lapack and libfftw3
sudo apt update -qq
sudo apt-get --yes install liblapack-dev
sudo apt-get --yes install libfftw3-dev

# Download Python2 or Python3 Miniconda.
if [[ $TOXENV == "2.7" ]]
then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
fi
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

# What does this do?
hash -r

# Not interactive, so always yes.
conda config --set always_yes yes --set changeps1 no

# Update conda?
conda update -q conda

# Useful for debugging any issues with conda
conda info -a

# Create the test environment including gcc compiler.
conda create -q -n test-environment python=$TOXENV gcc

# Need Python2 for SCons.
if [[ $TOXENV == "2.7" ]];
then
    scons
else
    conda create -q -n py2-env python=2.7
    source activate py2-env
    scons
fi

# Activate virtual environment test-environment.
source activate test-environment

# Install conda dependencies.
conda config --add channels conda-forge
conda install numpy scipy pytest matplotlib
conda install tifffile pillow
conda install shapely randomcolor pywavelets

# Install nose2 testing package.
pip install nose2

# Install the storm-analysis project.
python setup.py install
