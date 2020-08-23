## storm-analysis ##
This is a repository of code developed in the [Zhuang lab](http://zhuang.harvard.edu/) for analysis of STORM movies.

Some algorithms were developed in other groups and ported to Python. In this case the license applies only to our implementation of the code. If you plan to use the algorithm in a commercial product you should discuss this with the original developers ([IANAL](https://en.wikipedia.org/wiki/IANAL)).

The code should work with both Python2.7 and Python3, but Python3 is recommended. Python2 support will be deprecated in January 2019.

[![Linux Build Status](https://travis-ci.org/ZhuangLab/storm-analysis.svg?branch=master)](https://travis-ci.org/ZhuangLab/storm-analysis)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/nr6aha5hsn2g84j1?svg=true)](https://ci.appveyor.com/project/HazenBabcock/storm-analysis)
[![Documentation Status](https://readthedocs.org/projects/storm-analysis/badge/?version=latest)](https://readthedocs.org/projects/storm-analysis/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3528330.svg)](https://doi.org/10.5281/zenodo.3528330)

[Discussion Group](https://groups.google.com/d/forum/storm-analysis)

## Documentation ##

The documentation for the latest release is [here](http://storm-analysis.readthedocs.io/en/stable/).

The documentation for git head is [here](http://storm-analysis.readthedocs.io/en/latest/).


## Jupyter Notebooks ##

[Jupyter](http://jupyter.org/) notebooks that document how to use this project are available in the `jupyter_notebooks` directory.

Additional notebooks are available [here](https://drive.google.com/drive/folders/1k5vkzisz_I3XwXIw-2G1iOJLe996y_Wu).

Thanks to [Google Colaboratory](https://colab.research.google.com/notebooks/welcome.ipynb) you can install and run storm-analysis on a virtual machine in Google's cloud, no local install necessary. Please see the notebooks in the `colab` directory [here](https://drive.google.com/drive/folders/1k5vkzisz_I3XwXIw-2G1iOJLe996y_Wu). You will need a Google account to do this.


## Docker ##

[Docker](https://www.docker.com) images are available [here](https://hub.docker.com/u/zhuanglab/).

## Dependencies ##

### C ###

* [FFTW3](http://www.fftw.org/)
* [LAPACK](http://www.netlib.org/lapack/)

### Python ###

* [numpy](http://www.numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)
* [pillow](https://python-pillow.org/)
* [tifffile](https://pypi.python.org/pypi/tifffile)
* [Shapely](https://pypi.python.org/pypi/Shapely)
* [randomcolor](https://pypi.python.org/pypi/randomcolor)
* [PyWavelets](https://pypi.python.org/pypi/PyWavelets)
* [PyQt5](https://pypi.python.org/pypi/PyQt5)
* [h5py](https://www.h5py.org/)
* [astropy](http://www.astropy.org/)
