## storm-analysis ##
This is a repository of code developed in the [Zhuang lab](http://zhuang.harvard.edu/) and the  [Babcock lab](https://hazenbabcock.github.io/) for analysis of STORM movies.

Some algorithms were developed in other groups and ported to Python. In this case the license applies only to our implementation of the code. If you plan to use the algorithm in a commercial product you should discuss this with the original developers.

The code should work with both Python2.7 and Python3, but Python3 is recommended.

Settings file and splines for the 2016 SMLM Challenge are [here](http://zhuang.harvard.edu/smlm2016_settings.zip) (261MB). These settings files are only compatible with versions up to the 1.0 release.

[![Linux Build Status](https://travis-ci.org/ZhuangLab/storm-analysis.svg?branch=master)](https://travis-ci.org/ZhuangLab/storm-analysis)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/nr6aha5hsn2g84j1?svg=true)](https://ci.appveyor.com/project/HazenBabcock/storm-analysis)
[![Documentation Status](https://readthedocs.org/projects/storm-analysis/badge/?version=latest)](https://readthedocs.org/projects/storm-analysis/)

[Discussion Group](https://groups.google.com/d/forum/storm-analysis)

## Documentation ##

The documentation for the latest release is [here](http://storm-analysis.readthedocs.io/en/stable/).

The documentation for git head is [here](http://storm-analysis.readthedocs.io/en/latest/).


## Jupyter Notebooks ##

[Jupyter](http://jupyter.org/) notebooks that document how to use the code are available in the `jupyter_notebooks` directory.


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

## General Notes ##
Questions should be addressed to Hazen Babcock (hbabcock _at_ fas.harvard.edu)
