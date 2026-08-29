## storm-analysis ##
This is a repository of code developed in the [Zhuang lab](http://zhuang.harvard.edu/) for analysis of STORM movies.

Some algorithms were developed in other groups and ported to Python. In this case the license applies only to our implementation of the code. If you plan to use the algorithm in a commercial product you should discuss this with the original developers ([IANAL](https://en.wikipedia.org/wiki/IANAL)).

The code has most recently been tested with Python3.10.

[![Tests](https://github.com/ZhuangLab/storm-analysis/actions/workflows/python-package.yml/badge.svg)](https://github.com/ZhuangLab/storm-analysis/actions/workflows/python-package.yml)
[![Documentation Status](https://readthedocs.org/projects/storm-analysis/badge/?version=latest)](https://readthedocs.org/projects/storm-analysis/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3528330.svg)](https://doi.org/10.5281/zenodo.3528330)


## Documentation ##

The documentation for the latest release is [here](http://storm-analysis.readthedocs.io/en/stable/).

The documentation for git head is [here](http://storm-analysis.readthedocs.io/en/latest/).


## Installation ##

Installation instructions are [here](http://storm-analysis.readthedocs.io/en/stable/install.html).

On Windows you can skip building the C libraries by installing the wheel attached to the [latest release](https://github.com/ZhuangLab/storm-analysis/releases/latest). It contains the compiled DLLs, so no C compiler is needed:

```
pip install storm_analysis-<version>-py3-none-win_amd64.whl
```


## Jupyter Notebooks ##

[Jupyter](http://jupyter.org/) notebooks that document how to use this project are available in the `jupyter_notebooks` directory.

Additional notebooks are available [here](https://drive.google.com/drive/folders/1k5vkzisz_I3XwXIw-2G1iOJLe996y_Wu).


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
