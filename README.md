## storm-analysis ##
This is a repository of code developed in the [Zhuang Lab](http://zhuang.harvard.edu/) for analysis of STORM movies.

Some algorithms were developed in other groups and ported to Python. In this case the license applies only to our implementation of the code. If you plan to use the algorithm in a commercial product you should discuss this with the original developers.

The code is now located in the storm_analysis folder.

The code should work with both Python2.7 and Python3.

Settings file and splines for the 2016 SMLM Challenge are [here](http://zhuang.harvard.edu/smlm2016_settings.zip) (261MB).

## Dependencies ##

### C ###

* [fftw3](http://www.fftw.org/)
* [lapack](http://www.netlib.org/lapack/)

### Python ###

* [numpy](http://www.numpy.org/)
* [scipy](https://www.scipy.org/)
* [matplotlib](http://matplotlib.org/)
* [pillow](https://python-pillow.org/)
* [tifffile](https://pypi.python.org/pypi/tifffile)
* [shapely](https://pypi.python.org/pypi/Shapely)
* [randomcolor](https://pypi.python.org/pypi/randomcolor)
* [pywavelets](https://pypi.python.org/pypi/PyWavelets)


## Installation ##

The use of [virtualenv](https://virtualenv.pypa.io/en/stable/) is highly recommended.

### Using wheels: ###

Wheels for 64 bit Windows are [here](http://zhuang.harvard.edu/storm_analysis/).

Install the above Python dependencies first. Christoph Gohlke is one good source of [wheels](http://www.lfd.uci.edu/~gohlke/pythonlibs/) for Windows.

```sh
pip install storm_analysis-XX.whl
```

### Using setup.py: ###

```sh
git clone https://github.com/ZhuangLab/storm-analysis.git
cd storm-analysis
python setup.py build_c
python setup.py install
```

This uses the [SCons](http://scons.org/) to build the C libraries. Note that SCons does not work with Python3 so you'll need to Python2 to build the C libraries.

The easiest way to install scons (be sure to use Python2) is:
```sh
pip install scons
```

#### Windows notes ####

An example build_c command using mingw64 from [nuwen](https://nuwen.net/mingw.html).

```sh
python setup.py build_c --scons-exe=C:\Python27\Scripts\scons.bat --compiler=mingw
```

The mingw64 gcc compiler must be in your path for this to work.

64 bit Windows LAPACK and FFTW3 DLLs are included in the project so you will not need to install them yourself.

### Using Anaconda: ###

Use the Anaconda Python distribution which makes installation and dependencies management very easy : https://www.continuum.io/downloads.
This only sort of works as the fftw package has not yet made it into Anaconda, but will hopefully be there soon.

(Optional) create an environment to keep your main Python installation clean : 

```sh
conda create -n my_env python
source activate my_env  # or activate my_env under Windows
```

Install dependencies : 

```sh
# Linux / OSX
conda config --add channels conda-forge 
conda install numpy fftw lapack pytest pytest-runner gcc tifffile matplotlib pillow shapely randomcolor pywavelets

# Windows
conda config --add channels conda-forge 
conda install numpy fftw lapack pytest pytest-runner m2w64-toolchain tifffile matplotlib pillow shapely randomcolor pywavelets
```

Get the `storm_analysis` source code using git:

```sh
git clone https://github.com/ZhuangLab/storm-analysis.git
cd storm-analysis
python setup.py build_c
python setup.py install
```

Test the installation (this will take a few minutes to run):

```
python -c "import storm_analysis"
```

## Testing ##

Test the installation (this will take a few minutes to run):
```sh
cd storm_analysis/test
nose2
```

Note: Due to issues with creating pickle files that are compatible across multiple OSs and versions of Python some of the tests may fail.

## General Notes ##
Questions should be addressed to Hazen Babcock (hbabcock _at_ fas.harvard.edu)
