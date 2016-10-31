## storm-analysis ##
This is a repository of code developed in the [Zhuang Lab](http://zhuang.harvard.edu/) for analysis of STORM movies.

Some algorithms were developed in other groups and ported to Python. In this case the license applies only to our implementation of the code. If you plan to use the algorithm in a commercial product you should discuss this with the original developers.

The code is now located in the storm_analysis folder.

The code should work with both Python2.7 and Python3.

Settings file and splines for the 2016 SMLM Challenge are [here](http://zhuang.harvard.edu/smlm2016_settings.zip) (261MB).

## Installation ##

### Option 1 ###
Use the Anaconda Python distribution which makes the installation and dependencies management very easy : https://www.continuum.io/downloads

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
python setup.py install
```

Test the installation (this will take a few minutes to run):

```
python -c "import storm_analysis"
```

### Option 2 ###
This is straight forward on Linux but is likely difficult on Windows

Get the `storm_analysis` source code using git: 

```sh
git clone https://github.com/ZhuangLab/storm-analysis.git
cd storm-analysis
```

(Optional) Create an environment to keep your main Python installation clean:

```sh
virtualenv venv
source venv/bin/activate
```

Install dependencies and storm-analysis:
```sh
pip install numpy scipy matplotlib pillow tifffile shapely randomcolor
python setup.py install
```

Test the installation (this will take a few minutes to run):
```sh
cd storm_analysis/test
nose2
```

## General Notes ##
Questions should be addressed to Hazen Babcock (hbabcock _at_ fas.harvard.edu)
