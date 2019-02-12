Installation
============

Dependencies
------------

C libraries
~~~~~~~~~~~

* `FFTW3 <http://www.fftw.org/>`_
* `LAPACK <http://www.netlib.org/lapack/>`_

Python
~~~~~~

* `numpy <http://www.numpy.org/>`_
* `scipy <https://www.scipy.org/>`_
* `matplotlib <http://matplotlib.org/>`_
* `pillow <https://python-pillow.org/>`_
* `tifffile <https://pypi.python.org/pypi/tifffile>`_
* `Shapely <https://pypi.python.org/pypi/Shapely>`_
* `randomcolor <https://pypi.python.org/pypi/randomcolor>`_
* `PyWavelets <https://pypi.python.org/pypi/PyWavelets>`_
* `PyQt5 <https://pypi.python.org/pypi/PyQt5>`_
* `h5py <http://www.h5py.org/>`_
* `astropy <http://www.astropy.org/>`_

Installing using wheels
-----------------------

Wheels for 64 bit Windows are `here <https://github.com/ZhuangLab/storm-analysis/releases>`_.
Despite their names, these will not work with 32 bit Python as the C libraries are 64bit.

Wheels for the latest version and 64 bit Windows are also available on
`appveyor <https://ci.appveyor.com/project/HazenBabcock/storm-analysis>`_. Under "Job Name",
click on the job with the appropriate PYVERSION, then click on "ARTIFACTS".

`Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_ is one source for Windows
wheels for this project's Python dependencies.

.. highlight:: none
	       
An example 64 bit Windows / Python3 installation ::

  # Create a virtual environment:

  $ \path\to\Python3.X\python.exe -m venv pyenv\sa_python3
  $ pyenv\sa_python3\Scripts\activate

  
  # Install dependencies, wheel files downloaded from Gohlke.

  $ pip install \path\to\wheels\xyz_amd64.whl
  $ ...
  $ pip install randomcolor
  $ pip install pyqt5

  
  # Install the storm-analysis wheel.

  $ pip install \path\to\wheels\storm_analysis-X.X-py3-none-any.whl
  

  # (Optional) test the install.

  $ pip install pytest
  $ cd pyenv\sa_python3\Lib\site-packages\storm_analysis\test
  $ pytest

.. note:: In this install the code in the storm_analysis project can be found in ``pyenv\sa_python3\Lib\site-packages\storm_analysis`` so if you would like to run something from the command line this is the place to look.
  
Installing from source
----------------------

Virtual environments
~~~~~~~~~~~~~~~~~~~~

In order to isolate the storm-analysis from other Python projects, or if you are attempting
to install this package on a computer where you do not have root access the use of Python
virtual environments is highly recommended. Two good resources are:

1. `Guide to Python / Virtual Environments <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_
2. `PyQt5 / mailing list <https://www.riverbankcomputing.com/pipermail/pyqt/2017-March/039032.html>`_

.. note:: Link (2) above describes the best way to create the virtual environment on Windows in a way that they will work with PyQt5.

Using setup.py
~~~~~~~~~~~~~~

The C libraries are built using `SCons <http://scons.org/>`_.

Basic installation ::
  
   $ git clone https://github.com/ZhuangLab/storm-analysis.git
   $ cd storm-analysis
   $ python setup.py build_c
   $ python setup.py install

You may find that this does not work because ``build_c`` fails. This step is just a
wrapper for SCons, so you may have better luck running the SCons by itself, then using
``python setup.py install`` to install the project.

Linux / OS-X example ::
  
  $ cd storm-analysis
  $ scons
  $ python setup.py install
  
Windows (mingw64) example ::

  $ cd storm-analysis
  $ C:\path\to\scons.bat -Q compiler=mingw
  $ python setup.py install

`nuwen <https://nuwen.net/mingw.html>`_ is one source for mingw64.

.. note:: The OS-X build assumes that the lapack and fftw libraries are installed in the standard homebrew location, /usr/local/. If this is not the case you may need to edit storm-analysis/SConstruct.

.. note:: The OS-X build requires a fairly recent version of XCode, v8.1+? v8.3.3 is known to work.
   
Using `Anaconda <https://www.anaconda.com/downloads>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(Optional) create an environment to keep your main Python installation clean ::

  $ conda create -n my_env python=X.Y
  $ conda activate my_env  # or activate my_env under Windows

Install dependencies (Linux / OS-X) ::

  $ conda config --add channels conda-forge 
  $ conda install gcc numpy pytest pytest-runner tifffile scipy matplotlib h5py astropy pillow shapely randomcolor pywavelets scons

.. note:: Anaconda gcc and XCode gcc may clash. If you have XCode installed you may have better luck not using gcc from Anaconda.

Install dependencies (Windows) ::

  $ conda config --add channels conda-forge 
  $ conda install numpy pytest pytest-runner m2w64-toolchain tifffile scipy h5py astropy matplotlib pillow shapely randomcolor pywavelets scons
  
Get the ``storm-analysis`` source code using git ::

  $ git clone https://github.com/ZhuangLab/storm-analysis.git
  $ cd storm-analysis

Install storm-analysis ::

  # Windows / mingw
  $ scons -Q compiler=mingw
  $ python setup.py install

  # Linux / OS-X
  $ scons
  $ python setup.py install
 
Testing
~~~~~~~

Test the (source) installation (this will take a few minutes to run).

Option 1 ::
    
  $ cd storm-analysis
  $ python setup.py test

Option 2 ::
  
  $ cd storm-analysis
  $ pytest

.. note:: Due to issues with creating pickle files that are compatible between Python2
	  and Python3 all of the tests that involve pickles (Spliner mostly) are skipped
	  on Python2.

Also
----

If you are modifying the code in the storm-analysis project you may find it more convenient
to add a .pth file to your pythonX.Y/site-packages directory. Then you won't have to
run ``python setup.py install`` after every change.
