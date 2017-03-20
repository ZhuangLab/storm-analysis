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

Installing using wheels
-----------------------

Wheels for 64 bit Windows are `here <http://zhuang.harvard.edu/storm_analysis/>`_.
Despite their names, these will not work with 32 bit Python as the C libraries are 64bit.

Wheels for the latest version and 64 bit Windows are also available on
`appveyor <https://ci.appveyor.com/project/HazenBabcock/storm-analysis>`_. Under "Job Name",
click on the job with the appropriate PYVERSION, then click on "ARTIFACTS".

`Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_ is one source for Windows
wheels for this project's Python dependencies.

Installing from source
----------------------

Using setup.py
~~~~~~~~~~~~~~

The C libraries are built using `SCons <http://scons.org/>`_. Note that SCons does not
work with Python3 so you'll need to Python2 to build the C libraries.

Basic installation ::
  
   git clone https://github.com/ZhuangLab/storm-analysis.git
   cd storm-analysis
   python setup.py build_c
   python setup.py install

You may find that this does not work because ``build_c`` fails. This step is just a
wrapper for SCons, so you may have better luck running the SCons by itself, then using
``python setup.py install`` to install the project.

Linux example ::
  
  cd storm-analysis
  scons
  python setup.py install
  
Windows (mingw64) example ::

  cd storm-analysis
  C:\path\to\scons.bat -Q compiler=mingw
  python setup.py install

`nuwen <https://nuwen.net/mingw.html>`_ is one source for mingw64.

Using `Anaconda <https://www.continuum.io/downloads>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(Optional) create an environment to keep your main Python installation clean ::

  conda create -n my_env python=X.Y
  source activate my_env  # or activate my_env under Windows

Install dependencies (Linus/ OS-X) ::

  conda config --add channels conda-forge 
  conda install numpy pytest pytest-runner gcc
  conda install tifffile scipy matplotlib
  conda install pillow shapely randomcolor pywavelets

Install dependencies (Windows) ::

  conda config --add channels conda-forge 
  conda install numpy pytest pytest-runner
  conda install m2w64-toolchain tifffile scipy
  conda install matplotlib pillow shapely randomcolor pywavelets

Get the ``storm-analysis`` source code using git ::

  git clone https://github.com/ZhuangLab/storm-analysis.git
  cd storm-analysis

Python2 ::

  conda install scons

  # Windows / mingw
  scons -Q compiler=mingw
  python setup.py install

  # Linux
  scons
  python setup.py install

Python3 (this requires that you also have Python2 installed for SCons) ::

  # Windows / mingw	
  C:\path\to\scons.bat -Q compiler=mingw
  python setup.py install

  # Linux
  scons                                   
  python setup.py install
 
Testing
-------

Test the installation (this will take a few minutes to run)

Option 1 ::
  
  cd storm-analysis/storm_analysis/test
  nose2

Option 2 ::
    
  cd storm-analysis
  python setup.py test

.. note:: Due to issues with creating pickle files that are compatible across multiple OSs and versions of Python some of the tests may fail on Windows. They are all expected to pass on Linux.
