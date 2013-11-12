# storm-analysis #
This is a repository of code developed in the [Zhuang Lab](http://zhuang.harvard.edu/) for analysis of STORM movies.

# Installation #
You will need Python and you will need the storm-analysis directory in your Python path.

Pre-compiled libraries (using MinGW) are provided for 64-bit Windows. For other operating systems you will need to compile the libraries and executables yourself as described in the README files in the relevant folders. There are also compile_bat.bat files which have examples of how to compile everything.

To use the software on 64-bit Windows you will need the windows_dll directory to be in your PATH.

# Directory Layout #
3d_daostorm - The core code to perform 3D-DAOSTORM analysis as described in this publication [3D-DAOSTORM](http://dx.doi.org/10.1186/2192-2853-1-6)

L1H - The core code to perform l1H analysis as described in a forthcoming publication.

decon-storm - Matlab, C and Python code for DeconSTORM analysis as described in this publication [DeconSTORM](http://dx.doi.org/10.1016/j.bpj.2012.03.070)

sCMOS - The core code to perform analysis of data from a sCMOS camera as described in this publication [Huang et al](http://dx.doi.org/10.1038/nmeth.2488)

sa_library - A collection of Python (and C) libraries that provide functionality whose utility is not limited to just a support role for 3D-DAOSTORM (and other analysis approaches).

sa_utilities - A collection of Python and C programs that perform functions such as tracking single molecules across multiple frames or determining and applying drift correction. These are used by 3D-DAOSTORM, but as with the functions in the library directory they are considered to be of more general utility.

visualizer - A PyQt4 based application that draws the found localizations on frame in which they were found.

windows_dll - The other DLLs that you will need to get the 3D-DAOSTORM analysis to work on windows.
