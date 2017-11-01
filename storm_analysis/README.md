
## Directory Layout ##

L1H - The code to perform l1H analysis as described in this publication [L1H](http://dx.doi.org/10.1364/OE.21.028583).

c_libraries - The location of the C libraries once they have been built. This also contains DLLs that you may need to get this package to work on windows (LAPACK, FFTW).

daostorm_3D - The code to perform 3D-DAOSTORM analysis as described in this publication [3D-DAOSTORM](http://dx.doi.org/10.1186/2192-2853-1-6).

dbscan - Density-based spatial clustering of applications with noise (DBSCAN) following [Ester et al](http://www.aaai.org/Papers/KDD/1996/KDD96-037).

decon-storm - Matlab, C and Python code for DeconSTORM analysis as described in this publication [DeconSTORM](http://dx.doi.org/10.1016/j.bpj.2012.03.070).

diagnostics - Code for evaluating the performance of different aspects of the different fitters, drift correction, etc..

fista - The code to perform FISTA image deconvolution following [Beck and Teboulle](http://dx.doi.org/10.1137/080716542).

frc - The code to perform FRC analysis following [Nieuwenhuizen et al](http://dx.doi.org/10.1038/nmeth.2448).

micrometry - Automatically find the affine transform between two localization files using geometric hashing.

multi_plane - Analysis of one or more planes of data from sCMOS camera(s), using one of the following PSF models (1) measured PSF, (2) pupil function or (3) 3D cubic spline.

psf_fft - The core code to analyze SMLM movies by fitting the measured PSF (using a FFT based approach).

pupilfn - The core code to analyze SMLM movies by fitting pupil functions.

rcc - The core code to perform RCC drift correction following [Wang et al](http://dx.doi.org/10.1364/OE.22.015982).

rolling_ball_bgr - The core code to perform rolling ball based background estimation.

sCMOS - The core code to perform analysis of data from a sCMOS camera following [Huang et al](http://dx.doi.org/10.1038/nmeth.2488).

sa_library - A collection of Python (and C) libraries that provide functionality whose utility is not limited to just a support role for 3D-DAOSTORM (and other analysis approaches).

sa_utilities - A collection of Python and C programs that perform functions such as tracking single molecules across multiple frames or determining and applying drift correction. These are used by 3D-DAOSTORM, but as with the functions in the library directory they are considered to be of more general utility.

simulator - A simple simulator for generating test data.

spliner - The core code to perform C-Spline analysis as described in [Babcock and Zhuang](http://dx.doi.org/10.1038/s41598-017-00622-w).

tests - Simple tests of 3D-DAOSTORM, sCMOS and other programs to verify that they work. This might also be a good place to look to get an idea of how different programs are run.

visualizer - A PyQt5 based application that draws the found localizations on frame in which they were found.

voronoi - Voronoi diagram based clustering following [Levet et al](http://dx.doi.org/10.1038/nmeth.3579).

wavelet_bgr - The core code to perform wavelet based background estimation following [Galloway et al](http://www.opticsinfobase.org/as/abstract.cfm?URI=as-63-12-1370).


## General Notes ##

Questions should be addressed to Hazen Babcock (hbabcock _at_ fas.harvard.edu)
