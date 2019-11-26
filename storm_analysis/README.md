
## Directory Layout ##

L1H - The code to perform l1H analysis as described in this publication [Babcock et al](http://dx.doi.org/10.1364/OE.21.028583).

admm - The code to perform ADMM image deconvolution following [Boyd et al](http://dx.doi.org/10.1561/2200000016), [ADMM](http://stanford.edu/~boyd/admm.html). The math for 3D ADMM follows [Ovesny et al](https://doi.org/10.1364/OE.22.031263).

c_libraries - The location of the C libraries once they have been built. This also contains DLLs that you may need to get this package to work on windows (LAPACK, FFTW).

daostorm_3D - The code to perform 3D-DAOSTORM analysis as described in this publication [Babcock et al](http://dx.doi.org/10.1186/2192-2853-1-6).

dbscan - Density-based spatial clustering of applications with noise (DBSCAN) following [Ester et al](http://www.aaai.org/Papers/KDD/1996/KDD96-037).

decon-storm - Matlab, C and Python code for DeconSTORM analysis as described in this publication [Mukamel et al](http://dx.doi.org/10.1016/j.bpj.2012.03.070).

densestorm - The code to perform 3denseSTORM image deconvolution following [Ovesny et al](https://doi.org/10.1364/OE.22.031263).

diagnostics - Code for evaluating the performance of different aspects of the different fitters, drift correction, etc..

fista - The code to perform FISTA image deconvolution following [Beck and Teboulle](http://dx.doi.org/10.1137/080716542).

frc - The code to perform FRC analysis following [Nieuwenhuizen et al](http://dx.doi.org/10.1038/nmeth.2448).

jupyter_examples - The code for creating sample data for the [Jupyter](http://jupyter.org/) notebook examples.

micrometry - Automatically find the affine transform between two localization files using geometric hashing following [Lang et al](http://dx.doi.org/10.1088/0004-6256/139/5/1782).

multi_plane - The code to analyze multiple image plane data as described in this publication [Babcock](http://dx.doi.org/doi:10.1038/s41598-018-19981-z).

psf_fft - The code to analyze SMLM movies by fitting the measured PSF (using a FFT based approach).

pupilfn - The code to analyze SMLM movies by fitting pupil functions.

rcc - The code to perform RCC drift correction following [Wang et al](http://dx.doi.org/10.1364/OE.22.015982).

rolling_ball_bgr - The code to perform rolling ball based background estimation.

sCMOS - The code to perform analysis of data from a sCMOS camera following [Huang et al](http://dx.doi.org/10.1038/nmeth.2488).

sa_library - A collection of Python (and C) libraries that provide functionality whose utility is not limited to just a support role for 3D-DAOSTORM (and other analysis approaches).

sa_utilities - A collection of Python and C programs that perform functions such as tracking single molecules across multiple frames or determining and applying drift correction. These are used by 3D-DAOSTORM, but as with the functions in the library directory they are considered to be of more general utility.

simulator - A simulator for generating test data.

slurm - Code for running storm-analysis in a distributed computing environment that uses [SLURM](https://slurm.schedmd.com/) for job management.

spliner - The code to perform C-Spline analysis as described in [Babcock and Zhuang](http://dx.doi.org/10.1038/s41598-017-00622-w).

test - Unit tests for the project.

visualizer - A PyQt5 based application that draws the found localizations on frame in which they were found.

voronoi - Voronoi diagram based clustering following [Levet et al](http://dx.doi.org/10.1038/nmeth.3579).

wavelet_bgr - The code to perform wavelet based background estimation following [Galloway et al](http://www.opticsinfobase.org/as/abstract.cfm?URI=as-63-12-1370).


## General Notes ##

Questions should be addressed to Hazen Babcock (hbabcock _at_ fas.harvard.edu)
