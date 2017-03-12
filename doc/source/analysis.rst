Analysis Programs
=================

These are the different localization finding and fitting approaches
in this project.

3D-DAOSTORM
-----------

This approach performs maximum likelihood estimation (MLE) guassian fitting.
It can be used to analyze 2D and 3D astigmatism STORM movies.

`Babcock et al <http://dx.doi.org/10.1186/2192-2853-1-6>`_

sCMOS
-----

This is many ways very similar to 3D-DAOSTORM, but it is designed to handle
the analysis of data from sCMOS cameras. In order for this to work well
you will need to have a calibration file containing the offset, gain
and variance for each camera pixel.

`Huang et al <http://dx.doi.org/10.1038/nmeth.2488>`_

Spliner
-------

This approach performs MLE fitting using a cubic spline approximation of
the microscope PSF. It can be used to analyze both 2D and 3D STORM movies
with arbitrary PSF shapes. In order to use it you will need to have
a fairly accurate measurement of your microscope PSF.

`Babcock and Zhuang <http://dx.doi.org/10.1101/083402>`_

L1H
---

This is a compressed sensing approach. It is substantially slower than
all of the above approaches and only works with 2D STORM movies. If your
localization data is very high it may be a better choice.

`Babcock et al <http://dx.doi.org/10.1364/OE.21.028583>`_
