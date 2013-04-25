
Description:

STORM image analysis using Maximum Likelihood Deconvolution
(deconSTORM).

At the moment, this does not include the official version
that was used for the paper (and written in matlab). This is a
C & Python implementation of the image deconvolution procedure
following Lingenfelter 2009. It enables you to do deconvolution
with a fixed compression parameter (i.e. not spatially varying
across the image). The machinery is in place for spatially
(and temporarily) varying compression parameters but that
remains un-implemented.


Files:

fixed_compression.py - A python script for performing deconvolution
   with fixed compression on spe format STORM movie.

mlem_c.py - A python front-end for the mlem_sparse C library.

mlem_sparse.c - C library that does the actual deconvolution.


Notes:

You need to have the library directory in your python-path for
this to work.
