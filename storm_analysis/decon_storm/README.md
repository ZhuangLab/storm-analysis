# DeconSTORM #
This software is presented as supplementary material for the publication, "Statistical Deconvolution for Super-Resolution Fluorescence Microscopy," by Eran A. Mukamel, Hazen Babcock and Xiaowei Zhuang, Biophysical Journal (2012). Please see the paper for a full explanation of the method.

Documentation in HTML format is in the folder "doc".

Copyright, 2012
Eran Mukamel, Hazen Babcock and Xiaowei Zhuang

# Description #
STORM image analysis using Maximum Likelihood Deconvolution
(deconSTORM).

# Files #
The DeconSTORM software as described in the publication is contained in the Matlab files.

fixed_compression.py - A python script for performing deconvolution with fixed compression on spe format STORM movie.

mlem_c.py - A python front-end for the mlem_sparse C library.

mlem_sparse.c - C library that does the actual deconvolution.

# Notes #
You need to have the library directory in your python-path for
this to work.
