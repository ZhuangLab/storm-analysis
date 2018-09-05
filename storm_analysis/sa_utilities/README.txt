
Python Programs:

align_and_merge.py - Combine to localization files into a single localization file.

batch_analysis.py - A utility script for running multiple instances of 3d_daostorm
   (or sCMOS) at once.

bin_to_hdf5.py - Convert a .bin file to the storm-analysis HDF5 format.

bin_to_image.py - Create 2D and 3D images from a .bin file.

bin_to_lmchallenge_format.py - Convert a .bin file to the Single-Molecule Localization
   Microscopy challenge format as described here:
   http://bigwww.epfl.ch/smlm/methods/index.html?p=format

bin_to_PYME_h5r_format.py - Convert a .bin file to the PYME h5r format. These files
   can be rendered with the VisGUI program in PYME.

bin_to_lmchallenge_format.py - Convert a .bin file to the format used in the 2013
   localization microscopy challenge.
   
bin_to_tagged_spot_file.py - Convert a .bin file to the Micro-Manager tagged spot file
   (tsf) format. Tagged spot file format files can then be rendered with programs
   such as Micro-Manager and ViSP (1).

fiducials.py - Functions for fiducial tracking.

finding_fitting_error.py - Calculate the localization error for simulations where the
   the true locations are known.

fitz_c.py - This program is used to determine the z position from the localization x and y
   widths and a previously determined calibration curve.

hdf5_to_bin.py - Convert a storm-analysis HDF5 file to the Insight3 binary format.

hdf5_to_image.py - Create an image from a HDF5 file.

hdf5_to_txt.py - Convert a storm-analysis HDF5 file to a comma separated text file.

merge_bin.py - Merge two or more bin files into a single file.

merge_hdf5.py - Merge two or more (tracked) HDF5 files into a single file.

mortensen.py - Calculate X/Y localization accuracy Cramer-Rao bound as in Mortensen,
   Nature Methods, 2010.
   
read_tagged_spot_file.py - Read .tsf format file. This is useful mostly as a debugging
   aid to make sure that the .tsf file gotten written properly (1).

recall_fraction.py - Calculate the recall fraction (or noise fraction) of measured
   localizations compared to ground truth localizations.

reduce_mlist.py - Remove localizations from a .bin file that are outside of an AOI
   and/or minimum and maximum frame number.

std_analysis.py - This attempts to encapsulate the basic analysis functions such
   as averaging, driftCorrection, standardAnalysis, tracking and zFitting. All of
   these functions are used by both 3D-DAOSTORM and the sCMOS analysis.

tiffs_to_dax.py - This will create a .dax format file from a directory containing
   a .tif file for each frame in a movie. The .tif files are sorted by name and
   then added to the .dax file.

track_drift_correct.py - This does tracking, averaging and drift correction on a 
   localizations file. It can be useful if you want to redo these parts without
   redoing peak finding and fitting.

tracker.py - This program tracks objects across multiple frames & assigns the
   appropriate category to each object (i.e. specific or non-specific activation, etc.)
   
xyz_drift_correction.py - This will determine (but not apply) the drift correction
   for a STORM movie.


C Programs:

apply_drift_correction.c - C library for apply_drift_correction_c.py

avemlist.c - C library for avemlist_c.py

fitz.c - C library for fitz_c.py


(1) These programs require Google protocol buffers. I have also included the 
   (generated) TSFProto_pb2.py file that describes the tsf format. However, if the
   tsf format gets changed in Micro-Manager this may fall out of sync.

