
Python Programs:

apply_drift_correction_c.py - Applies a previously determined drift correction to each
   of the localizations.

avemlist_c.py - This program averages all the objects in a track into a single object.

batch_analysis.py - A utility script for running multiple instances of 3d_daostorm
   (or sCMOS) at once.

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

finding_fitting_error.py - Calculate the localization error for simulations where the
   the true locations are known.

fitz_c.py - This program is used to determine the z position from the localization x and y
   widths and a previously determined calibration curve.
   
read_tagged_spot_file.py - Read .tsf format file. This is useful mostly as a debugging
   aid to make sure that the .tsf file gotten written properly (1).

std_analysis.py - This attempts to encapsulate the basic analysis functions such
   as averaging, driftCorrection, standardAnalysis, tracking and zFitting. All of
   these functions are used by both 3D-DAOSTORM and the sCMOS analysis.

tiffs_to_dax.py - This will create a .dax format file from a directory containing
   a .tif file for each frame in a movie. The .tif files are sorted by name and
   then added to the .dax file.

track_average_correct.py - This does tracking, averaging and drift correction on a 
   molecule list file. It can be useful if you want to redo these parts without
   redoing the localization step, which takes most of the time.

tracker_c.py - This program tracks objects across multiple frames & assigns the
   appropriate category to each object (i.e. specific or non-specific activation, etc.)
   
xyz-drift-correction.py - This will determine (but not apply) the drift correction
   for a STORM movie.


C Programs:

apply_drift_correction.c - C library for apply_drift_correction_c.py

avemlist.c - C library for avemlist_c.py

fitz.c - C library for fitz_c.py

tracker.c - C library for tracker_c.py


(1) These programs require Google protocol buffers. I have also included the 
   (generated) TSFProto_pb2.py file that describes the tsf format. However, if the
   tsf format gets changed in Micro-Manager this may fall out of sync.

