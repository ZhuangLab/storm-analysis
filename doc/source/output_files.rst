.. highlight:: none

Output files
============

These are the output file formats.

XX_list.bin
-----------

This is a binary file format containing the localizations. Typically
there will be two of these, the XX_mlist.bin file contains all the
localizations in every frame. In the XX_alist.bin file however
localizations that were identified as spanning multiple frames are
averaged together into a single localization.

Format
~~~~~~

It has the following format: ::

  # header (16 bytes)
  version - 4 byte string, this should be "M425".
  frames - int32, this is not the number of frames, it should be 1.
  status - int32, this should be 6..
  molecules - int32, the number of localizations in the file.

  # data record (there are molecules number of these in the file).
  x - float32, x location
  y - float32, y location
  xc - float32, drift corrected x location
  yc - float32, drift corrected y location
  h - float32, fit height
  a - float32, fit area
  w - float32, fit width
  phi - float32, fit angle (for unconstrained elliptical gaussian)
  ax - float32, peak aspect ratio
  bg - float32, fit background
  i - float32, sum - baseline for pixels included in the peak
  c - int32, peak category ([0..9] for STORM images)
  fi - int32, fit iterations
  fr - int32, frame
  tl - int32, track length
  lk - int32, link (id of the next molecule in the trace)
  z - float32, original z coordinate
  zc - float32, drift corrected z coordinate

  # footer (4 bytes)
  na - int32, this is 0.

  # meta-data as an XML string (optional, not all files will have this).
  <xml>..</xml>

It is important to note that the analysis programs may not set all of these
use fields and may use some of them for different purposes. In particular,
3D-DAOSTORM, sCMOS and Spliner make the following changes. ::

  fi - The fit status.
  i - The fit error.
  
DBSCAN and Voronoi cluster identification make the following changes. ::

  a - The number of localizations in the cluster.
  lk - The cluster ID.
  fr - This is also the cluster ID.

Input / Output
~~~~~~~~~~~~~~

Reading and writing of these files is handled by:

``storm_analysis/sa_library/readinsight3.py``
``storm_analysis/sa_library/writeinsight3.py``

The numpy data type for these files is defined here:

``storm_analysis/sa_library/i3dtype.py``

This file can be converted to more standard format using:

``storm_analysis/sa_utilities/bin_to_lmchallenge_format.py`` - to CSV text.

``storm_analysis/sa_utilities/bin_to_PYME_h5r_format.py`` - to `Python Microscopy Environment <http://www.python-microscopy.org/>`_.

``storm_analysis/sa_utilities/bin_to_tagged_spot_file.py`` - to `Tagged Spot File (tsf) <https://micro-manager.org/wiki/Tagged_Spot_File_(tsf)_format>`_.


XX_drift.txt
------------

This is a text file containing the estimated x, y and z drift correction
values for each frame.

Format
~~~~~~

The file is tab delimited with the following columns: frame number (1 indexed),
x offset (pixels), y offset (pixels), z offset (nanometers).

An example: ::
  
  1	-0.047	-0.056	0.000
  2	-0.047	-0.056	0.000
  3	-0.047	-0.056	0.000
  4	-0.047	-0.055	0.000
  5	-0.047	-0.055	0.000
  6	-0.046	-0.055	0.000
  7	-0.046	-0.055	0.000
  8	-0.046	-0.055	0.000
  9	-0.046	-0.054	0.000
  10	-0.046	-0.054	0.000

Input / Output
~~~~~~~~~~~~~~

These files are created by:

``storm_analysis/sa_utilities/xyz_drift_correction.py``
``storm_analysis/rcc/rcc_drift_correction.py``

And used by:

``storm_analysis/sa_utilities/apply_drift_correction_c.py``


XX.hres
-------

This is a binary output file created by L1H. It is a compressed version of the
high resolution image that L1H creates. Only pixels with non-zero values are
recorded. ::

  # header (100 bytes)
  x size - int32, image x size.
  y size - int32, image y size.

  # data record (12 bytes, repeats to the end of the file).
  fr - int32, frame number.
  i - pixel offset in the frame (as if the frame was a 1D array).
  z - pixel intensity.

Input / Output
~~~~~~~~~~~~~~

Reading of these files is handled by:

``storm_analysis/sa_library/readhres.py``
