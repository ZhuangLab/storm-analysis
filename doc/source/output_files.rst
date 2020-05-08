.. highlight:: none

Output files
============

These are the output file formats.

XX.hdf5
-------

This is a HDF5 format file containing the localizations. The internal
structure is one group per frame analyzed. In the group each
localization property is saved as a separate dataset.

Localization Fields: ::

  background - Estimated background under the localization (e-).
  error - Fitting error chi-square value.
  height - Localization height (e-).
  iterations - Number of iterations to fit the localization.
  significance - Probability the localization is not noise in
                 units of sigma.
  sum - A sum over the counts in the localization (e-).
  track_id - The id of the track this localization was assigned
             to.
  x - X location in the image (pixels).
  xsigma - X sigma (pixels).
  y - Y location in the image (pixels).
  ysigma - Y sigma (pixels).
  z - Z location (microns).

  Some of the fields may not be present depending on the analysis
  algorithm. For example, only 3D-DAOSTORM and sCMOS will estimate
  xsigma and ysigma, and track_id won't be present if tracking was
  not performed.

  e- is photo-electrons. This assumes that the gain, etc.. for the
  camera have been entered correctly in the analysis XML.

Localizations that have been tracked and averaged together are
stored as tracks in groups with a maximum dataset size of 
'track_block_size'.

Track Fields: ::

  Tracks have the same fields as localizations and also the
  following.

  category - A number such as '0', '1' that is usually only relevant
             for multi-activator STORM. A value of '9' indicates
	     that the z value was out of the expected range.
  frame_number - The frame that the track started in.
  track_id - A unique ID for the track. This can be used to find all
             the localizations that were in the track.
  track_length - The length of the track in frames.

  In the track, most of the fields are the sum of the fields in the
  localizations that were assigned to the track. For example the
  track height will be the total of the heights of all the
  localizations in the track.

The metadata is stored in the 'metadata.xml' root dataset as a 
variable length unicode string. This is just a copy of the parameters
XML file that was used for the analysis.

Accessing the raw localizations. ::

  import storm_analysis.sa_library.sa_h5py as saH5Py
  
  with saH5Py.SAH5Reader(h5_name) as h5:
     for fnum, locs in h5.localizationsIterator():
        ...

Accessing the tracked localizations. ::

  with saH5Py.SAH5Reader(h5_name) as h5:
     for t in h5.tracksIterator():
        ...

Both iterators return the localizations as Python dictionaries
containing a numpy array for each property. If you only want
some of the properties you can specify this with the ``fields``
keyword argument. ::
  
  with saH5Py.SAH5Reader(h5_name) as h5:
     for fnum, locs in h5.localizationsIterator(fields = ['x', 'y']):
        ...

By default the ``localizationsIterator`` will skip empty frames
but you can disable this by setting the ``skip_empty`` keyword to
``False``. ::

  with saH5Py.SAH5Reader(h5_name) as h5:
     for fnum, locs in h5.localizationsIterator(skip_empty = False):
        ...

You can request the non-drift corrected positions by setting the
``drift_corrected`` keyword to ``False``. ::
  
  with saH5Py.SAH5Reader(h5_name) as h5:
     for fnum, locs in h5.localizationsIterator(drift_corrected = False):
        ...
	
You can get the movie size with ``getMovieInformation``. ::
  
  with saH5Py.SAH5Reader(h5_name) as h5:
     [movie_x, movie_y, movie_l, hash_value] = h5.getMovieInformation()
     ...
	
You can get the movie length with ``getMovieLength``. ::
  
  with saH5Py.SAH5Reader(h5_name) as h5:
     movie_l = h5.getMovieLength()
     ...

And you can get the pixel size (in nanometers) with ``getPixelSize``. ::

  with saH5Py.SAH5Reader(h5_name) as h5:
     pixel_size = h5.getPixelSize()
     ...
	
Input / Output
~~~~~~~~~~~~~~

Reading and writing of these files is handled by:

``storm_analysis/sa_library/sa_h5py.py``

This file can be converted to other formats using:

``storm_analysis/sa_utilities/hdf5_to_bin.py`` - to Insight3 format binary file.

``storm_analysis/sa_utilities/hdf5_to_txt.py`` - to CSV text file.


XX_drift.txt
------------

This is a text file containing the estimated x, y and z drift correction
values for each frame. It is not actually used as this information is also
saved in the localization XX.hdf5 file.

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

XX_list.bin (deprecated)
------------------------

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

This file can be converted to more standard formats using:

``storm_analysis/sa_utilities/bin_to_lmchallenge_format.py`` - to CSV text.

``storm_analysis/sa_utilities/bin_to_PYME_h5r_format.py`` - to `Python Microscopy Environment <http://www.python-microscopy.org/>`_.

``storm_analysis/sa_utilities/bin_to_tagged_spot_file.py`` - to `Tagged Spot File (tsf) <https://micro-manager.org/wiki/Tagged_Spot_File_(tsf)_format>`_.


