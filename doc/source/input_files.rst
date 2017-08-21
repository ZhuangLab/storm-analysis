.. highlight:: none
	       
Input (STORM movie) file formats
================================

The following are the image/movie formats that are supported. Loading
these files is handled by ``storm_analysis.sa_library.datareader.py``.

Dax
---

This is a Zhuang lab custom format. It is basically a raw format. Each
movie consists of two files, a binary **.dax** file that contains the data
as 16 bit integers and a text **.inf** file that includes the information
necessary to properly load the file.

A sample .inf file ::

  binning = 1 x 1
  data type = 16 bit integers (binary, big endian)
  frame dimensions = 256 x 256
  number of frames = 2
  x_start = 129
  x_end = 384
  y_start = 129
  y_end = 384  

Tiff
----

The standard tiff image format. Images should be gray-scale. `BigTIFF <http://bigtiff.org/>`_
images are also supported. IO is handled by the `tifffile <https://pypi.python.org/pypi/tifffile>`_
package.

Spe
---

This is a Roper Scientific binary format. We believe that our reader works
and that it supports the 16 and 32 bit integer formats and the 32 bit float
format, but we rarely test it.
