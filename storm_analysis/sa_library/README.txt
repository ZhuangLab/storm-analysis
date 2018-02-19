
Python Programs:

 affine_transform_c.py - A python interface to the C affine_transform library.

 analysis_io.py - Classes used for analysis IO by 3D-DAOSTORM, sCMOS,
   Spliner and Multiplane.
 
 arraytoimage.py - (Deprecated) For creating images from numpy arrays (using
   the PIL image library).

 batch_run.py - For running multiple Python instances in parallel

 cs_decon.py - Base class for compressed sensing deconvolution.
 
 cs_decon_utilities_c.py - A Python interface to the C cs_decon_utilities
   library.
 
 dao_fit_c.py - Python interface to the C fitting library.

 datareader.py - For reading various kinds of STORM movie data. This can
   read the Zhuang lab .dax format as well as .tiff and .spe (Roper
   Scientific) files.

 datawriter.py - For writing .dax files and .tif file. This is mostly
   used by the simulator.

 driftutilities.py - Functions that are common among the drift correction
   approaches.

 fitting.py - Contains the basic functions for peak fitting as well as the
   base class for storing the parameters, arrays and peaks used during the
   fitting. This provides the core functionality of the localization based
   analysis methods such as 3D-DAOSTORM and sCMOS.

 gaussfit.py - For gaussian (and lorentzian) fitting. Since this is pure
   Python it is a bit slow and is not generally used to fit localizations
   directly.

 grid_c.py - A python interface to the C grid library.

 i3dtype.py - (Deprecated) Definitions of the Insight3 localization storage
   structure as well as some utility functions.

 i3togrid.py - (Deprecated) For histogramming (or gridding) localizations
   that are stored in Insight3 format binary files.

 ia_utities_c.py - A python interface to the C ia_utilities library.

 imagecorrelation.py - For correlating images (or image stacks) to determine
   their offset in XY (and possibly Z).

 loadclib.py - Loads the correct C library depending on the OS, etc.

 matched_filter_c.py - A python interface to the C matched_filter library.

 parameters.py - For parsing xml files that describe how to perform the
   analysis. This includes a list of all the valid parameters and
   documentation of the purpose of each parameter.

 readinsight3.py - (Deprecated) For reading Insight3 format binary files.

 readhres.py - For reading hres format binary files.

 rebin.py - FFT based image resizing.

 recenterPSF.py - Recenters a PSF for FFT image convolution.
 
 regfilereader.py - (Deprecated) For reading transformation files generated
   with the ImageJ MultiStackReg plugin. These transformations can be applied to
   the localizations using I3GData class defined in the i3togrid.py file.

 sa_h5py.py - For reading and writing storm-analysis format HDF5 files.
 
 writeinsight3.py - (Deprecated) For writing Insight3 format binary files.


C Libraries:

 affine_transform.c - C functions for affine image transformation.

 cs_decon_utilities.c - C functions for post-analysis of images that have
   been deconvolved using compressed sensing approaches (FISTA, ADMM).
 
 dao_fit.c - C Gaussian fitting routines.
 
 grid.c - C functions for gridding 2D and 3D data.

 ia_utilities.c - A collection of C image analysis utility functions used by 3D-DAOSTORM.

 matched_filter.c - C code for FFT based image convolution.
 
 multi_fit.c - The core C localization fitting routines.
