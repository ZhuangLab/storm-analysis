
Python Programs:

 arraytoimage.py - For creating images from numpy arrays (using the PIL
   image library).

 datareader.py - For reading various kinds of STORM movie data. This can
   read the Zhuang lab .dax format as well as .tiff and .spe (Roper
   Scientific) files.

 daxwriter.py - For writing .dax files.

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

 i3dtype.py - Definitions of the Insight3 localization storage structure as
   well as some utility functions.

 i3togrid.py - For histogramming (or gridding) localizations that are stored
   in Insight3 format binary files.

 ia_utities_c.py - A python interface to the C ia_utilities library.

 imagecorrelation.py - For correlating images (or image stacks) to determine
   their offset in XY (and possibly Z).

 loadclib.py - Loads the correct C library depending on the OS, etc.

 matched_filter_c.py - A python interface to the C matched_filter library.
 
 multi_fit_c.py - A python interface to the C multi_fit library.

 parameters.py - For parsing simple xml files such as those that describe 
   how to perform the analysis.

 readinsight3.py - For reading Insight3 format binary files.

 readhres.py - For reading hres format binary files.

 rebin.py - FFT based image resizing.

 recenterPSF.py - Recenters a PSF for FFT image convolution.
 
 regfilereader.py - For reading transformation files generated with the
   ImageJ MultiStackReg plugin. These transformations can be applied to
   the localizations using I3GData class defined in the i3togrid.py file.

 writeinsight3.py - For writing Insight3 format binary files.


C Libraries:
 grid - C functions for gridding 2D and 3D data.

 ia_utilities - A collection of C image analysis utility functions used by 3D-DAOSTORM.

 matched_filter.c - C code for FFT based image convolution.
 
 multi_fit - The core C localization fitting routines.
