Modifying the Analysis
======================

In the process of developing this project we have attempted
to make it modular, hopefully making it easier to experiment
with different analysis approaches. The analysis is done in
two steps, peak finding and peak fitting.

Key Modules
-----------

1. ``storm_analysis/sa_library/dao_fit_c.py`` - Base class
   for interfacing with the C fitting program.

2. ``storm_analysis/sa_library/fitting.py`` - Base classes
   for both peak finding and fitting.

3. ``storm_analysis/sa_library/ia_utilities_c.py`` - A collection
   image analysis utilities such as identifying maxima in an image.

4. ``storm_analysis/sa_utilities/std_analysis.py`` - Performs all
   the steps in the analysis pipeline.
      
Other Modules
-------------

1. ``storm_analysis/sa_library/batch_run.py`` - A framework for
   running multiple python scripts in parallel.

2. ``storm_analysis/sa_library/imagecorrelation.py`` - Image cross
   correlation, primarily used for the purpose of drift correction.

3. ``storm_analysis/sa_library/matched_filter_c.py`` - Image convolution
   using an FFT based approach.

4. ``storm_analysis/sa_utilities/xyz_drift_correction.py`` - Drift
   correction using image cross-correlation.
