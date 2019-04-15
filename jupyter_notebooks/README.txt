
Notebook descriptions.

algorithm_resolution_limit - Analyze the resolution limit of an SLML analyis
                             algorithm.

dao3d_2d - Basic analysis with 3D-DAOSTORM. A good place to start as the other
           analysis programs work in a similar manner.

dao3d_zcal - How to measure z defocus curves for 3D astigmatism imaging
             using 3D-DAOSTORM.

dask_multiprocessing - Using Dask to parallelize 3D-DAOSTORM.

dbscan_clustering - Using DBSCAN to cluster localizations.

estimating_2D_fitting_precision - How to estimate the average 2D localization
                                  precision.

estimating_camera_parameters - One approach to estimating the gain and offset
                               given data from an uncalibrated camera.

fiducial_tracking - How to do fiducial tracking for drift correction.

gauss_fitting_cramer_rao - Calculating resolution for 2D Gaussian fitting
                           based on Cramer-Rao bounds.

micrometry_mapping - Using micrometry to find the best affine transform between
                     two images.

multiplane_mapping - How to measure a mapping file for multiplane analysis.

multiplane_measure_psf - How to measure PSFs for multiplane analysis.

multiplane_psfs_to_splines - Converting PSFs to splines for multiplane analysis.

scmos_cal - How to calibrate your sCMOS camera.

smlm_challenge_2013_dao3d_2dfixed - Analyze a 2013 SMLM challenge dataset with
                                    3D-DAOSTORM.

spliner_measure_psf - How to measure a PSF for Spliner.

voronoi_clustering - Using Voronoi diagrams to cluster localizations.


A Docker image with Jupyter, Python3 and storm-analysis is also available.

$ docker pull zhuanglab/jupyter-sa:latest
