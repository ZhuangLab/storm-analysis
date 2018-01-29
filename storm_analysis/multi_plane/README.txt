
This can used to fit biplane or multi-plane SMLM data. It includes the
ability to analyze data from a sCMOS camera. It can use any of the
following PSF models : splines, pupil functions, measured PSF
(multi_plane.py) and Gaussians (multi_plane_dao.py).

If you want to use it to analyze data that is not from an sCMOS camera
then the likely easiest thing to do is create a dummy calibration
file. For better performance the gain term should be in units of
ADU / photo-electron as this analysis is designed to work in units of
photo-electrons.


Python Programs:

channel_color.py - Calculate the first moment of the signal (height) as
   a function of color channel.

check_plane_offsets.py - Plots the PSF maximums as a function of z.

find_offset.py - Estimate the frame offset between movies from
   different cameras.

find_peaks_std.py - Configuration functions for multiplane peak finding
   and fitting using the Spliner, Pupil Function or PSF FFT PSF models.

find_peaks_std_dao.py - Configuration functions for multiplane peak
   finding and fitting using Gaussians (2D, variable width).

fitting_mp.py - sa_library.fitting sub-classed for multiple planes,
   this does the peak finding and fitting.

kmeans_classifer.py - Given a codebook, classify the localizations into
   categories based on their relative signals (heights) in different
   color channels.

kmean_measure_codebook.py - Create a codebook for k-means classification
   from the localization heights from one or more data-sets and the
   known number of different dyes.

mapper.py - A PyQt5 GUI for identifying the mapping between different image
   planes. (Deprecated, use micrometry.micrometry).

mapperView.py - A PyQt5 QGraphicsView specialized for mapper. (Deprecated, use
   micrometry.micrometry).

measure_psf.py - Used to measure the PSF given the average z_stack and
   a text file with the z-offsets of each frame.

mp_fit_c.py - The base interface between Python and the C fitting library.

mp_fit_arb_c.py - The interface between Python and the C fitting library
   for "arbitrary" PSF fitting (Spliner, etc..).
   
mp_fit_dao_c.py - The interface between Python and the C fitting library
   for Gaussian fitting.

mp_utilities.py - A collection of utility functions used for multiplane
   analysis.

multi_plane.py - Run multi-plane "arbitrary" PSF analysis.

multi_plane_dao.py - Run multi-plane Gaussian PSF analysis.

normalize_psfs.py - Normalizes the PSFs for each plane so that they have
   the correct (relative) height.

plane_weighting.py - Calculates how to weight the updates from each plane
   when deciding how to update the localization position during fitting.

plot_heights.py - Plot a (sampling) of the heights from a single data
   set.

print_mapping.py - Prints out a mapping file.

psf_localizations.py - Used to select localizations for PSF determination.
   These are localizations whose positions are inside all of the data
   planes, and that are not too close to each other.

psf_zstack.py - Used to create the average z_stack of the PSF localizations
   that is used for PSF measurement.

separate_channels.py - Used mostly for debugging. Splits the localizations
   in a HDF5 file out by channel.

zstack_xydrift.py - Estimate the XY drift that occurred during the
   acquisition of a PSF bead z stack measurement.


