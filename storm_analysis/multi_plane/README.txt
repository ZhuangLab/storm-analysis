
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


Analysis steps:

1. Acquire movie(s) with reasonably bright, small and well separated beads.

   a. Alignment movie - Hop back and forth between the different focal
      planes. This is for calculating the plane to plane mapping.

   b. PSF measurement movie - This should be a scan from the minimum to
      to the maximum z value that you want to allow in the analysis.

2. Analyze bead movies with sCMOS / 3D-DAOSTORM.

3. Identify (image) plane to (image) plane mapping using
   multi_plane/mapper.py.

3a. Alternatively use micrometry/micrometry.py and micrometry/merge_maps.py
   (if necessary).
   
4. Select good localizations (beads) to use for PSF measurement using
   multi_plane/psf_localizations.py.

5. Create 2x up-sampled average z stacks using for each plane using
   multi_plane/psf_zstack.py. This will also correct for pixel-wise
   differences in gain and offset. AOI size ~12.

6. Create z offsets file text file, possibly using spliner/offset_to_z.py.

7. Measure PSF for each plane using multi_plane/measure_psf.py. For multi-color
   analysis use --normalize True and skip the next step (step 8).

8. Normalize PSFs relative to each other using multi_plane/normalize_psfs.py.

9. (Optional) check plane z offsets using multi_plane/check_plane_offsets.py.
   If the offsets are not well centered it can be adjusted using the --deltaz
   argument to spliner/offset_to_z.py and restarting at step 6.

10. Create splines for each plane using spliner/psf_to_spline.py. Spline size ~20.

11. Create XML with multi_plane analysis parameters. Use a value of 1 for
    independent_heights for multi-color analysis.

12. Calculate weights for parameters as a function of plane and z using
    multi_plane/plane_weighting.py.

13. Acquire movie(s) of the sample of interest.

14. Determine frame offsets between movies from different cameras (if any)
    using multi_plane/find_offsets.py

15. Perform multi-plane analysis with multi_plane/multi_plane.py.


(Perform the following additional steps for multi-color analysis only).

16. Create files with height information from each channel using
    multi_plane/batch_heights.

17. Create a localization file with the z value replaced by the first
    moment of the signal (height) as a function of the color channel
    with multi_plane/ch_mean_as_z.py.

18. Create a codebook for k-means classification of localizations using
    multi_plane/kmeans_measure_codebook.py.

19. Categorize localizations with a codebook using
    multi_plane/kmeans_classifier.py


Also:

1. For multi-color measurements it may be helpful to have data from
relatively sparse single color antibodies. This is useful for creating
a codebook for k-means classification as well as developing some
sense of the amount of cross-talk between channels.

2. If you need to merge the results from multiple movies you can do this
with sa_utilities/align_and_merges.py (corrects for x,y,z offsets), or
sa_utilities/merge_bin.py (does not correct for x,y,z offset).
