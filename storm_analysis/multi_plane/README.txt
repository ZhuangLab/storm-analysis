
This can used to fit biplane or multi-plane SMLM data, as well as
single plane data. It includes the ability to analyze data from a
sCMOS camera. It uses the fitting approach of spliner (i.e. c-splines)
adapted for multi-plane imaging.

If you want to use it to analyze data that is not from an sCMOS camera
then the likely easiest thing to do is create a dummy calibration
file. For better performance the gain term should be in units of
ADU / photo-electron as this analysis is designed to work in units of
photo-electrons.


Python Programs:

batch_heights - Performs all the steps necessary to create the
   heights.npy file for a single data set.

ch_mean_as_z - Calculate the first moment of the signal (height) as
   a function of color channel and store in the z field of the bin
   file. Note that the per channel data should be in the correct order,
   i.e. lowest to highest wavelength or vice-versa. Also, the data
   from the first channel will be used in the output file for all
   the other fields, so if necessary this file should have been
   drift corrected.

check_plane_offsets - Plots the PSF maximums as a function of z.

copy_tracking - Merge tracking data from channel 0 and localizations
   from channel N into a new file.

find_offset - Estimate the frame offset between movies from
   different cameras.

find_peaks_std - This does the peak finding and fitting.

kmeans_classifer - Given a codebook, classify the localizations into
   categories based on their relative signals (heights) in different
   color channels.

kmean_measure_codebook - Create a codebook for k-means classification
   from the localization heights from one or more data-sets and the
   known number of different dyes.

map_binfile - Splits the localizations in an i3 file into localizations
   for each plane based on a mapping. This is primarily a debugging tool,
   but it is also useful for mapping locations in other channels to
   channel 0.

mapper - A PyQt5 GUI for identifying the mapping between different image
   planes.

mapperView - A PyQt5 QGraphicsView specialized for mapper.

measure_psf - Used to measure the PSF given the average z_stack and
   a text file with the z-offsets of each frame.

merge_heights - Merge height information from localization files for
   different channels into a single (numpy) file.

mp_fit_c - The interface between Python and the C fitting library.

mp_utilities_c - A collection of utility functions used for multiplane
   analysis.

multi_plane - Run multi-plane / multi-color analysis.

normalize_psfs - Normalizes the PSFs for each plane so that they have
   the correct (relative) height.

plane_weighting - Calculates how to weight the updates from each plane
   when deciding how to update the localization position during fitting.

plot_heights - Plot a (sampling) of the heights from a single data
   set.

print_mapping - Prints out a mapping file.

psf_localizations - Used to select localizations for PSF determination.
   These are localizations whose positions are inside all of the data
   planes, and that are not too close to each other.

psf_zstack - Used to create the average z_stack of the PSF localizations
   that is used for PSF measurement.

split_peaks - Used mostly for debugging.

zstack_xydrift - Estimate the XY drift that occurred during the
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
