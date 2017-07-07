
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

check_plane_offsets - Plots the PSF maximums as a function of z.

find_peaks_std - This does the peak finding and fitting.

map_binfile - Splits the localizations in an i3 file into localizations
   for each plane based on a mapping. This is primarily a debugging tool.

mapper - A PyQt5 GUI for identifying the mapping between different image
   planes.

mapperView - A PyQt5 QGraphicsView specialized for mapper.

measure_psf - Used to measure the PSF given the average z_stack and
   a text file with the z-offsets of each frame.

mp_fit_c - The interface between Python and the C fitting library.

mp_utilities_c - A collection of utility functions used for multiplane
   analysis.

normalize_psfs - Normalizes the PSFs for each plane so that they have
   the correct (relative) height.

plane_weighting - Calculates how to weight the updates from each plane
   when deciding how to update the localization position during fitting.

print_mapping - Prints out a mapping file.

psf_localizations - Used to select localizations for PSF determination.
   These are localizations whose positions are inside all of the data
   planes, and that are not too close to each other.

psf_zstack - Used to create the average z_stack of the PSF localizations
   that is used for PSF measurement.


Analysis steps:

1. Acquire movie(s) with reasonably bright, small and well separated beads.

   a. Alignment movie - Hop back and forth between the different focal
      planes. This is for calculating the plane to plane mapping.

   b. PSF measurement movie - This should be a scan from the minimum to
      to the maximum z value that you want to allow in the analysis.

2. Analyze bead movies with sCMOS / 3D-DAOSTORM.

3. Identify (image) plane to (image) plane mapping using
   multi_plane/mapper.py.

4. Select good localizations (beads) to use for PSF measurement using
   multi_plane/psf_localizations.py.

5. Create 2x up-sampled average z stacks using for each plane using
   multi_plane/psf_zstack.py. This will also correct for pixel-wise
   differences in gain and offset. AOI size ~12.

6. Create z offsets file text file, possibly using spliner/offset_to_z.py.

7. Measure PSF for each plane using multi_plane/measure_psf.py.

8. Normalize PSFs relative to each other using multi_plane/normalize_psfs.py.

9. Create splines for each plane using spliner/psf_to_spline.py. Spline size ~20.

10. Create XML with multi_plane analysis parameters.

11. Calculate weights for parameters as a function of plane and z using
    multi_plane/plane_weighting.py.

12. Acquire movie(s) of the sample of interest.

13. Perform multi-plane analysis with multi_plane/multi_plane.py.
