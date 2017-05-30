
This can used to fit biplane or multi-plane SMLM data. It includes
the ability to analyze data from a sCMOS camera. It uses the fitting
approach of spliner (i.e. c-splines) adapted for multi-plane imaging.

If you want to use it to analyze data that is not from an sCMOS camera
then the likely easiest thing to do is create a dummy calibration
file. For better performance the gain term should be in units of
ADU / photo-electron as this analysis is designed to work in units of
photo-electrons.

Python Programs:

mapper - A PyQt5 GUI for identifying the mapping between different image
   planes.

measure_psf - Used to measure the PSF given the average z_stack and
   a text file with the z-offsets of each frame.

psf_localizations - Used to select localizations for PSF determination.
   These are localizations whose positions are inside all of the data
   planes, and that are not too close to each other.

psf_zstack - Used to create the average z_stack of the PSF localizations
   that is used for PSF measurement.


Analysis steps:

1. Acquire movie(s) with reasonably bright, small and well separated beads.

2. Analyze bead movies with sCMOS / 3D-DAOSTORM.

3. Identify (image) plane to (image) plane mapping using
   multi_plane/mapper.py.

4. Select good localizations (beads) to use for PSF measurement using
   multi_plane/psf_localizations.py.

5. Create 2x up-sampled average z stacks using for each plane using
   multi_plane/psf_zstack.py. This will also correct for pixel-wise
   differences in gain and offset.

6. Create z offsets file text file, possibly using spliner/offset_to_z.py.

7. Measure PSF for each plane using multi_plane/measure_psf.py.

8. Normalize PSFs relative to each other using multi_plane/normalize_psfs.py.

9. Create splines for each plane using spliner/psf_to_spline.py.

10. Acquire movie(s) of the sample of interest.

11. Create XML with multi_plane analysis parameters.

12. Perform multi-plane analysis with multi_plane/multi_plane.py.
