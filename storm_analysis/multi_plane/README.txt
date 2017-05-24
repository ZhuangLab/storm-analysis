
This can used to fit biplane or multi-plane SMLM data. It includes
the ability to analyze data from a sCMOS camera.


Python Programs:

mapper - A PyQt5 GUI for identifying the mapping between different image
   planes.

psf_localizations - Used to select localizations for PSF determination.
   These are localizations whose positions are inside all of the data
   planes, and that are not too close to each other.

psf_zstack - Used to create the average z_stack of the PSF localizations
   that is used for PSF measurement.


Analysis steps:

1. Acquire movie(s) with reasonably bright & well separated beads.

2. Analyze bead movies with sCMOS / 3D-DAOSTORM.

3. Identify (image) plane to (image) plane mapping using
   multi_plane/mapper.py.

4. Select good localizations (beads) to use for PSF measurement using
   multi_plane/psf_localizations.py.

5. Create 2x up-sampled average z stacks using for each plane using
   multi_plane/psf_zstack.py.

6. Create z offsets file text file, possibly using spliner/offset_to_z.py.

7. Measure PSF for each plane using multi_plane/measure_psf.py.

8. Create splines for each plane using spliner/psf_to_spline.py.

9. Acquire movie(s) of the sample of interest.

10. Create XML with multi_plane analysis parameters.

11. Perform multi_plane analysis with multi_plane/multi_plane.py.

