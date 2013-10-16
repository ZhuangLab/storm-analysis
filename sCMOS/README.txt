
This directory includes a sparse peak finder designed for analysis of
data from a sCMOS camera, along with several utility programs. The
algorithms are based on those described in this reference:

"Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms"
F. Huang et al. Nature Methods, 10, p653-658.

The sCMOS code has the same dependencies as the 3D-DAOSTORM code. You
will need to follow the instructions in the 3D-DAOSTORM README.txt file
in order to get everything set up.

This is still a work in progress, and may not even be implemented 
correctly as I'm still puzzling over the exact meaning of equation 3.4
in the supplement.


Python Programs:

camera_calibration - Used to generate a camera calibration file given the
   appropriate input data files.

check_dax_variance - Used to check if the mean and variance of the current
   movie match what is expected from the camera calibration.

dax_to_calib_format - Convert a dax file into a format that can be input
   into the camera_calibration program. Note that this is not the most
   efficient way to calibrate your camera.

find_peaks - This does the actual peak finding. It is pretty similar to
   the file of the same name in the 3D-DAOSTORM directory.

reslice_calibration - If your camera calibration data (as output by the
   the camera_calibration program) covers a larger area of the chip than
   your movie you can use this slice it down.

scmos_analysis - Use this to perform peak finding on your sCMOS movie.


Sample Data:

These files are located in sCMOS/sample_data.

2d_fit.xml - The parameters to use for the analysis. This is pretty
   similar to 3D-DAOSTORM.

calib.npy - The camera calibration data for the area of the camera
   where the STORM data was taken.

sample.dax - A 10 frame movie of some STORM data from a Hamamatsu sCMOS
   camera.

sample_smooth.dax - The same data convolved as described in section 3.1
   (equation 3.1) of the supplement.


A sample run:
(execture this command in the sCMOS/sample_data directory)
python ../scmos_analysis.py sample.dax sample_mlist.bin 2d_fit.xml

If this works correctly you see the following output:

> Peak finding
>  Removing negative values in frame 0
> Frame: 0 34 34
>  Removing negative values in frame 1
> Frame: 1 38 72
>  Removing negative values in frame 2
> Frame: 2 33 105
>  Removing negative values in frame 3
> Frame: 3 46 151
>  Removing negative values in frame 4
> Frame: 4 37 188
>  Removing negative values in frame 5
> Frame: 5 37 225
>  Removing negative values in frame 6
> Frame: 6 44 269
>  Removing negative values in frame 7
> Frame: 7 40 309
>  Removing negative values in frame 8
> Frame: 8 38 347
>  Removing negative values in frame 9
> Frame: 9 48 395
>
> Added 395
>
> Tracking
> Molecules: 395 (sample_mlist.bin)
> Processing molecule 0 in frame 0 (tracker)
> Finished processing
> Found 395 tracks
> 
> Analysis complete
