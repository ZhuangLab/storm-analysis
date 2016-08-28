
This directory includes the 3D-DAOSTORM peak finder modified for analysis 
of data from a sCMOS camera, along with several utility programs. The
algorithms are based on those described in this reference:

"Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms"
F. Huang et al. Nature Methods, 10, p653-658.

The sCMOS code has the same dependencies as the 3D-DAOSTORM code. You
will need to follow the instructions in the 3D-DAOSTORM README.txt file
in order to get everything set up.


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
(execute this command in the sCMOS/sample_data directory)
python ../scmos_analysis.py sample.dax sample_mlist.bin 2d_fit.xml

If this works correctly you will see the following output:

Peak finding
 Removing negative values in frame 0
Frame: 0 31 31
 Removing negative values in frame 1
Frame: 1 35 66
 Removing negative values in frame 2
Frame: 2 31 97
 Removing negative values in frame 3
Frame: 3 44 141
 Removing negative values in frame 4
Frame: 4 34 175
 Removing negative values in frame 5
Frame: 5 36 211
 Removing negative values in frame 6
Frame: 6 41 252
 Removing negative values in frame 7
Frame: 7 39 291
 Removing negative values in frame 8
Frame: 8 37 328
 Removing negative values in frame 9
Frame: 9 48 376

Added 376

Tracking
Molecules: 376 (sample_mlist.bin)
Descriptor: 1
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 376 tracks
Analysis complete
