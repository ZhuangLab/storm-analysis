
The algorithms are based on those described in this reference:

"Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms"
F. Huang et al. Nature Methods, 10, p653-658.


Python Programs:

batch_analysis - Used to analyze an entire directory of movies with the
   same parameters.
   
camera_calibration - Used to generate a camera calibration file given the
   appropriate input data files.

check_dax_variance - Used to check if the mean and variance of the current
   movie match what is expected from the camera calibration.

dax_to_calib_format - Convert a dax file into a format that can be input
   into the camera_calibration program. Note that this is not the most
   efficient way to calibrate your camera.

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

$ python ../scmos_analysis.py --movie sample.dax --bin sample_mlist.bin --xml 2d_fit.xml 
Peak finding
 Removing values < 1.0 in frame 0
Frame: 1 70 70
 Removing values < 1.0 in frame 1
Frame: 2 81 151
 Removing values < 1.0 in frame 2
Frame: 3 84 235
 Removing values < 1.0 in frame 3
Frame: 4 92 327
 Removing values < 1.0 in frame 4
Frame: 5 88 415
 Removing values < 1.0 in frame 5
Frame: 6 81 496
 Removing values < 1.0 in frame 6
Frame: 7 86 582
 Removing values < 1.0 in frame 7
Frame: 8 87 669
 Removing values < 1.0 in frame 8
Frame: 9 87 756
 Removing values < 1.0 in frame 9
Frame: 10 96 852

Added 852
   0 fits lost due to Cholesky failure
   43 fits lost to image margin
   0 fits lost to negative value fit function
   0 fits lost to negative height
   0 fits lost to negative width

Tracking
Molecules: 852 (sample_mlist.bin)
Descriptor: 1
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 852 tracks
Analysis complete
