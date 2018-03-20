
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

movie_to_calib_format - Convert a movie file into a format that can be input
   into the camera_calibration program. Note that this is not the most
   efficient way to calibrate your camera.

reslice_calibration - If your camera calibration data (as output by the
   the camera_calibration program) covers a larger area of the chip than
   your movie you can use this slice it down.

scmos_analysis - Use this to perform peak finding on your sCMOS movie.


