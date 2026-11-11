
This program lets you visualize the results of localization fitting.
It overlays the localizations on the frame of the movie in which they
were found. It does not render a STORM image, for that see
sa_utilities/hdf5_to_image.py.

To use the program you will need PyQt6, which you can install with
"pip install PyQt6", or with "pip install .[gui]" when you install
this project from source.

You can start the program from any directory with

  python -m storm_analysis.visualizer.visualizer

or with "python visualizer.py" in the directory where it is located.
Either form also takes an optional directory to start the file
dialogs in.

It understands .dax, .fits, .spe and .tif format movies, and
storm-analysis HDF5 (.hdf5) or Insight3 (.bin) localization files.
Two localization files can be shown at the same time, the first in
green and the second in red.

The following keys can be used to go through the frames of a movie:
"," - backward one frame.
"." - forward one frame.
"k" - backward 200 frames.
"l" - forward 200 frames.
"Home" - go to the start of the movie.
"End" - go to the end of the movie.
(You may need to click on the movie window first)

The mouse scroll wheel zooms in and out when the pointer is over the
movie, and changes frames when it is over most other parts of the
window.

Clicking on the movie will also give you the fitting results for the
nearest localization to where you clicked in the current frame.

The contrast slider runs from 0 to the smallest power of 2 (at least
256) that covers both the brightest pixel shown so far and the maximum
contrast setting. It grows when you reach a brighter frame or type in
a larger maximum, and starts again when you load a new movie.

"Frame Capture" in the File menu saves the movie view, as it is on
the screen, as a PNG file.
