
This program lets you visualize the results of localization fitting
using 3D-DAOSTORM (or Insight3). Both programs output localization
lists in the same format, but some of the fields have different
meanings. It overlays the identified localizations on a frame of
a movie. It doesn't not render a STORM image. You will need Insight3
or equivalent for that purpose. Insight3 is available by request
from the Zhuang Lab.

To use the program you will need to have PyQt4 installed on your
computer.

You can run the program by typing "python visualizer.py"
in the directory where the program is located (and assuming that
python is in your path).

At present it should understand .dax and .tif format movies, and 
.bin format localization fitting result files.

The following keys can be used to go through the frames of a movie:
"," - backward one frame.
"." - forward one frame.
"k" - backward 200 frames.
"l" - forward 200 frames.
"Home" - go to the start of the movie.
"End" - go to the end of the movie.
(You may need to click on the movie window first)

You can zoom in and out on the movie by clicking on it & then using
the mouse scroll wheel.

Clicking on the movie will also give you the fitting results for the
nearest localization to where you clicked in the current frame.
