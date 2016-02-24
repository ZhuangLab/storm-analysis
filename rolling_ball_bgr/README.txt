
Based on:
Stanley Sternberg, "Biomedical Image Processing", IEEE Computer, January 1983.

The rolling ball background removal code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/

This code can be used as a standalone background removal program as well as
a module that you can import.

A sample run in standalone mode:
python ./rolling_ball.py /path/to/movie ball_radius smoothing_sigma

ball_radius - The ball radius to use (in pixels), this should be something
   like 5x the feature size.
smoothing_sigma - Sigma of the gaussian to use for image smoothing prior
   to rolling ball background estimation.

This will create a file called subtracted.dax which contains the original movie 
minus the estimated background with a 100 offset.
