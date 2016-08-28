
To analyze low SNR data you have the option to convolve the image
with a 2D gaussian to improve the peak finding step. You can do
this by adding a "filter_sigma" parameter to your analysis xml
file. Please see tests/test_3d_2d_fixed_low_snr.xml for an example
and an explanation of this parameter.

Note that peak fitting is still done on the original non-smoothed
image.
