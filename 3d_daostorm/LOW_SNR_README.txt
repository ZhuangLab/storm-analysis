
For low signal to noise data from an EMCCD camera you may get better results
using the sCMOS analysis pathway as this uses an image convolution based
approach to identify possible localizations in the image. The trick is
to supply it with a fake calibration file, for example, if your images
are 128x128 then you can create the appropriate calibration file like this:

Python 2.7.3 (default, Apr 10 2012, 23:24:47) [MSC v.1500 64 bit (AMD64)] on win
32
Type "help", "copyright", "credits" or "license" for more information.
>>> import numpy
>>> offset = numpy.zeros((128,128))
>>> gain = numpy.ones((128,128))
>>> numpy.save("calib.npy", [offset, gain, gain])
>>> exit()

