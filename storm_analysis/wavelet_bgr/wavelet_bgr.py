#!/usr/bin/python
#
# Wavelet based background removal as described in this publication:
# Galloway et al., An Iterative Algorithm for Background Removal in Spectroscopy by Wavelet Transforms, Applied Spectroscopy, 2009
#
# Hazen 01/15
#

import numpy
import pywt

class WaveletBGR(object):

    def __init__(self, wavelet_type = 'db4', padding_mode = 'sp1'):
        """
        Create a Wavelet background remover.

        :param wavelet_type: The type of wavelet to use (see the pywt documentation).
        :type wavelet_type: string.
        :param padding_mode: How to extrapolate off the ends of the data (see the pywt documentation).
        :type padding_mode: string.
        """
        self.coeffs = None
        self.padding_mode = padding_mode
        self.shape = None
        self.wavelet_level = None
        self.wavelet_type = wavelet_type

    def _estimateBG_(self, image, wavelet_level):
        """
        Private function that estimates the background.
        """

        # Check if requested level is different from the saved level.
        if (self.wavelet_level is not None) and (self.wavelet_level != wavelet_level):
            self.coeffs = None
        self.wavelet_level = wavelet_level

        # Check if image size is different from the previous level.
        if (self.shape is not None):
            if (self.shape[0] != image.shape[0]) or (self.shape[1] != image.shape[1]):
                self.coeffs = None
        self.shape = image.shape
                
        # 2D Wavelet Transform.
        coeffs = pywt.wavedec2(image, 
                               self.wavelet_type, 
                               mode = self.padding_mode, 
                               level = wavelet_level)
        
        # Set all the details coefficients to zero (if necessary).
        if self.coeffs is None:
            self.coeffs = []
            self.coeffs.append(coeffs[0])
            for i in range(len(coeffs) - 1):
                temp = []
                for j in range(3):
                    temp.append(numpy.zeros(coeffs[i+1][j].shape))
                self.coeffs.append(temp)
        else:
            self.coeffs[0] = coeffs[0]

        # 2D Inverse Wavelet Transform.
        bg_estimate = pywt.waverec2(self.coeffs,
                                    self.wavelet_type,
                                    mode = self.padding_mode)

        return bg_estimate

    def estimateBG(self, image, iterations, threshold, wavelet_level):
        """
        Estimate the background.

        :param image: The image.
        :type image: 2D numpy array.
        :param iterations: Number of iterations to perform to do the estimate.
        :type iterations: integer.
        :param threshold: Threshold for maximum difference between the foreground and the background.
        :type threshold: float or integer.
        :param wavelet_level: The wavelet level.
        :type wavelet_level: integer.        
        """
        temp = numpy.copy(image)
        for i in range(iterations):
            background = self._estimateBG_(temp, wavelet_level)
            mask = (temp > (background + threshold))
            temp[mask] = background[mask]

        return background

    def removeBG(self, image, iterations, threshold, wavelet_level):
        """
        Estimate the background and subtract it from the image.

        :param image: The image.
        :type image: 2D numpy array.
        :param iterations: Number of iterations to perform to do the estimate.
        :type iterations: integer.
        :param threshold: Threshold for maximum difference between the foreground and the background.
        :type threshold: float or integer.
        :param wavelet_level: The wavelet level.
        :type wavelet_level: integer.        
        """
        return image - self.estimateBG(image, iterations, threshold, wavelet_level)


if (__name__ == "__main__"):

    import sys

    import sa_library.datareader as datareader
    import sa_library.daxwriter as daxwriter

    if (len(sys.argv) < 6):
        print "usage <movie> <wavelet_type> <wavelet_level> <iterations> <threshold> <baseline (optional, 100 default)>"
        exit()

    input_movie = datareader.inferReader(sys.argv[1])
    output_dax = daxwriter.DaxWriter("subtracted.dax", 0, 0)

    iterations = int(sys.argv[4])
    threshold = float(sys.argv[5])
    wavelet_level = int(sys.argv[3])    

    offset = 100.0
    if (len(sys.argv) == 7):
        offset = float(sys.argv[6])

    wbgr = WaveletBGR(wavelet_type = sys.argv[2])

    for i in range(input_movie.filmSize()[2]):

        if((i%10) == 0):
            print "Processing frame", i

        image = input_movie.loadAFrame(i) - offset
        sub = wbgr.removeBG(image,
                            iterations,
                            threshold,
                            wavelet_level)
        output_dax.addFrame(sub + offset)

    output_dax.close()
