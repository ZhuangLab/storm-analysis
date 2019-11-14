#!/usr/bin/env python
"""
Wavelet based background removal as described in this publication:

Galloway et al., An Iterative Algorithm for Background Removal in Spectroscopy by Wavelet Transforms, Applied Spectroscopy, 2009

Hazen 01/15
"""

import numpy
import pywt

import storm_analysis.sa_library.datareader as datareader
import storm_analysis.sa_library.datawriter as datawriter


class WaveletBGR(object):

    def __init__(self, wavelet_type = 'db4', padding_mode = 'smooth', **kwds):
        """
        Create a Wavelet background remover.

        :param wavelet_type: The type of wavelet to use (see the pywt documentation).
        :type wavelet_type: string.
        :param padding_mode: How to extrapolate off the ends of the data (see the pywt documentation).
        :type padding_mode: string.
        """
        super(WaveletBGR, self).__init__(**kwds)
        
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


class WaveletBGRStdAna(WaveletBGR):
    """
    WaveletBGR object that more easily plugs into a
    standard peak finding object.
    """
    def __init__(self, iterations = None, threshold = None, wavelet_level = None, **kwds):
        super(WaveletBGRStdAna, self).__init__(**kwds)
        
        self.iterations = iterations
        self.threshold = threshold
        self.wavelet_level = wavelet_level

    def estimateBG(self, image):
        return super(WaveletBGRStdAna, self).estimateBG(image,
                                                        self.iterations,
                                                        self.threshold,
                                                        self.wavelet_level)
        

def waveletBGRSub(movie_in, movie_out, wavelet_type, wavelet_level, iterations, threshold, offset = 100):

    input_movie = datareader.inferReader(movie_in)
    output_dax = datawriter.inferWriter(movie_out)

    wbgr = WaveletBGR(wavelet_type = wavelet_type)

    for i in range(input_movie.filmSize()[2]):

        if((i%10) == 0):
            print("Processing frame", i)

        image = input_movie.loadAFrame(i) - offset
        sub = wbgr.removeBG(image,
                            iterations,
                            threshold,
                            wavelet_level)
        output_dax.addFrame(sub + offset)

    output_dax.close()


if (__name__ == "__main__"):

    import argparse
    
    parser = argparse.ArgumentParser(description = 'Wavelet background reduction following Galloway, Applied Spectroscopy, 2009')

    parser.add_argument('--movie_in', dest='movie_in', type=str, required=True,
                        help = "The name of the movie to analyze, can be .dax, .tiff or .spe format.")
    parser.add_argument('--movie_out', dest='movie_out', type=str, required=True,
                        help = "The name of the movie to save the results. This will always be .dax format.")
    parser.add_argument('--wavelet_type', dest='wavelet_type', type=str, required=True,
                        help = "See the pywt documentation, typically something like 'db4'.")
    parser.add_argument('--wavelet_level', dest='wavelet_level', type=int, required=True,
                        help = "How many levels of wavelet decomposition to perform. The larger the number the less response to local changes in the background, usually something like 2.")
    parser.add_argument('--iterations', dest='iterations', type=int, required=True,
                        help = "The number of iterations of background estimation and foreground replacement to perform (see the Galloway paper), usually something like 2.")
    parser.add_argument('--threshold', dest='threshold', type=int, required=True,
                        help = "This should probably be something like 1x to 2x the estimated noise in the background.")
    parser.add_argument('--baseline', dest='baseline', type=bool, required=False, default=100,
                        help = "Camera baseline in ADU.")

    args = parser.parse_args()

    waveletBGRSub(args.movie_in,
                  args.movie_out,
                  args.wavelet_type,
                  args.wavelet_level,
                  args.iterations,
                  args.threshold,
                  args.baseline)
    


