#!/usr/bin/env python
"""
Base class for Compressed Sensing algorithms such as
FISTA and ADMM.

Hazen 11/19
"""
import numpy


class CSAlgorithm(object):

    def __init__(self, **kwds):
        super(CSAlgorithm, self).__init__(**kwds)

        self.image = None
        self.x = None

    def cleanup():
        pass

    def dwlsError(self):
        dd = (self.getAx() - self.image)
        return numpy.sum(dd * dd * self.weights)

    def getAx(self):
        pass

    def getXVector(self):
        return self.x
    
    def iterate(self, a_lambda):
        pass

    def l1Error(self):
        return numpy.sum(numpy.abs(self.getXVector()))

    def l2Error(self):
        return numpy.linalg.norm(self.getAx() - self.image)

    def newImage(self, image, background):
        """
        image - The image to deconvolve, includes the background estimate.
        background - An estimate of the background in the image.
        """
        pass

    def run(self, a_lambda, iterations):
        for i in range(iterations):
            self.iterate(a_lambda)

