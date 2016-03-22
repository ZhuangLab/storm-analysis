#!/usr/bin/env python
#
# Pure Python code for doing FISTA in 3D.
#
# Hazen 1/16
#

import math
import numpy

numpy.set_printoptions(precision = 4)

# Conceptually anyway it is easier to draw the PSF in the center of a 
# array, but this does not work well when combined with FFT convolution.
#
# FIXME: use sa_library.recenter_psf
#
def recenterPSF(psf):

    shape = psf.shape
    recentered = numpy.zeros(shape)
    # move ul to br
    recentered[0:shape[0]/2,0:shape[1]/2] = psf[shape[0]/2:,shape[1]/2:]
    # move br to ul
    recentered[shape[0]/2:,shape[1]/2:] = psf[0:shape[0]/2,0:shape[1]/2]
    # move bl to ur
    recentered[0:shape[0]/2,shape[1]/2:] = psf[shape[0]/2:,0:shape[1]/2]
    # move ur to bl
    recentered[shape[0]/2:,0:shape[1]/2] = psf[0:shape[0]/2,shape[1]/2:]

    return recentered

class FISTA(object):

    #
    # psfs is an array of psfs for different image planes (nx, ny, nz).
    # They must be the same size as the image that will get analyzed.
    #
    def __init__(self, psfs, timestep):

        if (psfs.shape[0] != psfs.shape[1]):
            print "The PSF must be square (in X-Y)!"
            exit()
        self.shape = psfs.shape

        self.a_mats = psfs
        self.a_mats_fft = []
        self.a_mats_transpose_fft = []
        self.image = None
        self.image_fft = None
        self.l_term = None
        self.nz = self.shape[2]
        self.t_step = timestep
        self.tk = None
        self.x = None
        self.y = None

        # Compute FFTs of PSFs.
        for i in range(self.nz):
            psf = recenterPSF(psfs[:,:,i])
            psf_fft = numpy.fft.fft2(psf)
            self.a_mats_fft.append(psf_fft)
            self.a_mats_transpose_fft.append(numpy.conj(psf_fft))

    def calcAx(self):
        Ax_fft = numpy.zeros((self.shape[0], self.shape[1]), dtype = numpy.complex)
        for i in range(self.nz):
            x_fft = numpy.fft.fft2(self.x[:,:,i])
            Ax_fft += self.a_mats_fft[i] * x_fft
        Ax = numpy.real(numpy.fft.ifft2(Ax_fft))
        return Ax

    def getXVector(self):
        return self.x
    
    def iterate(self, l_term, fista = True):
        self.l_term = l_term
        if fista:
            last_x = self.x
            self.x = self.plk(self.y)
            tk_p1 = 0.5*(1.0 + math.sqrt(1.0 + 4.0 * self.tk * self.tk))
            self.y = self.x + ((self.tk - 1.0)/tk_p1)*(self.x - last_x)
            self.tk = tk_p1
        else:
            self.x = self.plk(self.x)

    def l1Error(self):
        return numpy.sum(numpy.abs(self.x))

    def l2Error(self):
        return numpy.linalg.norm(self.calcAx() - self.image)

    def overallError(self):
        l2_error = self.l2Error()
        return l2_error * l2_error + self.l_term * self.l1Error()

    # New image
    def newImage(self, image):
        self.image = image.copy()
        self.image_fft = numpy.fft.fft2(self.image)
        self.tk = 1.0
        self.x = numpy.zeros(self.shape)
        self.y = numpy.zeros(self.shape)

    def plk(self, vec):

        # Ax is (n,n)
        D = numpy.zeros((self.shape[0], self.shape[1]), dtype = numpy.complex)
        for i in range(self.nz):
            vec_fft = numpy.fft.fft2(vec[:,:,i])
            D += self.a_mats_fft[i]*vec_fft
        D -= self.image_fft

        # At(Ax - b) is (n,n,z)
        Y = numpy.zeros((self.shape))
        for i in range(self.nz):
            Y[:,:,i] = vec[:,:,i] - 2.0 * self.t_step * numpy.real(numpy.fft.ifft2(self.a_mats_transpose_fft[i] * D))

        return self.shrink(Y)

    def shrink(self, vec):
        t1 = numpy.abs(vec) - self.l_term * self.t_step
        mask = (t1 <= 0.0)
        t1[mask] = 0.0
        
        return t1 * numpy.sign(vec)

    
#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
