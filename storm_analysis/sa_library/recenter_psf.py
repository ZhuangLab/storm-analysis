#!/usr/bin/env python
"""
Recenter a PSF for FFT image convolution.

Hazen 3/16
"""
import numpy

def recenterPSF(psf, use_roll = True):
    """
    Conceptually anyway it is easier to draw the PSF in the center of a 
    array, but this does not work well when combined with FFT convolution.
    """

    # At some point the we'd remove one of these as they both do the
    # same thing.
    #
    if use_roll:
        recentered = numpy.roll(psf, int(psf.shape[0]/2), axis = 0)
        recentered = numpy.roll(recentered, int(psf.shape[1]/2), axis = 1)

    else:
        shape = psf.shape
        recentered = numpy.zeros_like(psf)
        
        # move ul to br
        recentered[0:int(shape[0]/2),0:int(shape[1]/2)] = psf[int(shape[0]/2):,int(shape[1]/2):]
        # move br to ul
        recentered[int(shape[0]/2):,int(shape[1]/2):] = psf[0:int(shape[0]/2),0:int(shape[1]/2)]
        # move bl to ur
        recentered[0:int(shape[0]/2),int(shape[1]/2):] = psf[int(shape[0]/2):,0:int(shape[1]/2)]
        # move ur to bl
        recentered[int(shape[0]/2):,0:int(shape[1]/2)] = psf[0:int(shape[0]/2),int(shape[1]/2):]

    return recentered


if (__name__ == "__main__"):

    # Check that both algorithms give the same answer.
    #
    psf = numpy.random.uniform(size = (20,10))
    psf1 = recenterPSF(psf)
    psf2 = recenterPSF(psf, use_roll = False)
    print(numpy.allclose(psf1, psf2))


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
