#!/usr/bin/env python
"""
Classes and functions for image correlation.

Hazen 07/09
"""
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import scipy
import scipy.signal

import storm_analysis.sa_library.gaussfit as gaussfit


class Align3D(object):
    """
    This class is used to aligns two 3D image stacks (translation
    only). It makes extensive use of FFTs for image shifting and 
    calculation of derivatives. It optimizes the sum of the product 
    of the two images.
    """
    def __init__(self, ref_image, xy_margin = 1, z_margin = 1):
        """
        ref_image - A 3D image, with z as the first axis.
        xy_margin - Amount of zero padding to put around the image
                    in pixels so that when it is shifted it does
                    not wrap.
        z_margin - Same as xy_margin.
        """
        assert(len(ref_image.shape) == 3), "Image must be a 3D stack."

        self.random_corr = None  # This is the (average) expected correlation for two random images.
        self.random_dev = None   # This is the expected standard deviation of the correlation between two random images.
        self.other_image = None
        self.other_image_fft = None
        
        self.x_size = ref_image.shape[1]
        self.y_size = ref_image.shape[2]
        self.z_size = ref_image.shape[0]
        
        self.xy_margin = xy_margin
        self.z_margin = z_margin

        self.ref_image = self.padImage(ref_image)

        # Coefficents for FFT shifting and derivative calculation.
        kx = numpy.fft.fftfreq(self.x_size + 2*self.xy_margin)
        ky = numpy.fft.fftfreq(self.y_size + 2*self.xy_margin)
        kz = numpy.fft.fftfreq(self.z_size + 2*self.z_margin)

        [self.kz, self.kx, self.ky] = numpy.meshgrid(kz, kx, ky, indexing = 'ij')

    def dfn_dx(self, dx, dy, dz):
        # Translate.
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))

        # Take derivative.
        temp_fft = -1j * 2.0 * numpy.pi * self.kx * temp_fft

        temp = numpy.sum(self.ref_image * numpy.fft.ifftn(temp_fft))
        return numpy.real(temp)

    def dfn_dy(self, dx, dy, dz):
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))
        temp_fft = -1j * 2.0 * numpy.pi * self.ky * temp_fft

        temp = numpy.sum(self.ref_image * numpy.fft.ifftn(temp_fft))
        return numpy.real(temp)

    def dfn_dz(self, dx, dy, dz):
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))
        temp_fft = -1j * 2.0 * numpy.pi * self.kz * temp_fft

        temp = numpy.sum(self.ref_image * numpy.fft.ifftn(temp_fft))
        return numpy.real(temp)
        
    def fn(self, dx, dy, dz):
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))
        temp = numpy.abs(numpy.fft.ifftn(temp_fft))
        return numpy.sum(temp * self.ref_image)

    def padImage(self, image):
        """
        Pad the image to size with margins, edge values are filled in by replication.
        """
        assert(len(image.shape) == 3), "Image must be a 3D stack."
        assert(image.shape[0] == self.z_size), "Image must be the same size as the reference image."
        assert(image.shape[1] == self.x_size), "Image must be the same size as the reference image."
        assert(image.shape[2] == self.y_size), "Image must be the same size as the reference image."

        return numpy.pad(image,
                         ((self.z_margin, self.z_margin),
                          (self.xy_margin, self.xy_margin),
                          (self.xy_margin, self.xy_margin)),
                         'edge')

    def setOtherImage(self, image):
        self.other_image = self.padImage(image)
        self.other_image_fft = numpy.fft.fftn(self.other_image)

        # Calculate expected correlation for two random images. Average
        # value in the first image times the average value in the second
        # image times the number of non-zero pixels.
        #
        self.random_corr = numpy.sum(self.ref_image) * numpy.sum(self.other_image) / float(self.ref_image.size)

        # FIXME: Figure out correct value here.
        self.random_dev = 1.0


def absIntRound(num):
    return abs(int(round(num)))

def crop2DImages(image1, image2, dx, dy):
    [size_x, size_y] = image1.shape
    adx = absIntRound(dx)
    ady = absIntRound(dy)
    if dx >= 0 and dy >= 0:
        image1 = image1[:size_x-adx,:size_y-ady]
        image2 = image2[adx:       ,ady:       ]
    if dx >= 0 and dy < 0:
        image1 = image1[:size_x-adx,ady:       ]
        image2 = image2[adx:       ,:size_y-ady]
    if dx < 0 and dy >= 0:
        image1 = image1[adx:       ,:size_y-ady]
        image2 = image2[:size_x-adx,ady:       ]
    if dx < 0 and dy < 0:
        image1 = image1[adx:       ,ady:       ]
        image2 = image2[:size_x-adx,:size_y-ady]
    return [image1, image2]

def crop3DImages(image1, image2, dx, dy):
    [size_x, size_y, size_z] = image1.shape
    adx = absIntRound(dx)
    ady = absIntRound(dy)
    if dx >= 0 and dy >= 0:
        image1 = image1[:size_x-adx,:size_y-ady,:]
        image2 = image2[adx:       ,ady:       ,:]
    if dx >= 0 and dy < 0:
        image1 = image1[:size_x-adx,ady:       ,:]
        image2 = image2[adx:       ,:size_y-ady,:]
    if dx < 0 and dy >= 0:
        image1 = image1[adx:       ,:size_y-ady,:]
        image2 = image2[:size_x-adx,ady:       ,:]
    if dx < 0 and dy < 0:
        image1 = image1[adx:       ,ady:       ,:]
        image2 = image2[:size_x-adx,:size_y-ady,:]
    return [image1, image2]

def xyCorrelate(image1, image2):
    return scipy.signal.fftconvolve(image1, image2[::-1, ::-1], mode="same")

def xyOffset(image1, image2, scale, center = None):
    """
    Note that the search is limited to a X by X region
    in the center of the overlap between the two images.
    """
    image1 = image1 - numpy.median(image1)
    image2 = image2 - numpy.median(image2)

    result = xyCorrelate(image1, image2)
    
    if False:
        import tifffile
        tifffile.imsave("corr_image1.tif", image1.astype(numpy.float32))
        tifffile.imsave("corr_image2.tif", image2.astype(numpy.float32))
        tifffile.imsave("corr_result.tif", result.astype(numpy.float32))

    # These are the coordinates of the image center.
    mx = int(round(0.5 * result.shape[0]))
    my = int(round(0.5 * result.shape[1]))

    # This is the area to search, 30 pixels * scale.
    s_size = int(30 * int(scale))
    sx = s_size
    sy = s_size

    # Adjust if the image is really small.
    if mx < (s_size + 5):
        sx = mx - 5
    if my < (s_size + 5):
        sy = my - 5

    # Use center position provided by the user.
    if isinstance(center, list):
        rx = int(round(center[0]) + sx)
        ry = int(round(center[1]) + sy)
        if False:
            print(rx)
            print(numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 1)))
            print(ry)
            print(numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 0)))

    # Otherwise find local maximum near the image center.
    else:
        rx = numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 1))
        ry = numpy.argmax(numpy.max(result[mx-sx:mx+sx+1,my-sy:my+sy+1], axis = 0))

    # Adjust from maxima search area to full image size.
    rx += (mx - sx)
    ry += (my - sy)

    #
    # Fit a gaussian to the maxima.
    #
    # FIXME: Should use elliptical gaussian with variable axis?
    #
    # FIXME: Should fit a larger / smaller area?
    #
    [fit, success] = gaussfit.fitSymmetricGaussian(result[rx-5:rx+6,ry-5:ry+6], 2.0)

    if 0:
        fit_fn = gaussfit.symmetricGaussian(*fit)
        surf = numpy.zeros((11,11))
        fp = open("fit.txt", "w")
        fp.write("x, y, correlation, fit\n")
        for x in range(11):
            for y in range(11):
                surf[x,y] = fit_fn(x,y)
                fp.write(str(x) + ", " + str(y) + ", " + str(result[rx-5+x,ry-5+y]) + ", " + str(surf[x,y]) + "\n")
        fp.close()

    if (fit[0] == fit[1]) or (fit[2] < 2.0) or (fit[2] > 8.0) or (fit[3] < 2.0) or (fit[3] > 8.0):
        print("Bad fit center:", fit[2], fit[3])
        success = False
    return [result[mx-sx:mx+sx+1,my-sy:my+sy+1],
            rx + fit[2] - 5.0 - 0.5 * result.shape[0],
            ry + fit[3] - 5.0 - 0.5 * result.shape[1],
            success]

def xyOffsetWithDxDy(image1, image2, dx, dy, scale):
    [cimage1, cimage2] = crop2DImages(image1, image2, dx, dy)
    return xyOffset(cimage1, cimage2, scale)

def zOffset(image1, image2):

    # FIXME: It looks like this pretty much always returns True for fitting success
    #        even if there is no correlation between the two images. A more reliable
    #        metric is needed..
    #
    image1 = image1 - numpy.median(image1)
    image2 = image2 - numpy.median(image2)
    [size_x, size_y, size_z] = image1.shape
    corr = numpy.zeros(2*size_z-1)
    for i in range(size_z-1):
        corr[i] = numpy.sum(image1[:,:,size_z-1-i:size_z] * image2[:,:,:i+1])/float(i+1)
    for i in range(size_z):
        corr[size_z-1+i] = numpy.sum(image1[:,:,:size_z-i] * image2[:,:,i:size_z])/float(size_z-i)
    
    # This handles data that is actually 2D
    number_non_zero = 0
    which_non_zero = 0
    for i in range(2*size_z-1):
        if (corr[i] > 0.0):
            which_non_zero = i
            number_non_zero += 1

    if (number_non_zero > 1):
        [fit, success] = gaussfit.fitLorentzian(corr)
        lorentzian_fit = gaussfit.lorentzian(*fit)(*numpy.indices(corr.shape))
    else:
        fit = [0, 0, which_non_zero]
        success = True
        lorentzian_fit = corr

    # Debugging.
    if False:
        x = numpy.arange(corr.size)
        fig = pyplot.figure()
        pyplot.scatter(x, corr, color = "magenta")
        pyplot.plot(x, lorentzian_fit, color = "green")
        pyplot.show()
        
    return [corr, lorentzian_fit, fit[2] - size_z + 1, success]

def zOffsetWithDxDy(image1, image2, dx, dy):
    [cimage1, cimage2] = crop3DImages(image1, image2, dx, dy)
    return zOffset(cimage1, cimage2)



# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
