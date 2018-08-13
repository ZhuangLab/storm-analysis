#!/usr/bin/env python
"""
Classes and functions for image correlation.

Hazen 07/09
"""
import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import scipy
import scipy.optimize
import scipy.signal

import storm_analysis.sa_library.gaussfit as gaussfit


class ImageCorrelationException(Exception):
    pass


class Align3D(object):
    """
    This class is used to aligns two 3D image stacks (translation
    only). It makes extensive use of FFTs for image shifting and 
    calculation of derivatives. 
    """
    def __init__(self, ref_image, xy_margin = 1, z_margin = 1):
        """
        ref_image - A 3D image, with z as the first axis.
        xy_margin - Amount of zero padding to put around the image
                    in pixels so that when it is shifted it does
                    not wrap.
        z_margin - Same as xy_margin.
        """
        super(Align3D, self).__init__()
        
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

    def align(self, dx = 0.0, dy = 0.0, dz = 0.0):
        """
        Find the optimal alignment for 'other' and return a translated
        version of other, along with an estimate of the alignment quality.
        """
        if self.other_image is None:
            raise ImageCorrelationException("Other image not specified.")
        
        [disp, success, fun, status] = self.maximize(dx = dx, dy = dy, dz = dz)
        if not success:

            # Ignore status = 2 precision loss messages for now as it looks
            # like the optimization is just complaining that it can't refine
            # to a precision that we did not care about anyway. Presumably
            # there is some way to set the desired precision, but using a
            # larger value for 'xtol' does not seem to be it.
            #
            if (status == 2):
                print("Warning! Possible precision loss in alignment")
            else:
                raise ImageCorrelationException("Align3d.maximize failed.")

        temp = self.translate(disp[0], disp[1], disp[2])

        # Trim off padding.
        temp = temp[self.z_margin:-self.z_margin,
                    self.xy_margin:-self.xy_margin,
                    self.xy_margin:-self.xy_margin]
        
        # Quality score is how many sigma we are away from the expected correlation
        # between two Poisson random images.
        #
        q_score = (fun - self.random_corr)/self.random_dev

        return [temp, q_score, disp]

    def dx(self, dx, dy, dz, order = 1):
        # Translate.
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))

        # Take derivative.
        dd = numpy.power(-1j * 2.0 * numpy.pi * self.kx, order)
        temp_fft = dd * temp_fft

        # Return real part of inverse FFT. The complex part should be very small anyway..
        return numpy.real(numpy.fft.ifftn(temp_fft))

    def dy(self, dx, dy, dz, order = 1):
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))
        dd = numpy.power(-1j * 2.0 * numpy.pi * self.ky, order)
        temp_fft = dd * temp_fft
        return numpy.real(numpy.fft.ifftn(temp_fft))

    def dz(self, dx, dy, dz, order = 1):
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))
        dd = numpy.power(-1j * 2.0 * numpy.pi * self.kz, order)
        temp_fft = dd * temp_fft
        return numpy.real(numpy.fft.ifftn(temp_fft))
        
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

    def translate(self, dx, dy, dz):
        temp_fft = self.other_image_fft * numpy.exp(-1j * 2.0 * numpy.pi * (self.kx * dx + self.ky * dy + self.kz * dz))
        return numpy.abs(numpy.fft.ifftn(temp_fft))


class Align2D(Align3D):
    """
    This seemed like the simplest approach. We just add an extra dimension
    to the 2D array so that it works as a 3D array.
    """
    def __init__(self, ref_image, xy_margin = 1):
        """
        ref_image - A 3D image, with z as the first axis.
        xy_margin - Amount of zero padding to put around the image
                    in pixels so that when it is shifted it does
                    not wrap.
        z_margin - Same as xy_margin.
        """
        ref_image = numpy.expand_dims(ref_image, axis = 0)
        super(Align2D, self).__init__(ref_image, xy_margin = xy_margin, z_margin = 0)
        
    def setOtherImage(self, image):
        image = numpy.expand_dims(image, axis = 0)
        super(Align2D, self).setOtherImage(image)
        

class Align3DProduct(Align3D):
    """
    This class aligns the image based on the sum of the product of the
    two images.
    """
    def func(self, x, sign = 1.0):
        """
        Calculation function (scipy.optimize.minimize friendly form).
        """
        return sign * numpy.sum(self.ref_image * self.translate(x[0], x[1], x[2]))
        
    def hessian(self, x, sign = 1.0):
        """
        Calculation hessian (scipy.optimize.minimize friendly form).
        """
        hess = numpy.zeros((3,3))

        # Off diagonal terms.
        dd = [self.dx(x[0],x[1],x[2]),
              self.dy(x[0],x[1],x[2]),
              self.dz(x[0],x[1],x[2])]
        for i in range(3):
            for j in range(3):
                if (i != j):
                    hess[i,j] = -sign * numpy.sum(self.ref_image * dd[i] * dd[j])

        # Diagonal terms.
        hess[0,0] = sign * numpy.sum(self.ref_image * self.dx(x[0], x[1], x[2], order = 2))
        hess[1,1] = sign * numpy.sum(self.ref_image * self.dy(x[0], x[1], x[2], order = 2))
        hess[2,2] = sign * numpy.sum(self.ref_image * self.dz(x[0], x[1], x[2], order = 2))
        
        return hess

    def jacobian(self, x, sign = 1.0):
        """
        Calculation jacobian (scipy.optimize.minimize friendly form).
        """
        dfn_dx = numpy.sum(self.ref_image * self.dx(x[0], x[1], x[2]))
        dfn_dy = numpy.sum(self.ref_image * self.dy(x[0], x[1], x[2]))
        dfn_dz = numpy.sum(self.ref_image * self.dz(x[0], x[1], x[2]))
        return sign * numpy.array([dfn_dx, dfn_dy, dfn_dz])

    def setOtherImage(self, image):
        super(Align3DProduct, self).setOtherImage(image)

        # Calculate expected product for two Poisson random images. 
        #
        
        # Poisson lambda is the average pixel occupancy.
        #
        lambda_img1 = numpy.sum(self.ref_image) / float(self.ref_image.size)
        lambda_img2 = numpy.sum(self.other_image) / float(self.ref_image.size)

        # Mean product is the product of average pixel occupancy times the number of pixels.
        #
        self.random_corr = lambda_img1 * lambda_img2 * float(self.ref_image.size)

        # Variance per pixel is:
        #
        # Var(im1 x im2)) = E(im1^2 * im2^2) - E(im1*im2)^2
        #                 = Var(im1) * Var(im2) + Var(im1) * E(im2)^2 + Var(im2) * E(im1)^2
        #
        # Standard deviation is sqrt(per pixel variance x number of pixels).
        #
        self.random_dev = math.sqrt(lambda_img1 * lambda_img2 * (1.0 + lambda_img1 + lambda_img2) * float(self.ref_image.size))
    

class Align3DProductLM(Align3DProduct):
    """
    Align using a variant of the Levenberg-Marquadt algorithm.
    """
    def __init__(self, ref_image, xy_margin = 1, z_margin = 1, tolerance = 1.0e-6, max_reps = 100):
        super(Align3DProductLM, self).__init__(ref_image, xy_margin = xy_margin, z_margin = z_margin)

        self.fn_curr = None
        self.fn_old = None
        self.lm_lambda = 1.0
        self.max_reps = max_reps
        self.tolerance = tolerance
        
    def hasConverged(self):
        return (abs((self.fn_old - self.fn_curr)/self.fn_curr) < self.tolerance)

    def update(self, x):
        """
        Return the update vector at x.
        """
        jac = self.jacobian(x, sign = -1.0)
        hess = self.hessian(x, sign = -1.0)
        for i in range(jac.size):
            hess[i,i] += hess[i,i] * self.lm_lambda
        delta = numpy.linalg.solve(hess, jac)
        return delta

    def maximize(self, dx = 0.0, dy = 0.0, dz = 0.0):
        self.lm_lambda = 1.0
        xo = numpy.array([dx, dy, dz])

        self.fn_curr = self.func(xo, sign = -1.0)
        for i in range(self.max_reps):
            xn = xo - self.update(xo)
            fn = self.func(xn, sign = -1.0)

            # If we did not improve increase lambda and try again.
            if (fn > self.fn_curr):
                self.lm_lambda = 2.0 * self.lm_lambda
                continue

            self.lm_lambda = 0.9 * self.lm_lambda
            self.fn_old = self.fn_curr
            self.fn_curr = fn
            xo = xn
                
            if self.hasConverged():
                break

        success = (i < (self.max_reps - 1))
        return [xo, success, -self.fn_curr, 0]
        
    
class Align3DProductNewtonCG(Align3DProduct):
    """
    Align using the 'Newton-CG' algorithm in scipy.optimize.minimize.
    """

    def maximize(self, dx = 0.0, dy = 0.0, dz = 0.0):
        """
        Find the offset that optimizes the correlation of self and other.
        """
        x0 = numpy.array([dx, dy, dz])

        fit = scipy.optimize.minimize(self.func,
                                      x0,
                                      args=(-1.0,),
                                      method='Newton-CG',
                                      jac=self.jacobian,
                                      hess=self.hessian,
                                      options={'xtol': 1e-3, 'disp': False})

        if not fit.success:
            print("Maximization failed with:")
            print(fit.message)
            print("Status:", fit.status)
            print("X:", fit.x)
            print("Function value:", -fit.fun)
            print()
                        
        return [fit.x, fit.success, -fit.fun, fit.status]


class Align2DProduct(Align2D):
    """
    This class aligns the image based on the sum of the product of the
    two images.
    """
    def func(self, x, sign = 1.0):
        """
        Calculation function (scipy.optimize.minimize friendly form).
        """
        return sign * numpy.sum(self.ref_image * self.translate(x[0], x[1], 0.0))
        
    def hessian(self, x, sign = 1.0):
        """
        Calculation hessian (scipy.optimize.minimize friendly form).
        """
        hess = numpy.zeros((2,2))

        # Off diagonal terms.
        dd = [self.dx(x[0],x[1], 0.0),
              self.dy(x[0],x[1], 0.0)]
        for i in range(2):
            for j in range(2):
                if (i != j):
                    hess[i,j] = -sign * numpy.sum(self.ref_image * dd[i] * dd[j])

        # Diagonal terms.
        hess[0,0] = sign * numpy.sum(self.ref_image * self.dx(x[0], x[1], 0.0, order = 2))
        hess[1,1] = sign * numpy.sum(self.ref_image * self.dy(x[0], x[1], 0.0, order = 2))
        
        return hess

    def jacobian(self, x, sign = 1.0):
        """
        Calculation jacobian (scipy.optimize.minimize friendly form).
        """
        dfn_dx = numpy.sum(self.ref_image * self.dx(x[0], x[1], 0.0))
        dfn_dy = numpy.sum(self.ref_image * self.dy(x[0], x[1], 0.0))
        return sign * numpy.array([dfn_dx, dfn_dy])

    def setOtherImage(self, image):
        super(Align2DProduct, self).setOtherImage(image)

        # Calculate expected product for two Poisson random images. 
        #
        
        # Poisson lambda is the average pixel occupancy.
        #
        lambda_img1 = numpy.sum(self.ref_image) / float(self.ref_image.size)
        lambda_img2 = numpy.sum(self.other_image) / float(self.ref_image.size)

        # Mean product is the product of average pixel occupancy times the number of pixels.
        #
        self.random_corr = lambda_img1 * lambda_img2 * float(self.ref_image.size)

        # Variance per pixel is:
        #
        # Var(im1 x im2)) = E(im1^2 * im2^2) - E(im1*im2)^2
        #                 = Var(im1) * Var(im2) + Var(im1) * E(im2)^2 + Var(im2) * E(im1)^2
        #
        # Standard deviation is sqrt(per pixel variance x number of pixels).
        #
        self.random_dev = math.sqrt(lambda_img1 * lambda_img2 * (1.0 + lambda_img1 + lambda_img2) * float(self.ref_image.size))


class Align2DProductNewtonCG(Align2DProduct):
    """
    Align using the 'Newton-CG' algorithm in scipy.optimize.minimize.
    """

    def maximize(self, dx = 0.0, dy = 0.0):
        """
        Find the offset that optimizes the correlation of self and other.
        """
        x0 = numpy.array([dx, dy])

        fit = scipy.optimize.minimize(self.func,
                                      x0,
                                      args=(-1.0,),
                                      method='Newton-CG',
                                      jac=self.jacobian,
                                      hess=self.hessian,
                                      options={'xtol': 1e-3, 'disp': False})

        if not fit.success:
            print("Maximization failed with:")
            print(fit.message)
            print("Status:", fit.status)
            print("X:", fit.x)
            print("Function value:", -fit.fun)
            print()
                        
        return [fit.x, fit.success, -fit.fun, fit.status]


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
