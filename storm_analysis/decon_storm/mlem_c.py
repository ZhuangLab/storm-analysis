#!/usr/bin/python
#
# Python front end for mlem_sparse.so
#
# Hazen 01/11
#
# 
# Tweaked to work with the new mlem_sparse struct
# based approach.
#
# Hazen 03/11
#


from ctypes import *
import numpy
import os
import scipy
import scipy.ndimage.filters
import sys

from numpy.ctypeslib import ndpointer, as_ctypes

# define useful pointers
c_int_p = POINTER(c_int)
c_double_p = POINTER(c_double)

# Load mlem C library
directory = os.path.dirname(__file__)
if(os.path.exists(directory + "/mlem_sparse.so")):
    mlem = cdll.LoadLibrary(directory + "/mlem_sparse.so")
else:
    mlem = cdll.LoadLibrary(directory + "/mlem_sparse.dll")

# Define structures
class GAUSS(Structure):
    _fields_ = [("half_x_size", c_int),
                ("half_y_size", c_int),
                ("x_size", c_int),
                ("y_size", c_int),
                ("xc", c_double),
                ("yc", c_double),
                ("zc", c_double),
                ("width_x", c_double),
                ("width_y", c_double),
                ("coeff", c_double_p)]

class PSF(Structure):
    _fields_ = [("p_num_i", c_int),
                ("p_loc", c_int),
                ("p_xj", c_double),
                ("p_comp", c_double),
                ("p_i", c_int_p),
                ("p_coeff", c_double_p)]

# Define functions
mlem.backward.argtypes = [c_void_p]
mlem.backwardCompressed.argtypes = [c_void_p, c_double]
mlem.backwardCompressedFixedBg.argtypes = [c_void_p, c_double]
mlem.backwardVarCompressionFixedBg.argtypes = [c_void_p]
mlem.cleanup.argtypes = [c_void_p]
mlem.cull.argtypes = [c_void_p, c_double]
mlem.forward.argtypes = [c_void_p]
mlem.fractionLow.argtypes = [c_void_p, c_double]
mlem.fractionLow.restype = c_double
mlem.getAPsf.argtypes = [c_void_p, c_void_p, c_int]
mlem.getBackground.argtypes = [c_void_p, ndpointer(dtype=numpy.float64)]
mlem.getDiff.argtypes = [c_void_p]
mlem.getDiff.restype = c_double
mlem.getFit.argtypes = [c_void_p, ndpointer(dtype=numpy.float64)]
mlem.getForeground.argtypes = [c_void_p, ndpointer(dtype=numpy.float64)]
mlem.getGauss.argtypes = [c_void_p, c_int, c_int, c_int]
mlem.getPeaks.argtypes = [c_void_p, c_void_p, c_double, c_int, c_int]
mlem.getPeaks.restype = c_void_p
mlem.newBackground.argtypes = [c_void_p, ndpointer(dtype=numpy.float64), c_int]
mlem.newForeground.argtypes = [c_void_p, ndpointer(dtype=numpy.float64), c_double]
mlem.setCompression.argtypes = [c_void_p, ndpointer(dtype=numpy.float64)]
#mlem.setForeground.argtypes = [c_void_p, ndpointer(dtype=numpy.float64)]
mlem.setup2D.argtypes = [c_double, c_int, c_int]
mlem.setup2D.restype = c_void_p

# Helper functions
def getArray(fn, decon, size):
    arr = numpy.ascontiguousarray(numpy.zeros((size,size)))
    fn(decon, arr)
    return arr

# FIXME: This appears to be almost the same as subBackground
def getBackground(image, gridsize = 8):
    max_i = image.shape[0]/gridsize
    background = numpy.zeros(image.shape)
    for i in range(max_i):
        for j in range(max_i):
            tmp = image[i*gridsize:(i+1)*gridsize,
                        j*gridsize:(j+1)*gridsize]
            med_ij = numpy.median(tmp)
            mean_ij = numpy.mean(tmp)
#            background[i*gridsize:(i+1)*gridsize,
#                       j*gridsize:(j+1)*gridsize] = med_ij - (mean_ij - med_ij)
            background[i*gridsize:(i+1)*gridsize,
                       j*gridsize:(j+1)*gridsize] = med_ij - (mean_ij - med_ij)
    background = scipy.ndimage.filters.uniform_filter(background,
                                                      size = 2*gridsize)
    return background

def padImage(image, pad_size):
    [x_size, y_size] = image.shape
    lg_image = numpy.ones((x_size+2*pad_size,y_size+2*pad_size))
    lg_image[pad_size:(x_size+pad_size),pad_size:(y_size+pad_size)] = image
    lg_image[0:pad_size,:] = numpy.flipud(lg_image[pad_size:2*pad_size,:])
    lg_image[(x_size+pad_size):(x_size+2*pad_size),:] = numpy.flipud(lg_image[x_size:(x_size+pad_size),:])
    lg_image[:,0:pad_size] = numpy.fliplr(lg_image[:,pad_size:2*pad_size])
    lg_image[:,(y_size+pad_size):(y_size+2*pad_size)] = numpy.fliplr(lg_image[:,y_size:(y_size+pad_size)])
    return lg_image

def subBackground(image, gridsize = 16):
    background = getBackground(image, gridsize = gridsize)
    return image - background

# MLEM Deconvolution class
class Fitter():
    def __init__(self, background, sigma, scale, threshold, tol = 1.0e-4, gridsize = False):
        self.compression = None
        self.cull_threshold = 0.1
        self.foreground = None
        self.image = None
        self.init = True
        self.pad_size = 10
        self.scale = scale
        self.size = background.shape[0]
        self.threshold = threshold
        self.tol = tol

        self.decon_size = (self.size + 2*self.pad_size) * self.scale
        if gridsize and (not ((self.size % gridsize)==1)):
            print "warning: size modulo", gridsize, "must equal 1."
        else:
            gridsize = self.size + 2*self.pad_size - 1
        self.decon = mlem.setup2D(float(sigma), int(self.size + 2*self.pad_size), int(scale))
        self.newBackground(padImage(background, self.pad_size), gridsize)

    def cleanup(self):
        mlem.cleanup(self.decon)

    def cull(self, threshold):
        mlem.cull(self.decon, float(threshold))

    def __doFit__(self, backward_fn, verbose, max_iter):
        last_diff = None
        for i in range(max_iter):
            mlem.forward(self.decon)
            backward_fn()
            #self.cull(self.cull_threshold)

            diff = mlem.getDiff(self.decon)
            last_diff = diff
            if verbose:
                if((i%10)==0):
                    print i,diff
        return diff

    def doFit(self, verbose = False, max_iter = 1000):
        backward_fn = lambda : mlem.backward(self.decon)
        return self.__doFit__(backward_fn, verbose, max_iter)

    def doFitCompressed(self, compression, verbose = False, max_iter = 1000):
        backward_fn = lambda : mlem.backwardCompressed(self.decon, float(compression))
        return self.__doFit__(backward_fn, verbose, max_iter)

    def doFitCompressedFixedBg(self, compression, verbose = False, max_iter = 1000):
        backward_fn = lambda : mlem.backwardCompressedFixedBg(self.decon, float(compression))
        return self.__doFit__(backward_fn, verbose, max_iter)

    def doFitVarCompressionFixedBg(self, verbose = False, max_iter = 1000):
        backward_fn = lambda : mlem.backwardVarCompressionFixedBg(self.decon)
        return self.__doFit__(backward_fn, verbose, max_iter)

    def getAPsf(self, index):
        i_arr = numpy.ascontiguousarray(numpy.zeros((1000)),dtype = numpy.int32)
        c_arr = numpy.ascontiguousarray(numpy.zeros((1000)))
        psf = PSF(p_i = i_arr.ctypes.data_as(c_int_p),
                  p_coeff = c_arr.ctypes.data_as(c_double_p))
        psf_ptr = pointer(psf)
        mlem.getAPsf(self.decon, psf_ptr, index)
        locs = numpy.resize(i_arr, psf.p_num_i)
        coeff = numpy.resize(c_arr, psf.p_num_i)
        return psf, locs, coeff

    def getBackground(self):
        return getArray(mlem.getBackground, self.decon, self.size)

    def getFit(self):
        return getArray(mlem.getFit, self.decon, self.size)

    def getForeground(self, remove_padding = True):
        foreground = getArray(mlem.getForeground, self.decon, self.decon_size)
        if remove_padding:
            margin = self.pad_size * self.scale
            foreground = foreground[margin:-margin,margin:-margin]
        return foreground

    def getGauss(self, gx, gy, gz):
        arr = numpy.ascontiguousarray(numpy.zeros((1000)))
        gauss = GAUSS(coeff = arr.ctypes.data_as(c_double_p))
        g_ptr = pointer(gauss)
        mlem.getGauss(g_ptr, gx, gy, gz)
        coeff = numpy.resize(arr, (gauss.x_size, gauss.y_size))
        return gauss, coeff

    def getPeaks(self, threshold, guard_band, neighborhood):
        n = c_int(0)
        buffer = mlem.getPeaks(self.decon, byref(n),float(threshold),int(guard_band),int(neighborhood))
        buffer = cast(buffer, POINTER(c_double))
        b_size = 4*n.value
        peaks = []
        for i in range(b_size):
            peaks.append(buffer[i])
        return numpy.array(peaks)

    def newBackground(self, image, gridsize):
        self.image = numpy.ascontiguousarray(image).astype(numpy.float64)
        mlem.newBackground(self.decon, self.image, gridsize)

    def newForeground(self, image, threshold = None):
        if not threshold:
            threshold = self.threshold
        self.image = padImage(numpy.ascontiguousarray(image).astype(numpy.float64), self.pad_size)
        mlem.newForeground(self.decon, self.image, float(threshold))
        if 1:
            mlem.forward(self.decon)

    def setCompression(self, compression):
        if not (compression.shape[0] == self.decon_size):
            print "warning: compression size should equal decon size"
        self.compression = numpy.ascontiguousarray(compression).astype(numpy.float64)
        mlem.setCompression(self.decon, self.compression)

    def setForeground(self, high_res):
        if not (high_res.shape[0] == self.decon_size):
            print "warning: new foreground size should equal decon size"
        self.foreground = numpy.ascontiguousarray(high_res).astype(numpy.float64)

    def subBackground(self, image, gridsize = 16):
        max_i = image.shape[0]/gridsize
        print max_i, image.shape[0], gridsize
        background = numpy.zeros(image.shape)
        for i in range(max_i):
            for j in range(max_i):
                tmp = image[i*gridsize:(i+1)*gridsize,
                            j*gridsize:(j+1)*gridsize]
                med_ij = numpy.median(tmp)
                print i, j, med_ij
                mean_ij = numpy.mean(tmp)
                background[i*gridsize:(i+1)*gridsize,
                           j*gridsize:(j+1)*gridsize] = med_ij
        background = scipy.ndimage.filters.uniform_filter(background,
                                                          size = gridsize)
        print numpy.min(background), numpy.max(background), background.shape
        return image - background


    
