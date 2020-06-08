#!/usr/bin/env python
"""
Given a spline file create an object that generate PSFs
at requested z values.

Hazen 01/16
"""

import pickle
import numpy

import storm_analysis.sa_library.fitting as fitting

import storm_analysis.spliner.cubic_spline_c as cubicSplineC
import storm_analysis.spliner.spline2D as spline2D
import storm_analysis.spliner.spline3D as spline3D


class SplineToPSF(fitting.PSFFunction):

    def cleanup(self):
        if self.c_spline is not None:
            self.c_spline.cleanup()
        
    def getCPointer(self):
        return self.c_spline.getCPointer()
        
    def getMargin(self):
        return int(self.getSize()/2 + 2)

    def getSize(self):
        """
        This returns the X/Y size in pixels covered by the spline.
        """
        return self.spline_size
    
    def loadSplineFile(self, spline_file):
        """
        Load the spline_file if it has not already been loaded. Otherwise
        just return it under the assumption that is a unpickled spline file.
        """
        if isinstance(spline_file, str):
            with open(spline_file, 'rb') as fp:
                spline_data = pickle.load(fp)
            return spline_data
        else:
            return spline_file
            
    def upsize(self, psf, shape, up_sample):
        """
        Handles drawing into a large array (if requested).
        """
        if shape is not None:
            psf_size = psf.shape[0]
            im_size_x = shape[0] * up_sample
            im_size_y = shape[1] * up_sample

            if((psf_size%2) == 0):
                start_x = int(im_size_x/2 - psf_size/2)
                start_y = int(im_size_y/2 - psf_size/2)
            else:
                start_x = int(im_size_x/2 - psf_size/2) + 1
                start_y = int(im_size_y/2 - psf_size/2) + 1
                
            end_x = start_x + psf_size
            end_y = start_y + psf_size

            temp = numpy.zeros((im_size_x, im_size_y))
            temp[start_x:end_x,start_y:end_y] = psf

            psf = temp

        return psf


class SplineToPSF2D(SplineToPSF):

    def __init__(self, spline_file = None, **kwds):
        super(SplineToPSF2D, self).__init__(**kwds)
        
        spline_data = self.loadSplineFile(spline_file)

        # Check that this is not an old spline, which will be transposed.
        assert("version" in spline_data), "v0 spline file detected! Please re-measure!"
        assert(spline_data["version"] >= 2.0), "v0/v1 spline file detected! Please re-measure!"
        
        # These are used when we check the starting z-value(s)
        # provided by the user, if any.
        self.zmin = -1.0
        self.zmax = 1.0

        # The Python representation of the spline.
        self.spline = spline2D.Spline2D(spline_data["spline"], spline_data["coeff"])
        self.spline_size = self.spline.getSize()

        # The C representation of the spline. This class does not use
        # this, but it keeps track of it for the C fitting library.
        self.c_spline = cubicSplineC.CSpline2D(self.spline)

    def getPSF(self, z_value, shape = None, up_sample = 1, normalize = True):
        """
        This has the same arguments as the 3D version for convenience. 
        The z_value is ignored as long it is 0.0.
        """
        if (z_value != 0.0):
            print("Warning!! SplineToPSF2D got a non-zero z_value", z_value)

        assert(isinstance(up_sample, int)), "up_sample must be an integer."
        
        psf_size = up_sample * (self.spline_size - 1)
                
        psf = numpy.zeros((psf_size, psf_size))
        if((psf_size%2) == 0):
            for x in range(psf_size):
                for y in range(psf_size):
                    psf[y,x] = self.spline.f(float(y)/float(up_sample) + 0.5,
                                             float(x)/float(up_sample) + 0.5)
        else:
            for x in range(psf_size):
                for y in range(psf_size):
                    psf[y,x] = self.spline.f(float(y)/float(up_sample) + 1.0,
                                             float(x)/float(up_sample) + 1.0)

        if shape is not None:
            psf = self.upsize(psf, shape, up_sample)

        # Normalize if requested.
        if normalize:
            psf = psf/numpy.sum(psf)
            
        return psf

    def getScaledZ(self, z_value):
        return z_value * 0.0

    def getType(self):
        return "2D"
    

class SplineToPSF3D(SplineToPSF):

    def __init__(self, spline_file = None, **kwds):
        super(SplineToPSF3D, self).__init__(**kwds)

        spline_data = self.loadSplineFile(spline_file)

        # Check that this is not an old spline, which will be transposed.
        assert("version" in spline_data), "v0 spline file detected! Please re-measure!"
        assert(spline_data["version"] >= 2.0), "v0/v1 spline file detected! Please re-measure!"
        
        self.zmin = spline_data["zmin"]
        self.zmax = spline_data["zmax"]

        # The Python representation of the spline.        
        self.spline = spline3D.Spline3D(spline_data["spline"], spline_data["coeff"])
        self.spline_size = self.spline.getSize()

        # The C representation of the spline. This class does not use
        # this, but it keeps track of it for the C fitting library.
        self.c_spline = cubicSplineC.CSpline3D(self.spline)
        
    def getPSF(self, z_value, shape = None, up_sample = 1, normalize = True):
        """
        z_value needs to be inside the z range covered by the spline.
        z_value should be in nanometers.
        """
        assert(isinstance(up_sample, int)), "up_sample must be an integer."

        # Calculate PSF at requested z value.
        scaled_z = self.getScaledZ(z_value)
        
        psf_size = up_sample * (self.spline_size - 1)

        psf = numpy.zeros((psf_size, psf_size))

        #
        # Only extensively tested for psf_size even and up_sample = 1 as
        # this is what Spliner() uses.
        #
        if((psf_size%2) == 0):
            for x in range(psf_size):
                for y in range(psf_size):
                    psf[y,x] = self.spline.f(scaled_z,
                                             float(y)/float(up_sample) + 0.5,
                                             float(x)/float(up_sample) + 0.5)
        else:
            for x in range(psf_size):
                for y in range(psf_size):
                    psf[y,x] = self.spline.f(scaled_z,
                                             float(y)/float(up_sample) + 1.0,
                                             float(x)/float(up_sample) + 1.0)

        if shape is not None:
            psf = self.upsize(psf, shape, up_sample)

        # Normalize if requested.
        if normalize:
            psf = psf/numpy.sum(psf)
            
        return psf

    def getScaledZ(self, z_value):
        return float(self.spline_size) * (z_value - self.zmin) / (self.zmax - self.zmin)

    def getType(self):
        return "3D"

    def rescaleZ(self, z_value):
        spline_range = self.zmax - self.zmin
        return 1.0e-3*(z_value * spline_range / self.spline_size + self.zmin)


def loadSpline(spline_file):

    with open(spline_file, 'rb') as fp:
        spline_data = pickle.load(fp)
    if (spline_data["type"] == "3D"):
        return SplineToPSF3D(spline_data)
    else:
        return SplineToPSF2D(spline_data)
