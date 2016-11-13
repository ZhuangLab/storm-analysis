#!/usr/bin/python
#
# Handles parsing analysis xml files.
#
# Hazen 10/13
#

import numpy
import os

from xml.etree import ElementTree


class filename(str):
    pass


class ParametersException(Exception):
    pass


class Parameters(object):
    """
    Base parameters class.
    """
    def __init__(self):

        # This dictionary will contain all of the valid analysis parameters.
        # The format is "key" : [type, value].
        self.attr = {}

    def getAttr(self, name, default = None):
        try:
            if self.attr[name][1] is None:
                if default is None:
                    ParametersException(name, "is not initialized and no default was specified.")
                else:
                    return default
            else:
                return self.attr[name][1]
        except KeyError:
            ParametersException(name, "is not a valid parameter name.")
        
    def hasAttr(self, name):
        return name in self.attr

    def setAttr(self, name, value):
        if self.hasAttr(self, name):
            if isinstance(value, self.attr[name][0]):
                self.attr[name] = value
            else:
                ParametersException(value, "is the wrong type for", name, "expected", self.attr[name][0], "got", type(value))
        else:
            ParametersException(name, "is not a valid parameter name.")


class ParametersAnalysis(Parameters):
    """
    The parameters that are common to all of the STORM analysis programs.
    """
    def __init__(self):
        Parameters.__init__(self)
        
        self.attr.update({

            ##
            # Analysis parameters.
            ##

            # This is what the camera reads with the shutter closed.
            "baseline" : [float, None],

            # The frame to stop analysis on, -1 = analyze to the end of the film.
            "max_frame" : [int, None],

            # CCD orientation, generally you should use "normal", but if you want to compare
            # the analysis with older versions of Insight3 you'll sometimes find that
            # "inverted" works best.
            "orientation" : [str, None],
            
            # CCD pixel size (in nm).
            "pixel_size" : [float, None],
                        
            # The frame to start analysis on, -1 = start at the beginning of the film.
            "start_frame" : [int, None],

            # X start of the analysis AOI, leave as None to start at the edge of the image.
            "x_start" : [int, None],

            # X end of the analysis AOI, leave as None to end at the edge of the image.
            "x_stop" : [int, None],
            
            # Y start of the analysis AOI, leave as None to start at the edge of the image.
            "y_start" : [int, None],

            # Y end of the analysis AOI, leave as None to end at the edge of the image.
            "y_stop" : [int, None],
            

            ##
            # Tracking parameters
            ##

            # Tracking parameter, frame descriptor string
            # 0 - activation frame
            # 1 - non-specific frame
            # 2 - channel1 frame
            # 3 - channel2 frame
            # 4 - etc..
            "descriptor" : [str, None],

            # Radius for matching peaks from frame to frame. Localizations that are closer than
            # this value (in pixels) in adjacent frames (ignoring activation frames) are assumed
            # to come from the same emitter and are averaged together to create a (hopefully) 
            # more accurately localized emitter. If this is zero then no matching will be done.
            "radius" : [float, None],

            
            ##
            # Drift correction parameters.
            ##
            
            # This is the "scale" at which to render the sub-STORM images for drift correction.
            # Drift correction works by creating STORM images from frame_step sized groups 
            # of frames. These are rendered scaled by the d_scale parameter. For example, if
            # your data is 256x256 pixels then the drift-correction will create 512x512 sub-STORM 
            # images (for d_scale = 2) and then attempt to correlate these images to each other
            # to calculate the drift. Using a larger d_scale value creates higher resolution 
            # sub-STORM images, but they are also sparser so you might not see any improvement
            # in the drift correction.
            #
            # ... 2 is usually a good choice.
            "d_scale" : [int, None],

            # Do drift correction 0 = No.
            "drift_correction" : [int, None],

            # Number of frames in each (drift correction) sub-STORM image.
            "frame_step" : [int, None],

            })

    def getZRange(self):
        """
        Get z range.
        """
        return [self.getAttr("min_z", -0.5), self.getAttr("max_z", 0.5)]


class ParametersDAO(ParametersAnalysis):
    """
    Parameters that are specific to 3D-DAOSTORM analysis.
    """
    def __init__(self):
        ParametersAnalysis.__init__(self)

        self.attr.update({

            # Z fit cutoff (used when z is calculated later from wx, wy).
            "cutoff" : [float, None],
            
            # Do z fitting (or not), only relevant for "3d" fitting (see "model" parameter).
            "do_zfit" : [int, None],
            
            # Gaussian filter sigma, this is the sigma of a 2D gaussian to convolve the data with
            # prior to peak indentification. When your data has a low SNR this can help for peak
            # finding. For optimal sensitivity it should be the same as the expected sigma for your
            # peaks.
            #
            # You will need to adjust your threshold parameter as the threshold is now used for
            # peak finding in the convolved image and not the original image.
            #       
            # If you set it to zero (or comment it out) then this will not be performed, which can
            # make the analysis faster.
            #
            # Note: This is not relevant for sCMOS analysis.
            #
            "filter_sigma" : [float, None]

            # To be a peak it must be the maximum value within this radius (in pixels).
            "find_max_radius" : [int, None],
            
            # Maximum number of iterations for new peak finding.
            "iterations" : [int, None],

            # Model is one of 2dfixed, 2d, 3d, or Z.
            #
            #  2dfixed - fixed sigma 2d gaussian fitting.
            #  2d - variable sigma 2d gaussian fitting.
            #  3d - x, y sigma are independently variable, z will be fit after peak fitting.
            #  Z - x, y sigma depend on z, z is fit as part of peak fitting.
            "model" : [str, None],

            # Initial guess for sigma, this is in units of pixels. If you are using the 2dfixed
            # model then it needs to be pretty close to the correct value. For 2d it should be
            # close, probably within 50% or so of the average peak sigma or the fitting might fail
            # to converge on many peaks. 3d is similar to 2d. It should not effect fitting for Z
            # the model.
            "sigma" : [float, None],

            # Threshold for a maximum to considered a peak.
            #
            # Usually this is the same as the minimum height parameter for peak finding in
            # Insight3 (but see note above if you are also using "filter sigma"). You should
            # use a number roughly equal to the value of the brightest pixel (minus the CCD
            # baseline) in the dimmest peak that you want to detect. If this is too low more
            # background will be detected. If it is too high more peaks will be missed.
            #
            "threshold" : [float, None],
            
            # wx vs z parameters
            #
            # See Huang, Science 2008 for a more detailed explanation.
            #
            "wx_wo" : [float, None],
            "wx_c" : [float, None],
            "wx_d" : [float, None],
            "wxA" : [float, None],
            "wxB" : [float, None],
            "wxC" : [float, None],
            "wxD" : [float, None],

            # wy vs z parameters.
            "wy_wo" : [float, None],
            "wy_c" : [float, None],
            "wy_d" : [float, None],
            "wyA" : [float, None],
            "wyB" : [float, None],
            "wyC" : [float, None],
            "wyD" : [float, None],            
            })
        
    def getWidthParams(self, which, for_mu_Zfit = False):
        """
        Get "x" or "y" peak width versus z paremeters.
        """
        par = ["_wo", "_c", "_d", "A", "B", "C", "D"]
        np_par = numpy.zeros(len(par))
        for i,p in enumerate(par):
            attr = "w" + which + p
            if hasattr(parameters, attr):
                np_par[i] = getattr(parameters, attr)
        if for_mu_Zfit:
            np_par[0] = np_par[0]/parameters.pixel_size
            np_par[1] = np_par[1]*0.001
            np_par[2] = np_par[2]*0.001
        return np_par


class ParametersL1H(ParametersAnalysis):
    """
    Parameters that are specific to L1H analysis.
    """
    def __init__(self):
        ParametersAnalysis.__init__(self)
        
        self.attr.update({

            # A matrix file.
            "a_matrix" : [filename, None]

            # Epsilon, in Bo's paper he suggested 1.5 for poisson simulated data,
            # 2.1 for EMCCD data.
            "epsilon" : [float, None]

            })


class ParametersSCMOS(ParametersDAO):
    """
    Parameters that are specific to sCMOS analysis.
    """
    def __init__(self):
        ParametersDAO.__init__(self)

        self.attr.update({
            
            # This file contains the sCMOS calibration data for the region of the camera
            # that the movie comes from. It consists of 3 numpy arrays, [offset, variance, gain],
            # each of which is the same size as a frame of the movie that is to be analyzed.
            # This can be generated for a camera using camera_calibration.py and (if it needs
            # to be resliced), reslice_calibration.py.
            "camera_calibration" : [filename, None],

            # Initial guess for sigma, this is in units of pixels.
            #
            # For sCMOS analysis it is used as the sigma psf in image segmentation
            # (see Section 3.1 of the supplementary material of:
            #
            # "Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms"
            # Huang et al, Nature Methods, 2013.
            #
            # It also used to initialize fitting. If you are using the 2dfixed model then it
            # needs to be pretty close to the correct value. For 2d it should be close, probably 
            # within 50% or so of the average peak sigma or the fitting might fail to converge
            # on many peaks. 3d is similar to 2d. It should not effect fitting for Z the model.
            "sigma" : [float, None],

            # Threshold for a maximum to considered a peak. This has the same meaning as for
            # 3D-DAOSTORM, except that it is applied to the convolved image, so you will likely
            # need to use a smaller value.
            "threshold" : [float, None],
            
            })
    

class ParametersSpliner(ParametersAnalysis):
    """
    Parameters that are specific to Spliner analysis.
    """
    def __init__(self):
        ParametersAnalysis.__init__(self)

        self.attr.update({
        
            # This is the spline file to use for fitting. Based on the spline the analysis will
            # decide whether to do 2D or 3D spline fitting, 2D if the spline is 2D, 3D if the
            # spline is 3D.
            "spline" : [filename, None],

            # Use FISTA deconvolution for peak finding. If this is not set then the analysis
            # will be done using a matched filter for peak finding. This is much faster, but
            # possibly less accurate at higher densities.
            "use_fista" : [int, None],
            })


class ParametersSplinerSTD(ParametersSpliner):
    """
    Parameters that are specific to Spliner standard analysis.
    """
    def __init__(self):
        ParametersSpliner.__init__(self)

        self.attr.update({

            # To be a peak it must be the maximum value within this radius (in pixels).
            "find_max_radius" : [int, None],
            
            # Maximum number of iterations for new peak finding.
            "iterations" : [int, None],
            
            # Threshold for a maximum to considered a peak. This has the same meaning as for
            # 3D-DAOSTORM, except that it is applied to the image convolved with the PSF, so
            # you will likely need to use a smaller value.
            "threshold" : [float, None],
            
            # Z value(s) in nanometers at which we will perform convolution with the PSF for
            # the purposes of peak finding. If this is not specified the default value is z = 0.0. -->
            "z_value" : [list, None],
            
            })
        

class ParametersSplinerFISTA(ParametersSpliner):
    """
    Parameters that are specific to Spliner FISTA analysis.
    """
    def __init__(self):
        ParametersSpliner.__init__(self)

        self.attr.update({

            # Iterations of FISTA deconvolution to perform. The larger this value is the sharper
            # the peaks will be.
            "fista_iterations" : [int, None],

            # FISTA lambda value. Larger values will increase the sparsity of the deconvolved image.
            "fista_lambda" : [float, None],

            # The number of z-planes to use in the deconvolution, more planes will give higher
            # accuracy at the expense of running time.
            "fista_number_z" : [int, None],

            # Local maxima in the FISTA deconvolved image with values larger than this will input
            # into the fitter as localizations to be fit. This number should be roughly the minimum
            # peak height that would be considered real times the integral of a peak of this height.
            "fista_threshold" :  [float, None],

            # FISTA timestep. Larger values will cause FISTA to converge faster, but if the value is
            # too large FISTA will rapidly diverge.
            "fista_timestep" : [float, 0.1],

            # Radius of the rolling ball in pixels.
            "rb_radius" : [float, None],

            # Sigma in pixels of the gaussian smoothing to apply to the background estimate after
            # the rolling ball step.
            "rb_sigma" : [float, None]

            })
        
                         
def fromFile(params, parameters_file):
    """
    Set the attributes of the params Parameters object from an XML file.
    """
    settings = ElementTree.parse(parameters_file).getroot()
    for node in settings:
        slot = node.tag
        value = node.text
        ntype = node.attrib.get("type")
            
        if (ntype == "int"):
            params.setAttr(slot, int(value))
                
        elif (ntype == "int-array"):
            text_array = value.split(",")
            int_array = []
            for elt in text_array:
                int_array.append(int(elt))
                params.setAttr(slot, int_array)
                    
        elif (ntype == "float"):
            params.setAttr(slot, float(value))
                
        elif (ntype == "float-array"):
            text_array = value.split(",")
            float_array = []
            for elt in text_array:
                float_array.append(float(elt))
            params.setAttr(slot, float_array)
                
        elif (ntype == "string-array"):
            params.setAttr(slot, value.split(","))
                
        elif (ntype == "filename"):
            dirname = os.path.dirname(os.path.abspath(parameters_file))
            fname = filename(os.path.join(dirname, value))
            params.setAttr(slot, fname)
            
        else: # everything else is assumed to be a string
            params.setAttr(slot, value)

    params.setAttr("parameters_file", parameters_file)

    return params


#
# The MIT License
#
# Copyright (c) 2013 Zhuang Lab, Harvard University
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
