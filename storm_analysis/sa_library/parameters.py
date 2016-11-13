#!/usr/bin/python
#
# Handles parsing analysis xml files.
#
# Hazen 10/13
#

import numpy
import os

from xml.etree import ElementTree


class ParametersException(Exception):
    pass


class Parameters(object):
    """
    All the parameters that are understood by various storm-analysis programs.
    """
    def __init__(self):

        #
        # All valid parameters (in alphabetical order).
        #
        self.attr = {

            # This is what the camera reads with the shutter closed.
            "baseline" : [float, None],

            # This file contains the sCMOS calibration data for the region of the camera
            # that the movie comes from. It consists of 3 numpy arrays, [offset, variance, gain],
            # each of which is the same size as a frame of the movie that is to be analyzed.
            # This can be generated for a camera using camera_calibration.py and (if it needs
            # to be resliced), reslice_calibration.py.
            "camera_calibration" : [str, None],

            # Z fit cutoff (used when z is calculated later from wx, wy).
            "cutoff" : [float, None],

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
            
            # Tracking parameter, frame descriptor string
            # 0 - activation frame
            # 1 - non-specific frame
            # 2 - channel1 frame
            # 3 - channel2 frame
            # 4 - etc..
            "descriptor" : [str, None],

            # Do z fitting (or not), only relevant for "3d" fitting.
            "do_zfit" : [int, None],

            # Do drift correction 0 = No.
            "drift_correction" : [int, None],

            # Number of frames in each (drift correction) sub-STORM image.
            "frame_step" : [int, None],
            
            # Maximum number of iterations for new peak finding.
            "iterations" : [int, None],

            # The frame to stop analysis on, -1 = analyze to the end of the film.
            "max_frame" : [int, None],

            # Model is one of 2dfixed, 2d, 3d, or Z.
            #
            #  2dfixed - fixed sigma 2d gaussian fitting.
            #  2d - variable sigma 2d gaussian fitting.
            #  3d - x, y sigma are independently variable, z will be fit after peak fitting.
            #  Z - x, y sigma depend on z, z is fit as part of peak fitting.
            #
            # This is used by 3D-DAOSTORM and sCMOS.
            "model" : [str, None],

            # CCD orientation, generally you should use "normal", but if you want to compare
            # the analysis with older versions of Insight3 you'll sometimes find that
            # "inverted" works best.
            "orientation" : [str, None],
            
            # CCD pixel size (in nm).
            "pixel_size" : [float, None],

            # Radius for matching peaks from frame to frame. Localizations that are closer than
            # this value (in pixels) in adjacent frames (ignoring activation frames) are assumed
            # to come from the same emitter and are averaged together to create a (hopefully) 
            # more accurately localized emitter. If this is zero then no matching will be done.
            "radius" : [float, None],
            
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

            # The frame to start analysis on, -1 = start at the beginning of the film.
            "start_frame" : [int, None],

            # This is basically the same as the minimum height parameter for peak finding in
            # Insight3. You should use a number roughly equal to the value of the brightest
            # pixel (minus the CCD baseline) in the dimmest peak that you want to detect. If
            # this is too low more background will be detected. If it is too high more peaks
            # will be missed.
            "threshold" : [float, None],

            # wx vs z parameters.
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
            
            # X start of the analysis AOI, leave as None to start at the edge of the image.
            "x_start" : [int, None],

            # X end of the analysis AOI, leave as None to end at the edge of the image.
            "x_stop" : [int, None],
            
            # Y start of the analysis AOI, leave as None to start at the edge of the image.
            "y_start" : [int, None],

            # Y end of the analysis AOI, leave as None to end at the edge of the image.
            "y_stop" : [int, None],
            
            }

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

    def getZRange(self):
        """
        Get z range.
        """
        return [self.getAttr("min_z", -0.5), self.getAttr("max_z", 0.5)]

    def hasAttr(self, name):
        return name in self.attr

    def setAttr(self, name, value):
        if self.hasAttr(self, name):
            if (type(value) == self.attr[name][0]):
                self.attr[name] = value
            else:
                ParametersException(value, "is the wrong type for", name, "expected", self.attr[name][0], "got", type(value))
        else:
            ParametersException(name, "is not a valid parameter name.")
        

def fromFile(parameters_file):
    """
    Create a Parameters object from an XML file.
    """
    params = Parameters()

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
            fname = os.path.join(dirname, value)
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
