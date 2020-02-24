#!/usr/bin/env python
"""
Handles parsing analysis xml files & interaction with the object
created from same.

FIXME: The class heirarchy is too deep? I was trying to avoid redundant 
       parameters, but perhaps at the expense of clarity.

Hazen 10/13
"""
import copy
import numpy
import os

from xml.dom import minidom
from xml.etree import ElementTree

import storm_analysis


class ParametersException(storm_analysis.SAException):
    pass


class Parameters(object):
    """
    Base parameters class.
    """
    def __init__(self, **kwds):
        super(Parameters, self).__init__(**kwds)

        # This dictionary will contain all of the valid analysis parameters.
        # The format is "key" : [type, value, doc].
        self.attr = {}

        self.filename = None

    def changeAttr(self, name, value, node_type = None):
        """
        Use this for user interaction. Use setAttr() for machine interaction.
        """
        if not (name in self.attr):
            raise ParametersException("No such parameter", name)

        # Check that the user provided one of the available types
        # for a parameter with multiple possible types.
        #
        if isinstance(self.attr[name][0], tuple):
            if node_type is None:
                raise ParametersException("You must provide node_type for this attribute " + str(self.attr[name][0]))
            
            if not node_type in self.attr[name][0]:
                raise ParametersException("node_type must be one of " + str(self.attr[name][0]))

            self.attr[name][0] = node_type
        elif node_type is not None:
            if node_type != self.attr[name][0]:
                print("Warning! Ignoring incorrect node type", node_type, "for", name)
        
        self.attr[name][1] = value
        
    def getAttr(self, name, default = None):
        if not (name in self.attr):
            raise ParametersException("No such parameter", name)
        if self.attr[name][1] is not None:
            return self.attr[name][1]
        if default is None:
            raise ParametersException(name, "is not initialized and no default was specified.")
        else:
            return default
        
    def hasAttr(self, name):
        if not (name in self.attr):
            return False
        return (self.attr[name][1] is not None)

    def helpAttr(self, name):
        if not (name in self.attr):
            raise ParametersException("No such parameter", name)
        data = self.attr[name]
        h_string = "'" + name + "' type is " + str(self.attr[name][0]) + ".\n\n"
        if (len(data) > 2):
            # Strip leading spaces off multi-line help strings.
            return h_string + "\n".join(map(lambda x: x.lstrip(' '), data[2].split("\n")))
        else:
            return h_string + "No help available for " + name + "."
        
    def initAttr(self, nodes, warnings = True):
        """
        Set the attributes of the params Parameters 
        object using an ElementTree object.
        """
        for node in nodes:
            slot = node.tag
            value = node.text
            ntype = node.attrib.get("type")
            
            if (ntype == "int"):
                self.setAttr(slot, ntype, int(value), warnings)
                
            elif (ntype == "int-array"):
                text_array = value.split(",")
                int_array = []
                for elt in text_array:
                    int_array.append(int(elt))
                self.setAttr(slot, ntype, int_array, warnings)

            elif (ntype == "float"):
                self.setAttr(slot, ntype, float(value), warnings)
                
            elif (ntype == "float-array"):
                text_array = value.split(",")
                float_array = []
                for elt in text_array:
                    float_array.append(float(elt))
                self.setAttr(slot, ntype, float_array, warnings)
                
            elif (ntype == "string-array"):
                self.setAttr(slot, ntype, value.split(","), warnings)
                
            elif (ntype == "filename"):
                dirname = os.path.dirname(os.path.abspath(self.filename))
                fname = os.path.join(dirname, value)
                self.setAttr(slot, ntype, fname, warnings)
            
            else: # everything else is assumed to be a string
                self.setAttr(slot, ntype, value, warnings)

        return self

    def initFromFile(self, parameters_file, warnings = True):
        """
        Set the attributes of the params Parameters object from an XML file.
        """
        self.filename = parameters_file
        
        # Save the parameters file name.
        self.attr.update({"parameters_file" : ["parameters_filename", None]})
        self.setAttr("parameters_file", "parameters_filename", self.filename, warnings)
    
        nodes = ElementTree.parse(parameters_file).getroot()
        return self.initAttr(nodes, warnings = warnings)
        
    def initFromString(self, parameters_string, warnings = True):
        """
        Set the attributes of the params Parameters object from an XML string.
        """
        nodes = ElementTree.fromstring(parameters_string)
        return self.initAttr(nodes, warnings = warnings)

    def prettyPrint(self):
        for fname in sorted(self.attr):
            print(fname + " (" + str(self.attr[fname][0]) + ") - " + str(self.attr[fname][1]))
        
    def setAttr(self, name, node_type, value, warnings = True):
        if name in self.attr:

            # Check that this was parsed as the correct type.
            good = False
            if isinstance(self.attr[name][0], tuple):
                if node_type in self.attr[name][0]:
                    good = True
                    self.attr[name][0] = node_type
            elif (node_type == self.attr[name][0]):
                good = True
                
            if good:
                self.attr[name][1] = value
            else:
                raise ParametersException(value, "is the wrong type for", name, "expected", self.attr[name][0], "got", node_type)
        elif warnings:
            # Exception for v1.0 type parameters.
            if (name == "baseline"):
                raise ParametersException(name + " is a version 1.0 parameter, please update the XML file!")

            # Not sure whether an Exception or a warning is the best choice here..
            print("Warning!!", name, "is not a relevant parameter!!")

    def update(self, other):
        """
        Add parameters in other to self. Duplicates are not allowed.
        """
        for key in other.attr:
            if key in self.attr:
                if not (key == "parameters_file"):
                    raise ParametersException("'" + key + "' exists in both parameter objects!")
            else:
                self.attr[key] = copy.deepcopy(other.attr[key])

    def toXMLElementTree(self, remove_paths):
        """
        Convert back to ElementTree object.
        """
        etree = ElementTree.Element("settings")
        for fname in sorted(self.attr):
            if self.attr[fname][1] is not None:
                field = ElementTree.SubElement(etree, fname)
                if (self.attr[fname][0] == "filename"):
                    if remove_paths:
                        field.text = os.path.basename(self.attr[fname][1])
                    else:
                        field.text = self.attr[fname][1]
                elif "array" in self.attr[fname][0]:
                    field.text = ",".join(map(str, self.attr[fname][1]))
                else:
                    field.text = str(self.attr[fname][1])
                field.set("type", str(self.attr[fname][0]))

        return etree

    def toXMLFile(self, filename, pretty = False, remove_paths = True):
        """
        Write to a XML file. This file will not be nicely formatted..
        """
        
        # Is this a string?
        if isinstance(filename, str):
            with open(filename, "wb") as fp:
                fp.write(self.toXMLString(pretty = pretty, remove_paths = remove_paths))

        # If not, assume it is a file pointer.
        else:
            filename.write(self.toXMLString(pretty = pretty, remove_paths = remove_paths))
    
    def toXMLString(self, pretty = False, remove_paths = True):
        """
        Convert back to an XML string.
        """
        unpretty_string = ElementTree.tostring(self.toXMLElementTree(remove_paths), 'ISO-8859-1')
        if pretty:
            reparsed = minidom.parseString(unpretty_string)
            return reparsed.toprettyxml(indent="  ", encoding = "ISO-8859-1")
        return unpretty_string


class ParametersCommon(Parameters):
    """
    The parameters that are common to all of the STORM analysis programs.
    """
    def __init__(self, **kwds):
        super(ParametersCommon, self).__init__(**kwds)
        
        self.attr.update({
            
            ##
            # Analysis parameters.
            ##
            "aoi_radius" : ["int", None,
                            """Radius in pixels for a circular analysis AOI.
                            'x_center' and 'y_center' are also required."""],

            "max_frame" : ["int", None,
                           "The frame to stop analysis on, -1 = analyze to the end of the film."],
            
            "max_z" : ["float", None,
                       "Maximum z value for z fitting, specified in um."],
            
            "min_z" : ["float", None,
                       "Minimum z value for z fitting, specified in um."],
            
            "pixel_size" : ["float", None,
                            "CCD pixel size (in nm)."],
            
            "start_frame" : ["int", None,
                             "The frame to start analysis on, -1 = start at the beginning of the film."],

            "static_background_estimate" : ["int", None,
                                            """If this is set, and set to a number greater than 0, then the analysis will
                                            estimate the background by using the average over this number of frames.
                                            
                                            If this is not set, or set to 0, the background is estimated separately
                                            for each frame."""],

            "verbosity" : ["int", 1,
                           """Frequency in frames at which to provide feedback. The default is every 
                           frame (verbosity = 1). This should be an integer >= 1."""],

            "x_center" : ["int", None,
                          """X center for a circular analysis AOI.
                          'y_center' and 'aoi_radius' are also required."""],
            
            "x_start" : ["int", None,
                         "X start of the analysis AOI, leave as None to start at the edge of the image."],

            "x_stop" : ["int", None,
                        "X end of the analysis AOI, leave as None to end at the edge of the image."],
            
            "y_center" : ["int", None,
                          """Y center for a circular analysis AOI.
                          'x_center' and 'aoi_radius' are also required."""],
            
            "y_start" : ["int", None,
                         "Y start of the analysis AOI, leave as None to start at the edge of the image."],

            "y_stop" : ["int", None,
                        "Y end of the analysis AOI, leave as None to end at the edge of the image."],
            
            ##
            # Tracking parameters
            ##
            
            "descriptor" : ["string", None,
                            """Tracking parameter, frame descriptor string
                            0 - activation frame
                            1 - non-specific frame
                            2 - channel1 frame
                            3 - channel2 frame
                            4 - etc.."""],

            "max_gap" : ["int", None,
                         """The maximum gap in frames between two localizations that within radius of each
                         other in X/Y to be linked together in a track."""],

            "radius" : ["float", None,
                        """Radius for matching peaks from frame to frame. Localizations that are closer than
                        this value (in pixels) in adjacent frames (ignoring activation frames) are assumed
                        to come from the same emitter and are averaged together to create a (hopefully) 
                        more accurately localized emitter. If this is zero then no matching will be done."""],
            
            ##
            # Drift correction parameters.
            ##
            
            "d_scale" : ["int", None,
                         """This is the "scale" at which to render the sub-STORM images for drift correction.
                         Drift correction works by creating STORM images from frame_step sized groups 
                         of frames. These are rendered scaled by the d_scale parameter. For example, if
                         your data is 256x256 pixels then the drift-correction will create 512x512 sub-STORM 
                         images (for d_scale = 2) and then attempt to correlate these images to each other
                         to calculate the drift. Using a larger d_scale value creates higher resolution 
                         sub-STORM images, but they are also sparser so you might not see any improvement
                         in the drift correction.
                         
                         ... 2 is usually a good choice."""],

            "drift_correction" : ["int", None,
                                  "Do drift correction, 0 = No."],

            "frame_step" : ["int", None,
                            "Number of frames in each (drift correction) sub-STORM image."],

            "z_correction": ["int", None,
                             "Do z drift correction, 0 = No."],


            ##
            # File conversion.
            ##

            "convert_to" : ["string", None,
                            """Specify what, if any, formats to convert the output HDF5 file into upon completion
                            of the analysis.
                            
                            Options are .bin and .txt.
                            Use a comma separated list if you want both. i.e. ".bin, .txt"."""],
            })

    def getZRange(self):
        """
        Get z range.
        """
        return [self.getAttr("min_z", -0.5), self.getAttr("max_z", 0.5)]


class ParametersDecon(Parameters):
    """
    Parameters that are specific to deconvolution (for peak identification).
    """
    def __init__(self, **kwds):
        super(ParametersDecon, self).__init__(**kwds)

        self.attr.update({

            "background_estimator" : ["string", None,
                                      """Method to use for background estimation, either 'RollingBall' or 'Wavelet'."""],
            })



        
sigma_doc_string = """This is the estimated sigma of the PSF in pixels.

It serves several purposes:
\t(1) It is used in most of the analysis approaches as a measure of the
\t\tpeak to peak distance at which peak fits do not substantially
\t\teffect each other.

\t(2) In most of the analysis approaches, if two peaks are closer than
\t\tthis distance then the dimmer one will be discarded.

\t(3) In 3D-DAOSTORM, sCMOS and Multi-plane 3D-DAOSTORM analysis it is
\t\talso used as the initial guess for the peak sigma.

\t(4) In 3D-DAOSTORM, sCMOS and the peak widths are constrained to the
\t\trange 0.5x - 5x sigma.

\t(5) In Multi-plane 3D-DAOSTORM the peak widths are constrained to the
\t\trange 0.5x - 3x sigma.

Also, if you are using the '2dfixed' model then it needs to be pretty close 
to the correct value. For '2d' it should be close, probably within 50% or 
so of the average peak sigma or the fitting might fail to converge on many 
peaks. '3d' is similar to '2d'.
"""


class ParametersFitters(ParametersCommon):
    """
    The parameters that are common to the fitting based approaches, i.e. 3D-DAOSTORM, 
    sCMOS, Spliner and Multi-plane.
    """
    def __init__(self, **kwds):
        super(ParametersFitters, self).__init__(**kwds)
        
        self.attr.update({

            "background_sigma" : ["float", None,
                                  """Background filter sigma, this is the sigma of a 2D gaussian to convolve the
                                  data in order to estimate the background."""],

            "fftw_estimate" : ["int", 0,
                               """ FFTW should estimate the best FFT plan instead of measuring which is best. This 
                               can help speed the analysis of short movies that are very large in XY. 0 = (default) 
                               FFTW will measure the best FFT plan. 1 = FFTW will estimate the best FFT plan."""],

            "fit_error_model" : ["string", "MLE",
                                 """Specify which fitting error model to use.

                                 The default model is 'MLE' (Maximum Likelihood Estimation). Other options include  
                                 'ALS' (Anscombe Least Squares), 'LS' (Least Squares), 'DWLS' (Data Weighted Least 
                                 Squares) and 'FWLS' (Fit Weighted Least Squares).
                           
                                 Note that most fitters only support the 'MLE' fitting error model."""],

            "find_max_radius" : [("int", "float"), None,
                                 "To be a peak it must be the maximum value within this radius (in pixels)."],

            "iterations" : ["int", None,
                            "Maximum number of iterations for new peak finding."],

            "no_fitting" : ["int", None,
                            """If this is non-zero then we won't do any fitting iterations. This is useful for
                            testing the finder, as well as how accurately we're initializing the peak
                            parameter values."""],

            "peak_locations" : ["filename", None,
                                """This is for is you already know where your want fitting to happen, as
                                for example in a bead calibration movie and you just want to use the
                                approximate locations as inputs for fitting.
                                
                                peak_locations can be a text file with the peak x, y, height and background
                                values as white spaced columns (x and y positions are in pixels as
                                determined using visualizer).
                                
                                1.0 2.0 1000.0 100.0
                                10.0 5.0 2000.0 200.0

                                peak_locations can also be a HDF5 file. This is treated in a similar
                                fashion to the text file in that all of the locations are loaded.
                                If the fields 'xsigma' or 'ysigma' exist they will be used for the
                                initial X/Y sigma values of the localization.
                                ..."""],
            
            "sigma" : ["float", None, sigma_doc_string],

            "threshold" : ["float", None,
                           """Threshold for a maximum to considered a peak.
                           
                           This is the threshold for peak finding in units of signal to background. A
                           value of 3 for example corresponds to only selecting peaks with an (estimated)
                           signal to background ratio of 3.
                           
                           You probably want a value of at least 5."""],
            })

        
ww_doc_string = """wx/wy vs z parameters. Units are nanometers or dimensionless.

See Huang, Science 2008 for a more detailed explanation.
"""


class ParametersDAOsCMOS(ParametersFitters):
    """
    The parameters that are common to 3D-DAOSTORM and sCMOS.
    """
    def __init__(self, **kwds):
        super(ParametersDAOsCMOS, self).__init__(**kwds)
        
        self.attr.update({

            "cutoff" : ["float", None,
                        """Z fit cutoff (used when z is calculated later from wx, wy)."""],

            "do_zfit" : ["int", None,
                         """Do z fitting (or not), only relevant for "3d" fitting (see "model" parameter)."""],

            "foreground_sigma" : ["float", None,
                                  """Foreground filter sigma, this is the sigma of a 2D gaussian to convolve the data with
                                  prior to peak identification. When your data has a low SNR this can help for peak
                                  finding. For optimal sensitivity it should be the same as the expected sigma for your
                                  peaks.
                                  
                                  If you set it to zero (or comment it out) then this will not be performed, which can
                                  make the analysis faster."""],

            "model" : ["string", None,
                       """Model is one of 2dfixed, 2d, 3d, or Z.
                       
                       2dfixed - fixed sigma 2d gaussian fitting.
                       2d - variable sigma 2d gaussian fitting.
                       3d - x, y sigma are independently variable, z will be fit after peak fitting.
                       Z - x, y sigma depend on z, z is fit as part of peak fitting."""],

            "roi_size" : ["int", None,
                          """This is the size of the fitting ROI in pixels. If it is not specified than the value
                          will be calculated based on the sigma parameter value and the fitting model type"""],

            "sigma" : ["float", None, sigma_doc_string],

            "sigma_range" : ["float-array", None,
                             """A two element array that specifies the minimum and maximum sigma values to
                             allow when fitting for the peak width. If this is not specified the default
                             is [0.5 * sigma, 5.0 * sigma]."""],

            "wx_wo" : ["float", None, ww_doc_string],
            "wx_c" : ["float", None, ww_doc_string],
            "wx_d" : ["float", None, ww_doc_string],
            "wxA" : ["float", None, ww_doc_string],
            "wxB" : ["float", None, ww_doc_string],
            "wxC" : ["float", None, ww_doc_string],
            "wxD" : ["float", None, ww_doc_string],

            "wy_wo" : ["float", None, ww_doc_string],
            "wy_c" : ["float", None, ww_doc_string],
            "wy_d" : ["float", None, ww_doc_string],
            "wyA" : ["float", None, ww_doc_string],
            "wyB" : ["float", None, ww_doc_string],
            "wyC" : ["float", None, ww_doc_string],
            "wyD" : ["float", None, ww_doc_string],

            "z_value" : ["float", None,
                         "The starting z value for fitting. If this is not specified it defaults to 0.0."],

            "z_step" : ["float", None,
                        """The z step size for finding the optimal z value when using the 3d model. If
                        this is not specified it defaults to 1 nanometer. Units are microns."""],            
            })

    def getWidthParams(self, for_mu_Zfit = False):
        """
        Get "x" or "y" peak width versus z paremeters.
        """
        par = ["_wo", "_c", "_d", "A", "B", "C", "D"]
        wx_params = numpy.zeros(len(par))
        wy_params = numpy.zeros(len(par))
        for i, p in enumerate(par):
            wx_params[i] = self.getAttr("wx" + p, 0.0)
            wy_params[i] = self.getAttr("wy" + p, 0.0)

        for np_par in [wx_params, wy_params]:
            if for_mu_Zfit:
                np_par[0] = np_par[0]/self.getAttr("pixel_size")
                np_par[1] = np_par[1]*0.001
                np_par[2] = np_par[2]*0.001

        return [wx_params, wy_params]

    
##
## Parameter definitions for the different types of analysis.   
##

class ParametersADMM(ParametersDecon):
    """
    Parameters for ADMM deconvolution and localization identification.
    """
    def __init__(self, **kwds):
        super(ParametersADMM, self).__init__(**kwds)

        self.attr.update({
            
            "admm_iterations" : ["int", None,
                                  """Iterations of ADMM deconvolution to perform. The larger this value is the sharper
                                  the peaks will be."""],
            
            "admm_lambda" : ["float", None,
                              "ADMM lambda value. Larger values will increase the sparsity of the deconvolved image."],

            "admm_number_z" : ["int", None,
                                """The number of z-planes to use in the deconvolution. More planes will give higher
                                accuracy at the expense of running time, but see the note about z_value in
                                ParametersSplinerSTD as that also applies here."""],

            "admm_rho" : ["float", None,
                          """ADMM rho value. A value like 0.1 seems to work well."""],
            
            "admm_threshold" :  ["float", None,
                                  """Local maxima in the ADMM deconvolved image with values larger than this will input
                                  into the fitter as localizations to be fit. This number should be roughly the minimum
                                  peak height that would be considered real times the integral of a peak of this height."""],

            })

        
class ParametersDAO(ParametersDAOsCMOS):
    """
    Parameters that are specific to 3D-DAOSTORM analysis.
    """
    def __init__(self, **kwds):
        super(ParametersDAO, self).__init__(**kwds)

        self.attr.update({
     
            "camera_gain" : ["float", None,
                             """Conversion factor to go from camera ADU to photo-electrons. Units are ADU/e-, so the
                             camera ADU values will be divided by this number to convert to photo-electrons."""],
            
            "camera_offset" : ["float", None,
                               "This is what the camera reads with the shutter closed."],

            })



class Parameters3denseSTORM(ParametersDecon):
    """
    Parameters for 3denseSTORM deconvolution and localization identification.
    """
    def __init__(self, **kwds):
        super(Parameters3denseSTORM, self).__init__(**kwds)

        self.attr.update({

            "ds3_beta" : ["float", None,
                         "3denseSTORM beta value. The minimum SNR for detected localization. 2.44 recommended by Ovesny."],
            
            "ds3_eta" : ["float", None,
                         "3denseSTORM eta value. 9.0e-3 recommended by Ovesny"],
            
            "ds3_iterations" : ["int", None,
                                """Iterations of 3denseSTORM deconvolution to perform. The larger this value is the
                                sharper the peaks will be."""],

            "ds3_micro" : ["float", None,
                          "3denseSTORM micro value. 0.1 recommended by Ovesny"],

            "ds3_number_z" : ["int", None,
                              """The number of z-planes to use in the deconvolution. More planes will give higher
                              accuracy at the expense of running time, but see the note about z_value in
                              ParametersSplinerSTD as that also applies here."""],
            
            "ds3_threshold" :  ["float", None,
                                """Local maxima in the 3denseSTORM deconvolved image with values larger than this will 
                                input into the fitter as localizations to be fit. Depending on the image this is
                                probably something like 500."""],
            })


camera_gain_doc_string = """Conversion factor to go from camera ADU to photo-electrons. Units are ADU/e-, so the
camera ADU values will be divided by this number to convert to photo-electrons (e-)."""

camera_offset_doc_string = """This is what the camera reads with the shutter closed."""


class ParametersFISTA(ParametersDecon):
    """
    Parameters for FISTA deconvolution and localization identification.
    """
    def __init__(self, **kwds):
        super(ParametersFISTA, self).__init__(**kwds)

        self.attr.update({
            
            "fista_iterations" : ["int", None,
                                  """Iterations of FISTA deconvolution to perform. The larger this value is the sharper
                                  the peaks will be."""],

            "fista_lambda" : ["float", None,
                              "FISTA lambda value. Larger values will increase the sparsity of the deconvolved image."],

            "fista_number_z" : ["int", None,
                                """The number of z-planes to use in the deconvolution. More planes will give higher
                                accuracy at the expense of running time, but see the note about z_value in
                                ParametersSplinerSTD as that also applies here."""],

            "fista_threshold" :  ["float", None,
                                  """Local maxima in the FISTA deconvolved image with values larger than this will input
                                  into the fitter as localizations to be fit. This number should be roughly the minimum
                                  peak height that would be considered real times the integral of a peak of this height."""],

            "fista_timestep" : ["float", None,
                                """FISTA timestep. Larger values will cause FISTA to converge faster, but if the value is
                                too large FISTA will rapidly diverge."""],
            })
        

class ParametersL1H(ParametersCommon):
    """
    Parameters that are specific to L1H analysis.
    """
    def __init__(self, **kwds):
        super(ParametersL1H, self).__init__(**kwds)
        
        self.attr.update({

            "a_matrix" : ["filename", None,
                          "A matrix file"],

            "camera_gain" : ["float", None, camera_gain_doc_string],
            
            "camera_offset" : ["float", None, camera_offset_doc_string],

            "epsilon" : ["float", None,
                         """Epsilon, in Bo's paper he suggested 1.5 for poisson simulated data,
                         2.1 for EMCCD data."""],
            })

        
pm_cal_doc_string = "These are the (sCMOS) camera calibration files for each camera."
pm_ext_doc_string = """These are the extension onto the base name for each of the movies. If
your multi-plane data is all in a single file you will need to split it
into separate files first. Also, channel0 is "special" in that the
analysis will use this channel to figure out the movie length, and this
is the movie that will be used to calculate the movie signature / hash ID."""
pm_offset_doc_string = """These are to deal with the problem of not being able to get all the
cameras to start at the same time in a multi-camera setup. They specify
the (relative) offset for each channel in frames."""


class ParametersMultiplane(ParametersFitters):
    """
    Parameters that are specific to multi-plane analysis. Currently this is
    limited to a maximum of 8 planes.
    """
    def __init__(self, **kwds):
        super(ParametersMultiplane, self).__init__(**kwds)

        self.attr.update({


            "channel0_cal" : ["filename", None, pm_cal_doc_string],
            "channel1_cal" : ["filename", None, pm_cal_doc_string],
            "channel2_cal" : ["filename", None, pm_cal_doc_string],
            "channel3_cal" : ["filename", None, pm_cal_doc_string],
            "channel4_cal" : ["filename", None, pm_cal_doc_string],
            "channel5_cal" : ["filename", None, pm_cal_doc_string],
            "channel6_cal" : ["filename", None, pm_cal_doc_string],
            "channel7_cal" : ["filename", None, pm_cal_doc_string],
            
            "channel0_ext" : ["string", None, pm_ext_doc_string],
            "channel1_ext" : ["string", None, pm_ext_doc_string],
            "channel2_ext" : ["string", None, pm_ext_doc_string],
            "channel3_ext" : ["string", None, pm_ext_doc_string],
            "channel4_ext" : ["string", None, pm_ext_doc_string],
            "channel5_ext" : ["string", None, pm_ext_doc_string],
            "channel6_ext" : ["string", None, pm_ext_doc_string],
            "channel7_ext" : ["string", None, pm_ext_doc_string],
            
            "channel0_offset" : ["int", None, pm_offset_doc_string],
            "channel1_offset" : ["int", None, pm_offset_doc_string],
            "channel2_offset" : ["int", None, pm_offset_doc_string],
            "channel3_offset" : ["int", None, pm_offset_doc_string],
            "channel4_offset" : ["int", None, pm_offset_doc_string],
            "channel5_offset" : ["int", None, pm_offset_doc_string],
            "channel6_offset" : ["int", None, pm_offset_doc_string],
            "channel7_offset" : ["int", None, pm_offset_doc_string],


            "find_max_radius" : [("int", "float"), None,
                                 "To be a peak it must be the maximum value within this radius (in pixels)."],

            "iterations" : ["int", None,
                            "Maximum number of iterations for new peak finding."],
            
            "mapping" :  ["filename", None,
                          """This file contains the mapping between each of the planes. Typically it
                          is created using multi_plane/mapper.py."""],
            
            "z_value" : ["float-array", None,
                         """Z value(s) in microns at which we will perform convolution with the PSF for
                         the purposes of peak finding. If this is not specified the default value is
                         z = [0.0]. These are also the starting z values for fitting.
                         
                         If you are using this analysis to analyze single plane data then see the note
                         in the ParametersSpliner section regarding PSF Z degeneracy."""],
            })

        
pma_psf_doc_string = """These are the PSF files to use for PSF FFT fitting. There should be one of
them for each plane. The PSFs should have the same numbering as the mappings,
i.e. 'psf0' should be the PSF for channel0, etc."""
pma_pfn_doc_string = """These are the pupil function files to use for PupilFn fitting. There should
be one of them for each plane. The pupil functions should have the same
numbering as the mappings, i.e. 'pupil_fn0' should be the pupil function for
channel0, etc."""
pma_spl_doc_string = """These are the spline files to use for fitting. There should be one of them
for each plane. The splines should have the same numbering as the mappings,
i.e. 'spline0' should be the spline for channel0, etc."""


class ParametersMultiplaneArb(ParametersMultiplane):
    """
    Parameters that are specific to multi-plane analysis with the spline,
    pupil function or PSF FFT model.

    Note: multi-plane analysis works with either splines, pupil functions
          or PSF FFT models for the PSF, but you have to choose one of them,
          they cannot be mixed and matched.
    """
    def __init__(self, **kwds):
        super(ParametersMultiplaneArb, self).__init__(**kwds)

        self.attr.update({
            
            "independent_heights" : ["int", None,
                                     """Channel heights are independent, 0 = No. For multi-plane fitting you want
                                     this to be 0, for multi-color fitting you want this to be 1."""],

            "psf0" :  ["filename", None, pma_psf_doc_string],
            "psf1" :  ["filename", None, pma_psf_doc_string],
            "psf2" :  ["filename", None, pma_psf_doc_string],
            "psf3" :  ["filename", None, pma_psf_doc_string],
            "psf4" :  ["filename", None, pma_psf_doc_string],
            "psf5" :  ["filename", None, pma_psf_doc_string],
            "psf6" :  ["filename", None, pma_psf_doc_string],
            "psf7" :  ["filename", None, pma_psf_doc_string],            

            "pupilfn0" :  ["filename", None, pma_pfn_doc_string],
            "pupilfn1" :  ["filename", None, pma_pfn_doc_string],
            "pupilfn2" :  ["filename", None, pma_pfn_doc_string],
            "pupilfn3" :  ["filename", None, pma_pfn_doc_string],
            "pupilfn4" :  ["filename", None, pma_pfn_doc_string],
            "pupilfn5" :  ["filename", None, pma_pfn_doc_string],
            "pupilfn6" :  ["filename", None, pma_pfn_doc_string],
            "pupilfn7" :  ["filename", None, pma_pfn_doc_string],
            
            "spline0" :  ["filename", None, pma_spl_doc_string],
            "spline1" :  ["filename", None, pma_spl_doc_string],
            "spline2" :  ["filename", None, pma_spl_doc_string],
            "spline3" :  ["filename", None, pma_spl_doc_string],
            "spline4" :  ["filename", None, pma_spl_doc_string],
            "spline5" :  ["filename", None, pma_spl_doc_string],
            "spline6" :  ["filename", None, pma_spl_doc_string],
            "spline7" :  ["filename", None, pma_spl_doc_string],
            
            "weights" : ["filename", None,
                         """This specifies the file that contains how to optimally weight the updates
                         for each parameter from each plane as a function of z. If this is not
                         specified all planes will get equal weight."""],
            })

        
pmd_ww_doc_string = """Width versus z parameters for each plane. Units are nanometers or dimensionless.

The format is [w_wo, w_c, w_d, wA, wB, wC, wD]. These have the same meaning as
the wx vs z parameters in ParametersDAOsCMOS."""


class ParametersMultiplaneDao(ParametersMultiplane):
    """
    Parameters that are specific to multi-plane analysis using 3D-DAOSTORM
    (Gaussian) fitting.
    """
    def __init__(self, **kwds):
        super(ParametersMultiplaneDao, self).__init__(**kwds)

        self.attr.update({
            "foreground_sigma" : ["float", None,
                                  """Foreground filter sigma, this is the sigma of a 2D gaussian to convolve the data with
                                  prior to peak identification."""],

            "roi_size" : ["int", None,
                          """This is the size of the fitting ROI in pixels. If it is not specified then the value
                          will be calculated based on the sigma parameter value."""],

            "w_vs_z_params_1" : ["float-array", None, pmd_ww_doc_string],
            "w_vs_z_params_2" : ["float-array", None, pmd_ww_doc_string],
            "w_vs_z_params_3" : ["float-array", None, pmd_ww_doc_string],
            "w_vs_z_params_4" : ["float-array", None, pmd_ww_doc_string],
            "w_vs_z_params_5" : ["float-array", None, pmd_ww_doc_string],
            "w_vs_z_params_6" : ["float-array", None, pmd_ww_doc_string],
            "w_vs_z_params_7" : ["float-array", None, pmd_ww_doc_string],
            "w_vs_z_params_8" : ["float-array", None, pmd_ww_doc_string],     
        })


camera_calibration_doc_string = """This file contains the sCMOS calibration data for the region of the camera
that the movie comes from. It consists of 4 numpy arrays, [offset, variance, gain, rqe],
each of which is the same size as a frame of the movie that is to be analyzed.
This can be generated for a camera using camera_calibration.py and (if it needs
to be resliced), reslice_calibration.py."""


class ParametersPSFFFT(ParametersFitters):
    """
    Parameters that are specific to PSF FFT analysis.

    Note: The XML file should have either the 'camera_calibration' file for sCMOS analysis
          or 'camera_gain' and 'camera_offset', but not both.
    """
    def __init__(self, **kwds):
        super(ParametersPSFFFT, self).__init__(**kwds)

        self.attr.update({

            "camera_calibration" : ["filename", None, camera_calibration_doc_string],
            
            "camera_gain" : ["float", None, camera_gain_doc_string],

            "camera_offset" : ["float", None, camera_offset_doc_string],

            "psf" : ["filename", None,
                     "This is the psf file to use for fitting."],

            "z_value" : ["float-array", None,
                         """Z value(s) in microns at which we will perform convolution with the PSF for
                         the purposes of peak finding. If this is not specified the default value is
                         z = [0.0]. These are also the starting z values for fitting.
                         
                         See note in ParametersSplinerSTD for this parameter."""],
        })
        

class ParametersPupilFn(ParametersFitters):
    """
    Parameters that are specific to Pupil function analysis.

    Note: The XML file should have either the 'camera_calibration' file for sCMOS analysis
          or 'camera_gain' and 'camera_offset', but not both.
    """
    def __init__(self, **kwds):
        super(ParametersPupilFn, self).__init__(**kwds)

        self.attr.update({

            "camera_calibration" : ["filename", None, camera_calibration_doc_string],
            
            "camera_gain" : ["float", None, camera_gain_doc_string],

            "camera_offset" : ["float", None, camera_offset_doc_string],

            "pupil_function" : ["filename", None,
                                "This is the pupil function file to use for fitting."],

            "z_value" : ["float-array", None,
                         """Z value(s) in microns at which we will perform convolution with the PSF for
                         the purposes of peak finding. If this is not specified the default value is
                         z = [0.0]. These are also the starting z values for fitting.
                         
                         See note in ParametersSplinerSTD for this parameter."""],
        })


class ParametersRollingBall(Parameters):
    """
    Parameters for Rolling Ball background estimation.
    """
    def __init__(self, **kwds):
        super(ParametersRollingBall, self).__init__(**kwds)

        self.attr.update({
            "rb_radius" : ["float", None,
                           "Radius of the rolling ball in pixels."],
            
            "rb_sigma" : ["float", None,
                          """Sigma in pixels of the gaussian smoothing to apply to the background estimate after
                          the rolling ball step."""],
            })

        
class ParametersSCMOS(ParametersDAOsCMOS):
    """
    Parameters that are specific to sCMOS analysis.
    """
    def __init__(self, **kwds):
        super(ParametersSCMOS, self).__init__(**kwds)

        self.attr.update({

            "camera_calibration" : ["filename", None, camera_calibration_doc_string],

            })
    

class ParametersSpliner(ParametersFitters):
    """
    Parameters that are specific to Spliner analysis.

    Note: The XML file should have either the 'camera_calibration' file for sCMOS analysis
          or 'camera_gain' and 'camera_offset', but not both.
    """
    def __init__(self, **kwds):
        super(ParametersSpliner, self).__init__(**kwds)

        self.attr.update({

            "camera_calibration" : ["filename", None, camera_calibration_doc_string],
            
            "camera_gain" : ["float", None, camera_gain_doc_string],

            "camera_offset" : ["float", None, camera_offset_doc_string],

            "decon_method" : ["string", None,
                              """Use a compressed sensing deconvolution method for peak finding. If this is
                              not specified then peaks are identified by convolving the image with the PSF
                              at one or more z values. Possible values are 'FISTA', 'ADMM' and '3denseSTORM'"""],

            "spline" : ["filename", None,
                        """This is the spline file to use for fitting. Based on the spline the analysis will
                        decide whether to do 2D or 3D spline fitting, 2D if the spline is 2D, 3D if the
                        spline is 3D."""],
        })


class ParametersSplinerSTD(ParametersSpliner):
    """
    Parameters that are specific to Spliner standard analysis. In this
    approach peaks are identified by convolving the image with the PSF
    at one or more z values.
    """
    def __init__(self, **kwds):
        super(ParametersSplinerSTD, self).__init__(**kwds)

        self.attr.update({
            
            "z_value" : ["float-array", None,
                         """Z value(s) in microns at which we will perform convolution with the PSF for
                         the purposes of peak finding. If this is not specified the default value is
                         z = [0.0]. These are also the starting z values for fitting.
                         
                         If your PSF is not degenerate* in Z then it could be helpful to try multiple z
                         starting values. However most common 3D PSFs such as astigmatism do not meet
                         this criteria. The only one that does meet this criteria that is in (sort of)
                         common use is the double-helix PSF.
                         
                         * By degenerate I mean that the PSF at one z value can be modeled (with reasonable
                         \taccuracy) by summing several PSFs with a different z value. For example, most
                         \tastigmatic PSFs z != 0 can be modeled by summing several z = 0 PSFs with
                         \tvariable x,y positions."""],                 
        })
        

class ParametersWaveletBGR(Parameters):
    """
    Parameters for Wavelet background estimation.
    """
    def __init__(self, **kwds):
        super(ParametersWaveletBGR, self).__init__(**kwds)

        self.attr.update({
            "wbgr_iterations" : ["int", None,
                                 """The number of iterations of background estimation and foreground replacement to
                                 perform (see the Galloway paper), usually something like 2."""],

            "wbgr_threshold" : ["float", None,
                                """This is the difference between the current estimate and the signal at which the
                                signal we be considered "foreground". This should probably be something like 1x
                                to 2x the estimated noise in the background."""],

            "wbgr_wavelet_level" : ["int", None,
                                    """How many levels of wavelet decomposition to perform. The larger the number the less
                                    response to local changes in the background, usually something like 2."""],  
        })
        

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
