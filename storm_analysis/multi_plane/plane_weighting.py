#!/usr/bin/env python
"""
Uses the Cramer-Rao bound formalism to determine how to best 
weight the updates from each image plane. 

Ideally perhaps we'd calculate this in the C library for
each peak at each update given it's current parameters, but
I think this would make things even slower for only a fairly
minor improvement in fitting performance.

Hazen 10/17
"""
import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy
import pickle

import storm_analysis.multi_plane.mp_utilities as mpUtil

import storm_analysis.sa_library.parameters as params

import storm_analysis.psf_fft.cramer_rao as psfFFTCramerRao
import storm_analysis.pupilfn.cramer_rao as pupilFnCramerRao
import storm_analysis.spliner.cramer_rao as splinerCramerRao


def planeVariances(cr_psf_objects, background, photons, verbose = True):
    """
    Calculates the variances for different image planes as a function of z.

    Notes: 
       1. For multiplane analysis with a single peak height parameter
          (independent_heights = 0), the expectation is that the psf objects
          are properly normalized, i.e. if one image plane receives less 
          photons than another plane this is already included in the psf object.

       2. The background parameter is the average estimated background
          for each plane (in photons). If only one value is specified
          it will be used for all of the planes. The photons parameter 
          is the total number of photons in all of the planes.
    """
    n_planes = len(cr_psf_objects)
    
    photons = photons/n_planes

    # If only one background value was supplied assume it is the same
    # for all of the planes.
    #
    if (len(background) == 1):
        for i in range(n_planes-1):
            background.append(background[0])

    # Calculate Cramer-Rao bounds at different z positions across
    # the range covered by the CRPSFObject.
    #
    z_min = cr_psf_objects[0].getZMin()
    z_max = cr_psf_objects[0].getZMax()
    n_zvals = cr_psf_objects[0].getNZValues()
    step = (z_max - z_min)/float(n_zvals - 1.0)

    z_vals = numpy.arange(z_min, z_max + 0.5*step, step)
    assert (n_zvals == z_vals.size), "number of z values does not match."
    
    v_bg = numpy.zeros((n_zvals, n_planes))
    v_h = numpy.zeros((n_zvals, n_planes))
    v_x = numpy.zeros((n_zvals, n_planes))
    v_y = numpy.zeros((n_zvals, n_planes))
    v_z = numpy.zeros((n_zvals, n_planes))

    for i in range(z_vals.size):
        z = z_vals[i]
        if verbose:
            print("z {0:.1f}".format(z))
        for j in range(n_planes):
            crbs = splinerCramerRao.calcCRBound3D(cr_psf_objects[j],
                                                  background[j],
                                                  photons,
                                                  z * 1.0e-3)
            v_bg[i,j] = crbs[4]
            v_h[i,j] = crbs[0]
            v_x[i,j] = crbs[1]
            v_y[i,j] = crbs[2]
            v_z[i,j] = crbs[3]

    return [v_bg, v_h, v_x, v_y, v_z]

def planeWeights(variances):
    """
    Weights are 1/variance normalized to sum to unity.
    """
    weights = numpy.zeros(variances.shape)
    for i in range(weights.shape[0]):
        weights[i,:] = 1.0/variances[i,:]
        weights[i,:] = weights[i,:]/numpy.sum(weights[i,:])
    return weights

def planeWeighting(parameters, background, photons):
    """
    This calculates and return the weights to use for each parameter at each z value.
    """
    pixel_size = parameters.getAttr("pixel_size")
    cr_psf_objects = []


    # Try PSF FFTs.
    #
    if (len(mpUtil.getPSFFFTAttrs(parameters)) > 0):

        # Create PSF FFT CR PSF objects.
        for psf_fft_attr in mpUtil.getPSFFFTAttrs(parameters):
            psf_fft_name = parameters.getAttr(psf_fft_attr)
            cr_psf_objects.append(psfFFTCramerRao.CRPSFFn(psf_filename = psf_fft_name,
                                                          pixel_size = pixel_size))

    # Try pupil functions.
    #
    if (len(mpUtil.getPupilFnAttrs(parameters)) > 0):

        # Create pupil function CR PSF objects.
        [zmin, zmax] = parameters.getZRange()
        for pfn_attr in mpUtil.getPupilFnAttrs(parameters):
            pfn_name = parameters.getAttr(pfn_attr)
            cr_psf_objects.append(pupilFnCramerRao.CRPupilFn(psf_filename = pfn_name,
                                                             pixel_size = pixel_size,
                                                             zmax = zmax,
                                                             zmin = zmin))
    
    # Try splines.
    #
    if (len(mpUtil.getSplineAttrs(parameters)) > 0):

        # Create Spline CR PSF objects.
        for spline_attr in mpUtil.getSplineAttrs(parameters):
            spline_name = parameters.getAttr(spline_attr)
            cr_psf_objects.append(splinerCramerRao.CRSplineToPSF3D(psf_filename = spline_name,
                                                                   pixel_size = pixel_size))
        

    assert (len(cr_psf_objects) > 0), "No CR PSF objects were found."
    
    print("Calculating Cramer-Rao bounds.")
    variances = planeVariances(cr_psf_objects,
                               background,
                               photons)

    # Clean up Cramer-Rao PSF objects (some use C libraries).
    #
    for cr_po in cr_psf_objects:
        cr_po.cleanup()

    # I believe based on simulation testing that X and Y should not be transposed here.
    #
    print("Correcting for mapping.")
    [variances[2], variances[3]] = xyMappingCorrect(parameters.getAttr("mapping"),
                                                    variances[2],
                                                    variances[3])
    
    weights = {"bg" : planeWeights(variances[0]),
               "h" : planeWeights(variances[1]),
               "x" : planeWeights(variances[2]),
               "y" : planeWeights(variances[3]),
               "z" : planeWeights(variances[4]),
               "background" : background,
               "photons" : photons}

    return [weights, variances]

def runPlaneWeighting(xml, output, background, photons, no_plots = False):
    """
    xml - The analysis XML file.
    output - File name to save the weights in.
    background - Per pixel background in e- for each plane (as a list).
    photons - Integrated peak intensity in e-.
    no_plots - Don't show any plots.
    """
    parameters = params.ParametersMultiplaneArb().initFromFile(xml)
    [weights, variances] = planeWeighting(parameters, background, photons)

    with open(output, 'wb') as fp:
        pickle.dump(weights, fp)

    #
    # Plot results.
    #
    if not no_plots:
        for i, name in enumerate(["bg", "h", "x", "y", "z"]):

            # Plot per channel standard deviation.
            sd = numpy.sqrt(variances[i])
            x = numpy.arange(sd.shape[0])
            fig = pyplot.figure()
            for j in range(sd.shape[1]):
                pyplot.plot(x, sd[:,j])

            sd = 1.0/numpy.sqrt(numpy.sum(1.0/variances[i], axis = 1))
            pyplot.plot(x, sd, color = "black")
                
            pyplot.title(name)
            pyplot.xlabel("scaled z")
            pyplot.ylabel("standard deviation (nm)")
            pyplot.ylim(ymin=0)
            pyplot.show()


def xyMappingCorrect(mapping_file, var_x, var_y):
    """
    Apply the mapping to the weights so that they are all correct
    in the channel 0 frame, which is what the analysis uses.
    """
    with open(mapping_file, 'rb') as fp:
        mappings = pickle.load(fp)

    if not "0_0_x" in mappings:
        mappings["0_0_x"] = numpy.array([0.0, 1.0, 0.0])
        mappings["0_0_y"] = numpy.array([0.0, 0.0, 1.0])
        
    vx_corr = numpy.zeros(var_x.shape)
    vy_corr = numpy.zeros(var_y.shape)
    for i in range(var_x.shape[0]):
        for j in range(var_x.shape[1]):
            xv = var_x[i,j]
            yv = var_y[i,j]

            # Correct X.
            mx = mappings[str(j) + "_0_x"]            
            mag = 1.0/math.sqrt(mx[1]*mx[1] + mx[2]*mx[2])
            vx_corr[i,j] = abs(mx[1])*xv*mag + abs(mx[2])*yv*mag

            # Correct Y.
            my = mappings[str(j) + "_0_y"]
            mag = 1.0/math.sqrt(my[1]*my[1] + my[2]*my[2])
            vy_corr[i,j] = abs(my[1])*xv*mag + abs(my[2])*yv*mag

    return [vx_corr, vy_corr]

    
if (__name__ == "__main__"):

    import argparse


    parser = argparse.ArgumentParser(description = 'Calculates how to weight the different image planes.')
    
    parser.add_argument('--background', dest='background', type=int, required=True, nargs = "*",
                        help = "An estimate of the image background in photons.")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the file to save the weights in.")
    parser.add_argument('--photons', dest='photons', type=int, required=True,
                        help = "An estimate of the average number of photons in the localization.")
    parser.add_argument('--xml', dest='xml', type=str, required=True,
                        help = "The name of the settings xml file.")
    parser.add_argument('--no_plots', dest='no_plots', action='store_true', default=False,
                        help = "Don't show a plot of the results.")
    
    args = parser.parse_args()

    runPlaneWeighting(args.xml,
                      args.output,
                      args.background,
                      args.photons,
                      no_plots = args.no_plots)

