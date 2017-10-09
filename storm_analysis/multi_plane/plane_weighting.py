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
import numpy
import pickle

import storm_analysis.multi_plane.mp_utilities_c as mpUtilC

import storm_analysis.sa_library.parameters as params

import storm_analysis.spliner.cramer_rao as cramerRao


def planeVariances(background, photons, pixel_size, spline_file_names):
    """
    Calculates the variances for different image planes as a function of z.

    Notes: 
       1. For multiplane analysis with a single peak height parameter
          (independent_heights = 0), the expectation is that the splines 
          are properly normalized, i.e. if one image plane receives less 
          photons than another plane this is already included in the spline.

       2. The background parameter is the average estimated background
          for each plane (in photons). If only one value is specified
          it will be used for all of the planes. The photons parameter 
          is the total number of photons in all of the planes.
    """
    n_planes = len(spline_file_names)
    
    photons = photons/n_planes

    #
    # If only one background value was supplied assume it is the same
    # for all of the planes.
    #
    if (len(background) == 1):
        for i in range(n_planes-1):
            background.append(background[0])
        
    # Create 3D Cramer-Rao bound objects.
    #
    # Results for each plane are weighted by the contribution of
    # the plane to the overall signal.
    #
    CRB3Ds = []
    for name in spline_file_names:
        with open(name, 'rb') as fp:
            weighting = pickle.load(fp)["maximum"]
        CRB3Ds.append(cramerRao.CRBound3D(name, pixel_size, weighting = weighting))

    #
    # Calculate Cramer-Rao bounds at each (integer) z
    # position in the spline. This is the granularity that
    # we are using for multi-plane analysis.
    #
    # FIXME: Is this granular enough?
    #
    n_zvals = CRB3Ds[0].getSize()
    
    v_bg = numpy.zeros((n_zvals, n_planes))
    v_h = numpy.zeros((n_zvals, n_planes))
    v_x = numpy.zeros((n_zvals, n_planes))
    v_y = numpy.zeros((n_zvals, n_planes))
    v_z = numpy.zeros((n_zvals, n_planes))

    for i in range(n_zvals):
        print("scaled z", i, ", max", n_zvals - 1)
        for j in range(n_planes):
            crbs = CRB3Ds[j].calcCRBoundScaledZ(background[j], photons, i)
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

def xyMappingCorrect(mapping_file, var_x, var_y):
    """
    Apply the mapping to the weights so that they are all correct
    in the channel 0 frame, which is what the analysis uses.
    """
    with open(mapping_file, 'rb') as fp:
        mappings = pickle.load(fp)

    vx_corr = numpy.zeros(var_x.shape)
    vy_corr = numpy.zeros(var_y.shape)
    for i in range(var_x.shape[0]):
        for j in range(var_x.shape[1]):
            xv = var_x[i,j]
            yv = var_y[i,j]

            # Correct X.
            mx = mappings[str(j) + "_0_x"]            
            mag = 1.0/math.sqrt(mx[1]*mx[1] + mx[2]*mx[2])
            vx_corr[i,j] = mx[1]*xv*mag + mx[2]*yv*mag

            # Correct Y.
            my = mappings[str(j) + "_0_y"]
            mag = 1.0/math.sqrt(my[1]*my[1] + my[2]*my[2])
            vy_corr[i,j] = my[1]*xv*mag + my[2]*yv*mag

    return [vx_corr, vy_corr]

    
if (__name__ == "__main__"):

    import argparse

    import matplotlib
    import matplotlib.pyplot as pyplot

    parser = argparse.ArgumentParser(description = 'Calculates how to weight the different image planes.')
    
    parser.add_argument('--background', dest='background', type=int, required=True, nargs = "*",
                        help = "An estimate of the image background in photons.")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the file to save the weights in.")
    parser.add_argument('--photons', dest='photons', type=int, required=True,
                        help = "An estimate of the average number of photons in the localization.")
    parser.add_argument('--xml', dest='xml', type=str, required=True,
                        help = "The name of the settings xml file.")
    parser.add_argument('--no_plots', dest='no_plots', type=bool, required=False, default=False,
                        help = "Don't show plot of the results.")
    
    args = parser.parse_args()

    parameters = params.ParametersMultiplane().initFromFile(args.xml)
    spline_file_names = []
    for spline_attr in mpUtilC.getSplineAttrs(parameters):
        spline_file_names.append(parameters.getAttr(spline_attr))

    print("Calculating Cramer-Rao bounds.")
    variances = planeVariances(args.background,
                               args.photons,
                               parameters.getAttr("pixel_size"),
                               spline_file_names)

    #
    # FIXME: In the C library x and y are transposed, so we are 
    #        transposing the variances here. I'm pretty sure that
    #        this is correct it has not been formally tested.
    #
    print("Correcting for mapping.")
    [variances[3], variances[2]] = xyMappingCorrect(parameters.getAttr("mapping"),
                                                    variances[2],
                                                    variances[3])
    
    weights = {"bg" : planeWeights(variances[0]),
               "h" : planeWeights(variances[1]),
               "x" : planeWeights(variances[2]),
               "y" : planeWeights(variances[3]),
               "z" : planeWeights(variances[4]),
               "background" : args.background,
               "photons" : args.photons}

    with open(args.output, 'wb') as fp:
        pickle.dump(weights, fp)

    #
    # Plot results.
    #
    if not args.no_plots:
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
            pyplot.show()
