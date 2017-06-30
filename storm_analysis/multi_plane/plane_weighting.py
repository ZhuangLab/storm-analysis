#!/usr/bin/env python
"""
Uses the Cramer-Rao bound formalism to determine how
to best weight the updates from each image plane.

Hazen 06/17
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
       1. The expectation is that the splines are properly normalized,
          i.e. if one image plane receives less photons than another 
          plane this is already included in the spline.

       2. The background parameter is the average estimated background
          for each plane (in photons). The photons parameter is the
          total number of photons in all of the planes.
    """
    n_planes = len(spline_file_names)
    
    photons = photons/n_planes

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
    n_zvals = CRB3Ds[0].getSize()
    
    v_bg = numpy.zeros((n_zvals, n_planes))
    v_h = numpy.zeros((n_zvals, n_planes))
    v_x = numpy.zeros((n_zvals, n_planes))
    v_y = numpy.zeros((n_zvals, n_planes))
    v_z = numpy.zeros((n_zvals, n_planes))

    for i in range(n_zvals):
        print("scaled z", i, ", max", n_zvals - 1)
        for j in range(n_planes):
            crbs = CRB3Ds[j].calcCRBoundScaledZ(background, photons, i)
            v_bg[i,j] = crbs[4]
            v_h[i,j] = crbs[0]
            v_x[i,j] = crbs[1]
            v_y[i,j] = crbs[2]
            v_z[i,j] = crbs[3]

    return [v_bg, v_h, v_x, v_y, v_z]

def planeWeights(variances):
    weights = numpy.zeros(variances.shape)
    for i in range(weights.shape[0]):
        weights[i,:] = 1.0/variances[i,:]
        weights[i,:] = weights[i,:]/numpy.sum(weights[i,:])
    return weights


if (__name__ == "__main__"):

    import argparse

    import matplotlib
    import matplotlib.pyplot as pyplot

    parser = argparse.ArgumentParser(description = 'Calculates how to weight the different image planes.')
    
    parser.add_argument('--background', dest='background', type=int, required=True,
                        help = "The image background in photons.")
    parser.add_argument('--photons', dest='photons', type=int, required=True,
                        help = "The number of photons in the localization.")
    parser.add_argument('--xml', dest='xml', type=str, required=True,
                        help = "The name of the settings xml file.")

    args = parser.parse_args()

    parameters = params.ParametersMultiplane().initFromFile(args.xml)
    spline_file_names = []
    for spline_attr in mpUtilC.getSplineAttrs(parameters):
        spline_file_names.append(parameters.getAttr(spline_attr))
        
    variances = planeVariances(args.background,
                               args.photons,
                               parameters.getAttr("pixel_size"),
                               spline_file_names)

    weights = list(map(planeWeights, variances))

    print(weights[0][0,:])
    print(variances[0][0,:])
    print(numpy.sqrt(variances[0][0,:]))
    print(1.0/numpy.sqrt(numpy.sum(1.0/variances[0], axis = 1))[0])


    #
    # Plot results.
    #
    if False:
        for i, name in enumerate(["bg", "h", "x", "y", "z"]):

            # Plot per channel standard deviation.
            sd = numpy.sqrt(variances[i])
            x = numpy.arange(sd.shape[0])
            fig = pyplot.figure()
            for i in range(sd.shape[1]):
                pyplot.plot(x, sd[:,i])

            sd = numpy.sum(numpy.sqrt(variances[i] * weights[i]), axis = 1)
            pyplot.plot(x, sd, color = "black")
                
            pyplot.title(name)
            pyplot.xlabel("scaled z")
            pyplot.ylabel("standard deviation (nm)")
            pyplot.show()
