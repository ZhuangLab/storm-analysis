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


def planeWeights(background, photons, pixel_size, spline_file_names):
    """
    Calculates how to weight the fit updates from different image planes
    as a function of z.

    Note: The expectation is that the splines are properly normalized,
          i.e. if one image plane receives less photons than another 
          plane this is already included in the spline.
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
    
    # If there is only one plane than the weighting is easy.
    if (n_planes == 1):
        w_bg = numpy.ones((n_zvals, n_planes))
        w_h = numpy.ones((n_zvals, n_planes))
        w_x = numpy.ones((n_zvals, n_planes))
        w_y = numpy.ones((n_zvals, n_planes))
        w_z = numpy.ones((n_zvals, n_planes))

    else:
        w_bg = numpy.zeros((n_zvals, n_planes))
        w_h = numpy.zeros((n_zvals, n_planes))
        w_x = numpy.zeros((n_zvals, n_planes))
        w_y = numpy.zeros((n_zvals, n_planes))
        w_z = numpy.zeros((n_zvals, n_planes))

        for i in range(n_zvals):
            print("z", i)
            for j in range(n_planes):
                crbs = CRB3Ds[j].calcCRBoundScaledZ(background, photons, i)
                w_bg[i,j] = math.sqrt(crbs[4])
                w_h[i,j] = math.sqrt(crbs[0])
                w_x[i,j] = math.sqrt(crbs[1])
                w_y[i,j] = math.sqrt(crbs[2])
                w_z[i,j] = math.sqrt(crbs[3])

    return [w_bg, w_h, w_x, w_y, w_z]


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
        
    [w_bg, w_h, w_x, w_y, w_z] = planeWeights(args.background,
                                              args.photons,
                                              parameters.getAttr("pixel_size"),
                                              spline_file_names)

    for elt in [["bg", w_bg], ["h", w_h], ["x", w_x], ["y", w_y], ["z", w_z]]:
        [name, ww] = elt
        x = numpy.arange(ww.shape[0])
        fig = pyplot.figure()
        for i in range(ww.shape[1]):
            pyplot.plot(x, ww[:,i])
        pyplot.title(name)
        pyplot.show()
