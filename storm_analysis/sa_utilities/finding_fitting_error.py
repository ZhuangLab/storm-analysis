#!/usr/bin/env python
"""
Calculate the error (for simulations) between the known
peak locations and the peak locations returned by the analysis.

Note: This measures the error and *not* the recall. If only 1
in 100 peaks is found but their locations agree exactly with
the known locations then this will be considered good.

One assumption here is that neither of these files has a
particularly large number of localizations.

Hazen 01/16
"""
import numpy

import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py


def findingFittingError(truth_h5, measured_h5, pixel_size = None, max_distance = None):
    """
    truth_h5 - A saH5Py.SAH5Py object with the ground truth localizations.
    measured_h5 - A saH5Py.SAH5Py object with the found localizations.
    pixel_size - The camera pixel size in nanometers. If not specified then the value
                 in the measured HDF5 file will be used.
    max_distance - If not none, found peaks that are greater than this distance from
                   a truth peak will be ignored. Units are nanometers.
    """
    if (measured_h5.getNLocalizations() == 0):
        return [None, None, None]

    md_in_pixels = None
    md_sqr = None
    if max_distance is not None:
        md_in_pixels = max_distance/pixel_size
        md_sqr = max_distance * max_distance

    if pixel_size is None:
        pixel_size = measured_h5.getPixelSize()
        
    all_dx = []
    all_dy = []
    all_dz = []
    for i in range(truth_h5.getMovieLength()):
        t_locs = truth_h5.getLocalizationsInFrame(i)
        m_locs = measured_h5.getLocalizationsInFrame(i)

        if not bool(t_locs) or not bool(m_locs):
            continue
        
        p_index = iaUtilsC.peakToPeakDistAndIndex(m_locs['x'], m_locs['y'],
                                                  t_locs['x'], t_locs['y'],
                                                  max_distance = md_in_pixels)[1]
        for i in range(m_locs['x'].size):
            if(p_index[i] < 0):
                continue
            dx = pixel_size * (m_locs['x'][i] - t_locs['x'][p_index[i]])
            dy = pixel_size * (m_locs['y'][i] - t_locs['y'][p_index[i]])
            
            if 'z' in m_locs:
                dz = 1000.0 * (m_locs['z'][i] - t_locs['z'][p_index[i]])
            else:
                dz = 0.0
                
            if md_sqr is not None:
                if ((dx*dx + dy*dy + dz*dz) < md_sqr):
                    all_dx.append(dx)
                    all_dy.append(dy)
                    all_dz.append(dz)
            else:
                all_dx.append(dx)
                all_dy.append(dy)
                all_dz.append(dz)

    return [numpy.array(all_dx), numpy.array(all_dy), numpy.array(all_dz)]

def findingFittingErrorHDF5File(truth_name, measured_name, pixel_size = None, max_distance = None):
    """
    A wrapper for findingFittingError to make it easier to use with HDF5 files.

    truth_name - The name of the HDF5 file containing the ground truth localizations.
    measured_name - The name of the HDF5 file containing the measured localizations.
    pixel_size - The camera pixel size in nanometers. If not specified then the value
                 in the measured HDF5 file will be used.
    max_distance - If not none, found peaks that are greater than this distance from
                   a truth peak will be ignored. Units are nanometers.
    """
    truth_h5 = saH5Py.SAH5Py(truth_name)
    measured_h5 = saH5Py.SAH5Py(measured_name)
    data = findingFittingError(truth_h5,
                               measured_h5,
                               pixel_size = pixel_size,
                               max_distance = max_distance)
    truth_h5.close()
    measured_h5.close()
    return data


if (__name__ == "__main__"):
    import argparse
    import matplotlib
    import matplotlib.pyplot as pyplot

    import storm_analysis.sa_library.gaussfit as gaussfit

    parser = argparse.ArgumentParser(description = 'Measure finding/fitting error.')

    parser.add_argument('--truth_bin', dest='truth_bin', type=str, required=True,
                        help = "Ground truth localization file.")
    parser.add_argument('--measured_bin', dest='measured_bin', type=str, required=True,
                        help = "Measured localization file.")
    parser.add_argument('--pixel_size', dest='pixel_size', type=float, required=False, default = 160.0,
                        help = "Camera pixel size in nanometers.")

    args = parser.parse_args()

    [all_dx, all_dy, all_dz] = findingFittingErrorHDF5File(args.truth_bin,
                                                           args.measured_bin,                           
                                                           pixel_size = args.pixel_size)
    
    print("means and standard deviations (in nm):")
    print("mean, std (dx)", numpy.mean(all_dx), numpy.std(all_dx))
    print("mean, std (dy)", numpy.mean(all_dy), numpy.std(all_dy))
    print("mean, std (dz)", numpy.mean(all_dz), numpy.std(all_dz))
    print("")

    h_range_xy = 100.0
    h_range_z = 200.0
    [hist_dx, bins_xy] = numpy.histogram(all_dx, bins = 30, range = (-h_range_xy, h_range_xy))
    [hist_dy, bins_xy] = numpy.histogram(all_dy, bins = 30, range = (-h_range_xy, h_range_xy))
    [hist_dz, bins_z] = numpy.histogram(all_dz, bins = 30, range = (-h_range_z, h_range_z))

    hist_dx = hist_dx.astype(numpy.float)/numpy.sum(hist_dx)
    hist_dy = hist_dy.astype(numpy.float)/numpy.sum(hist_dy)
    hist_dz = hist_dz.astype(numpy.float)/numpy.sum(hist_dz)

    centers_xy = bins_xy[:-1] + 0.5 * (bins_xy[1] - bins_xy[0])
    centers_z = bins_z[:-1] + 0.5 * (bins_z[1] - bins_z[0])

    print("gaussian fitting")
    bin_size_xy = bins_xy[1] - bins_xy[0]
    bin_size_z = bins_z[1] - bins_z[0]
    [fitx, goodx] =  gaussfit.fitSymmetricGaussian1D(hist_dx)
    [fity, goody] =  gaussfit.fitSymmetricGaussian1D(hist_dy)
    [fitz, goodz] =  gaussfit.fitSymmetricGaussian1D(hist_dz)

    print("")
    print("gaussian fit to error histogram width (in nm):")
    if goodx:
        print("x width", fitx[3]*bin_size_xy)
    if goody:
        print("y width", fity[3]*bin_size_xy)
    if goodz:
        print("z width", fitz[3]*bin_size_z)
    
    fig = pyplot.figure()
    pyplot.plot(centers_xy, hist_dx, color = "red")
    pyplot.plot(centers_xy, hist_dy, color = "green")
    pyplot.plot(centers_z, hist_dz, color = "blue")
    pyplot.xlabel("Error in nm")
    pyplot.ylabel("Density (AU)")
    pyplot.show()

    
#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
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
