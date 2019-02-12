#!/usr/bin/env python
"""
Collate analysis results.

Hazen 01/18
"""
import math
import numpy

import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.sa_utilities.finding_fitting_error as ffe
import storm_analysis.sa_utilities.recall_fraction as rfrac


def collateDAO(dirs, settings, calc_width_error = True):
    """
    Results collations for 3D-DAOSTORM and sCMOS analysis.
    """
    all_dx = []
    all_dy = []
    all_wx = []
    all_wy = []
    noise = 0
    noise_total = 0
    recall = 0
    recall_total = 0
    total_locs = 0
    total_time = 0.0
    for a_dir in dirs:
        print("Processing", a_dir)

        # Load timing information.
        with open(a_dir + "/timing.txt") as fp:
            total_time += float(fp.readline())
            
        # Load localizations.
        truth_h5 = saH5Py.SAH5Py(a_dir + "/test_ref.hdf5")
        measured_h5 = saH5Py.SAH5Py(a_dir + "/test.hdf5")
        total_locs += measured_h5.getNLocalizations()
    
        # Calculate fractional recall.
        [partial, total] = rfrac.recallFraction(truth_h5, measured_h5, settings.tolerance)
        recall += partial
        recall_total += total

        # Calculate noise fraction.
        [partial, total] = rfrac.noiseFraction(truth_h5, measured_h5, settings.tolerance)
        noise += partial
        noise_total += total

        # Calculate error in fitting width.
        if calc_width_error:

            for i in range(truth_h5.getMovieLength()):    
                t_locs = truth_h5.getLocalizationsInFrame(i)
                m_locs = measured_h5.getLocalizationsInFrame(i)

                # Widths for truth localizations.
                t_wx = t_locs["xsigma"]
                t_wy = t_locs["ysigma"]
    
                # Widths for found localizations.
                if ("xsigma" in m_locs):
                    m_wx = m_locs["xsigma"]
                else:
                    m_wx = 1.5*numpy.ones(m_locs['x'].size)
                                      
                if ("ysigma" in m_locs):
                    m_wy = m_locs["ysigma"]
                else:
                    m_wy = m_wx.copy()

                p_index = iaUtilsC.peakToPeakDistAndIndex(m_locs['x'], m_locs['y'],
                                                          t_locs['x'], t_locs['y'],
                                                          max_distance = settings.tolerance)[1]

                p_size = numpy.count_nonzero(p_index > -1)
                d_wx = numpy.zeros(p_size)
                d_wy = numpy.zeros(p_size)
                k = 0
                for j in range(m_locs["x"].size):
                    if(p_index[j] < 0):
                        continue
            
                    d_wx[k] = m_wx[j] - t_wx[p_index[j]]
                    d_wy[k] = m_wy[j] - t_wy[p_index[j]]
                    k += 1

                all_wx.append(d_wx)
                all_wy.append(d_wy)

        # Calculate fitting error in XY.
        max_distance = None
        if True:
            max_distance = 2.0 * settings.pixel_size
            print("Using max_distance", max_distance, "nm for error calcuations.")
        
        [dx, dy, dz] = ffe.findingFittingError(truth_h5,
                                               measured_h5,
                                               pixel_size = settings.pixel_size,
                                               max_distance = max_distance)
        if dx.size != 0:
            all_dx.append([numpy.std(dx), math.sqrt(numpy.mean(dx*dx))])
            all_dy.append([numpy.std(dy), math.sqrt(numpy.mean(dy*dy))])
        else:
            all_dx.append([0,0])
            all_dy.append([0,0])

        truth_h5.close()
        measured_h5.close()


    print()
    print("Analysis Summary:")
    print("Processed {0:0d} localizations in {1:.2f} seconds, {2:.2f}/sec".format(total_locs, total_time, float(total_locs)/float(total_time)))
    print("Recall {0:.5f}".format(float(recall)/float(recall_total)))
    print("Noise {0:.5f}".format(float(noise)/float(noise_total)))
    print("XY Error Standard Deviation (nm):")
    for i, a_dir in enumerate(dirs):
        print(a_dir + "\t{0:.2f}\t{1:.2f}".format(all_dx[i][0], all_dy[i][0]))
    print("")
    print("XY RMSE (nm):")
    for i, a_dir in enumerate(dirs):
        print(a_dir + "\t{0:.2f}\t{1:.2f}".format(all_dx[i][1], all_dy[i][1]))
        
    if calc_width_error:
        print("")
        print("XY Width Error, Mean difference with truth, Standard deviation (pixels):")
        for i, a_dir in enumerate(dirs):
            print(a_dir + "\t{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}".format(numpy.mean(all_wx[i]),
                                                                        numpy.std(all_wx[i]),
                                                                        numpy.mean(all_wy[i]),
                                                                        numpy.std(all_wy[i])))


def collateSpliner(dirs, settings):
    """
    Results collation for Spliner, Pupil-Function, PSF-FFT and Multiplane.

    Note: The maximum distance cutoff is in XYZ, so points with Z values
          that are way off will not be included in the error calculation
          even though their XY values might be very good.
    """
    all_dx = []
    all_dy = []
    all_dz = []
    noise = 0
    noise_total = 0
    recall = 0
    recall_total = 0
    total_locs = 0
    total_time = 0.0
    for a_dir in dirs:
        print("Processing", a_dir)

        # Load timing information.
        with open(a_dir + "/timing.txt") as fp:
            total_time += float(fp.readline())
        
        # Load localizations.
        truth_h5 = saH5Py.SAH5Py(a_dir + "/test_ref.hdf5")
        measured_h5 = saH5Py.SAH5Py(a_dir + "/test.hdf5")    
        total_locs += measured_h5.getNLocalizations()
    
        # Calculate fractional recall.
        [partial, total] = rfrac.recallFraction(truth_h5, measured_h5, settings.tolerance)
        recall += partial
        recall_total += total

        # Calculate noise fraction.
        [partial, total] = rfrac.noiseFraction(truth_h5, measured_h5, settings.tolerance)
        noise += partial
        noise_total += total

        # Calculate fitting error in XYZ.
        max_distance = None
        if True:
            max_distance = 2.0 * settings.pixel_size
            print("Using max_distance", max_distance, "nm for error calcuations.")
        
        [dx, dy, dz] = ffe.findingFittingError(truth_h5,
                                               measured_h5,
                                               pixel_size = settings.pixel_size,
                                               max_distance = max_distance)
    
        if dx.size != 0:
            all_dx.append([numpy.std(dx), math.sqrt(numpy.mean(dx*dx))])
            all_dy.append([numpy.std(dy), math.sqrt(numpy.mean(dy*dy))])
            all_dz.append([numpy.std(dz), math.sqrt(numpy.mean(dz*dz))])
        else:
            all_dx.append([0,0])
            all_dy.append([0,0])
            all_dz.append([0,0])


    print()
    print("Analysis Summary:")
    print("Processed {0:0d} localizations in {1:.2f} seconds, {2:.2f}/sec".format(total_locs, total_time, float(total_locs)/float(total_time)))
    print("Recall {0:.5f}".format(float(recall)/float(recall_total)))
    print("Noise {0:.5f}".format(float(noise)/float(noise_total)))
    print("XYZ Error Standard Deviation (nm):")
    for i, a_dir in enumerate(dirs):
        print(a_dir + "\t{0:.2f}\t{1:.2f}\t{2:.2f}".format(all_dx[i][0], all_dy[i][0], all_dz[i][0]))
    print("")
    print("XYZ RMSE Accuracy (nm):")
    for i, a_dir in enumerate(dirs):
        print(a_dir + "\t{0:.2f}\t{1:.2f}\t{2:.2f}".format(all_dx[i][1], all_dy[i][1], all_dz[i][1]))
    print("")
