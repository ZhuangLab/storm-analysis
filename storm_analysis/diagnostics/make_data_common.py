#!/usr/bin/env python
"""
That which is common to all of the XYZ.make_data diagnostic modules.

Hazen 06/19
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py


def makePeakFile(settings):
    # Create "peak_locations" file if needed.
    #
    if hasattr(settings, "peak_locations") and (settings.peak_locations is not None):
        with saH5Py.SAH5Py("test_01/test_ref.hdf5") as h5:
            locs = h5.getLocalizationsInFrame(0)
        
        if settings.peak_locations.endswith(".hdf5"):
            saH5Py.saveLocalizations(settings.peak_locations, locs)
        else:
            numpy.savetxt(settings.peak_locations,
                          numpy.transpose(numpy.vstack((locs['x'],
                                                        locs['y'],
                                                        locs['height'],
                                                        locs['background']))))
