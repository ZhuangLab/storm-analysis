#!/usr/bin/env python
"""
Uses dbscan to cluster the data in an insight3
file. Molecules in different categories are
clustered separately. The results are saved in 
the lk field of the output molecule list.

Hazen 11/11
"""

import numpy

import storm_analysis.dbscan.dbscan_c as dbscanC
import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def findClusters(mlist_name, clist_name, eps, mc, ignore_z = True, ignore_category = True):
    
    # Load the data.
    pix_to_nm = 160.0

    i3_data_in = readinsight3.loadI3GoodOnly(mlist_name)

    c = i3_data_in['c']
    x = i3_data_in['xc']*pix_to_nm
    y = i3_data_in['yc']*pix_to_nm

    if ignore_z:
        print("Warning! Clustering without using localization z value!")
        z = numpy.zeros(x.size)
    else:
        z = i3_data_in['zc']

    # Perform analysis without regard to category.
    if ignore_category:
        print("Warning! Clustering without regard to category!")
        c = numpy.zeros(c.size)

    # Cluster the data.
    labels = dbscanC.dbscan(x, y, z, c, eps, mc, z_factor=1.0)

    # Save the data.    
    i3_data_out = writeinsight3.I3Writer(clist_name)
    i3dtype.setI3Field(i3_data_in, 'lk', labels)
    i3_data_out.addMolecules(i3_data_in)
    i3_data_out.close()


