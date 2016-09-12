#!/usr/bin/env python
#
# Uses dbscan to cluster the data in an insight3
# file. Molecules in different categories are
# clustered separately. The results are saved in 
# the lk field of the output molecule list.
#
# Hazen 11/11
#

import numpy
import sys

import storm_analysis.db_scan.dbscan_c as dbscanC
import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

# Load the data.
pix_to_nm = 160.0

i3_data_in = readinsight3.loadI3GoodOnly(sys.argv[1])

c = i3_data_in['c']
x = i3_data_in['xc']*pix_to_nm
y = i3_data_in['yc']*pix_to_nm

if False:
    z = i3_data_in['zc']
else:
    print "Warning! Clustering without using localization z value!"
    z = numpy.zeros(x.size)

# Perform analysis without regard to category.
if True:
    print "Warning! Clustering without regard to category!"
    c = numpy.zeros(c.size)

# Cluster the data.
if(len(sys.argv)==4):
    print "Using eps =", sys.argv[2], "mc =", sys.argv[3]
    labels = dbscanC.dbscan(x,y,z,c,float(sys.argv[2]),int(sys.argv[3]),z_factor=1.0)
else:
    print "Using eps = 80, mc = 5"
    labels = dbscanC.dbscan(x,y,z,c,80.0,5,z_factor=1.0)

# Save the data.
i3_data_out = writeinsight3.I3Writer(sys.argv[1][:-8] + "clusters_list.bin")
i3dtype.setI3Field(i3_data_in, 'lk', labels)
i3_data_out.addMolecules(i3_data_in)
i3_data_out.close()
