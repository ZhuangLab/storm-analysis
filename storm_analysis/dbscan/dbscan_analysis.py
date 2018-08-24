#!/usr/bin/env python
"""
DBSCAN cluster analysis functions.

Hazen 08/18
"""
import numpy
import os

import storm_analysis.dbscan.clusters_sa_h5py as clSAH5Py
import storm_analysis.dbscan.dbscan_c as dbscanC


def dbscanAnalysis(h5_name, eps = 40, mc = 10, min_size = 50):
    """
    Runs findClusters() and clusterStats() on a storm-analysis format HDF5
    format file.

    h5_name - The name of the HDF5 file.
    eps - DBSCAN epsilon parameter (in nanometers).
    mc - DBSCAN mc parameter.
    min_size - Minimum size cluster for cluster statistics.

    Note: This uses the default behavior for findClusters() which is to
          ignore the z position and category.
    """
    findClusters(h5_name, eps, mc)
    

def findClusters(h5_name, eps, mc, ignore_z = True, ignore_category = True, z_factor = 1.0):
    """
    Perform DBSCAN clustering on an HDF5 localization file.

    h5_name - The name of the HDF5 file.
    eps - DBSCAN epsilon parameter (in nanometers).
    mc - DBSCAN mc parameter.
    ignore_z - Ignore localization z position when clustering.
    ignore_category - Ignore localization category when clustering.
    z_factor - Weighting of Z versus X/Y position. A value of 0.5 for
               example will make the clustering 1/2 as sensitive to
               Z position.

    Note: Because all the x/y/z location information must be loaded
          into memory for the DBSCAN algorithm there is a limit to
          the size of localization file that can be clustered.
    """
    with clSAH5Py.SAH5Clusters(h5_name) as cl_h5:
        [x, y, z, c, cl_dict] = cl_h5.getDataForClustering(ignore_z = ignore_z)

        if ignore_z:
            print("Warning! Clustering without using localization z value!")        

        # Perform analysis without regard to category.
        if ignore_category:
            print("Warning! Clustering without regard to category!")
            c = numpy.zeros(c.size)

        # Cluster the data.
        labels = dbscanC.dbscan(x, y, z, c, eps, mc, z_factor = z_factor)

        # Save the data.
        cl_h5.addClusters(labels, cl_dict)

        # Save clustering info.
        info = "dbscan,eps,{0:0.3f},mc,{1:d}".format(eps,mc)
        info += ",iz," + str(ignore_z)
        info += ",ic," + str(ignore_category)
        info += ",zf,{0:3f}".format(z_factor)
        cl_h5.setClusteringInfo(info)
                

if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'DBSCAN clustering following Ester, KDD-96, 1996')

    parser.add_argument('--bin', dest='mlist', type=str, required=True,
                        help = "The name of the localizations input file. This is a HDF5 format localization file.")
    parser.add_argument('--eps', dest='epsilon', type=float, required=False, default=40,
                        help = "The DBSCAN epsilon parameters in nanometers. The default is 40nm.")
    parser.add_argument('--mc', dest='mc', type=int, required=False, default=10,
                        help = "The DBSCAN mc parameter. The default is 10.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=50,
                        help = "The minimum cluster size to include when calculating cluster statistics. The default is 50.")

    args = parser.parse_args()

    dbscanAnalysis(args.mlist, args.epsilon, args.mc, args.min_size)
