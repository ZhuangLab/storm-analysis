#!/usr/bin/env python
"""
DBSCAN cluster analysis functions.

Hazen 08/18
"""
import math
import numpy
import os

import storm_analysis.dbscan.clusters_sa_h5py as clSAH5Py
import storm_analysis.dbscan.dbscan_c as dbscanC


def clusterStats(h5_name, min_size, verbose = True):
    """
    Creates a text file containing some common cluster statistics.
    """
    with clSAH5Py.SAH5Clusters(h5_name) as cl_h5:
        pix_to_nm = cl_h5.getPixelSize()

        stats_name = os.path.splitext(h5_name)[0] + "_stats.txt"
        stats_fp = open(stats_name, "w")
        header = ["cluster", "cat", "size",
                  "x-center(nm)", "y-center(nm)", "z-center(nm)",
                  "size-x(nm)", "size-y(nm)", "size-z(nm)", "rg"]
        stats_fp.write(" ".join(header) + "\n")

        # Calculate cluster stats.
        for index, cluster in cl_h5.clustersIterator(min_size = min_size, fields = ["category", "x", "y", "z"]):
            c = cluster['category']
            x = pix_to_nm * cluster['x']
            y = pix_to_nm * cluster['y']
            
            if 'z' in cluster:
                z = 1000.0 * cluster['z']
            else:
                z = numpy.zeros(x.size)

            # Calculate size in x, y, z.
            sx = numpy.max(x) - numpy.min(x)
            sy = numpy.max(y) - numpy.min(y)
            sz = numpy.max(z) - numpy.min(z)

            # Calculate radius of gyration.
            cx = numpy.mean(x)
            cy = numpy.mean(y)

            rx = x - cx
            ry = y - cy

            # 3D radius of gyration if we have 'z' data.
            if 'z' in cluster:
                cz = numpy.mean(z)
                rz = z - cz
                rg = math.sqrt(numpy.sum(rx*rx + ry*ry + rz*rz) / float(x.size))

            # Otherwise 2D.
            else:
                rg = math.sqrt(numpy.sum(rx*rx + ry*ry) / float(x.size))

            if verbose:
                print("Cluster:", index, x.size, "localizations")
            stats = map(str, [index, c[0], x.size, numpy.mean(x), numpy.mean(y), numpy.mean(z), sx, sy, sz, rg])
            stats_fp.write(" ".join(stats) + "\n")
            
        stats_fp.close()

    return stats_name

        
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
    clusterStats(h5_name, min_size)
    

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
        [x, y, z, c, cl_dict] = cl_h5.getDataForClustering()

        if ignore_z:
            print("Warning! Clustering without using localization z value!")

        # Perform analysis without regard to category.
        if ignore_category:
            print("Warning! Clustering without regard to category!")
            c = numpy.zeros(c.size)

        # Convert data to nanometers
        pix_to_nm = cl_h5.getPixelSize()
        x_nm = x * pix_to_nm
        y_nm = y * pix_to_nm

        if ignore_z:
            z_nm = numpy.zeros(z.size)
        else:
            z_nm = 1000.0 * z
        
        # Cluster the data.
        labels = dbscanC.dbscan(x_nm, y_nm, z_nm, c, eps, mc, z_factor = z_factor)

        # Save the data. As an optimization we also save the x,y,z and
        # category values for each cluster with the cluster. Note that
        # these are the original units x/y in pixels and z in microns,
        # not the nanometer values used for clustering.
        #
        cl_dict["x"] = x
        cl_dict["y"] = y
        cl_dict["z"] = z
        cl_dict["category"] = c
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

    parser.add_argument('--bin', dest='sah5', type=str, required=True,
                        help = "The name of the localizations input file. This is a HDF5 format localization file.")
    parser.add_argument('--eps', dest='epsilon', type=float, required=False, default=40,
                        help = "The DBSCAN epsilon parameters in nanometers. The default is 40nm.")
    parser.add_argument('--mc', dest='mc', type=int, required=False, default=10,
                        help = "The DBSCAN mc parameter. The default is 10.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=50,
                        help = "The minimum cluster size to include when calculating cluster statistics. The default is 50.")

    args = parser.parse_args()

    dbscanAnalysis(args.sah5, args.epsilon, args.mc, args.min_size)
