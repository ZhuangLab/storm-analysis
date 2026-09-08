#!/usr/bin/env python
"""
Does 2D Voroni cluster analysis.

Hazen 08/18
"""

import numpy
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import Polygon

import storm_analysis.dbscan.clusters_sa_h5py as clSAH5Py
import storm_analysis.dbscan.dbscan_analysis as dbscanAnalysis


def findClusters(h5_name, density_factor, min_size, verbose = True):
    """
    h5_name - The localizations HDF5 file.
    density_factor - Multiple of the median density to be a cluster member.
    min_size - The minimum number of localizations a cluster can have.
    """

    with clSAH5Py.SAH5Clusters(h5_name) as cl_h5:
        [x, y, z, c, cl_dict] = cl_h5.getDataForClustering()
        
        n_locs = x.size
        labels = numpy.zeros(n_locs, dtype = numpy.int32) - 1
        density = numpy.zeros(n_locs)
        
        # Convert data to nanometers
        pix_to_nm = cl_h5.getPixelSize()
        x_nm = x * pix_to_nm
        y_nm = y * pix_to_nm
        points = numpy.column_stack((x_nm, y_nm))

        if verbose:
            print("Creating Voronoi object.")
        vor = Voronoi(points)

        if verbose:
            print("Calculating 2D region sizes.")
        for i, region_index in enumerate(vor.point_region):
            if ((i%10000) == 0) and verbose:
                print("Processing point", i)

            vertices = []
            for vertex in vor.regions[region_index]:
        
                # I think these are edge regions?
                if (vertex == -1):
                    vertices = []
                    break

                vertices.append(vor.vertices[vertex])
            
            if (len(vertices) > 0):
                area = Polygon(vertices).area
                density[i] = 1.0/area

        # Used median density based threshold.
        ave_density = numpy.median(density)
        if verbose:
            print("Min density", numpy.amin(density))
            print("Max density", numpy.amax(density))
            print("Median density", ave_density)

        # Record the neighbors of each point. These are polygons so there shouldn't
        # be that many neighbors (sides). 40 is more than safe?
        #
        max_neighbors = 40
        neighbors = numpy.zeros((n_locs, max_neighbors), dtype = numpy.int32) - 1
        neighbors_counts = numpy.zeros((n_locs), dtype = numpy.int32)

        if verbose:
            print("Calculating neighbors")
        for ridge_p in vor.ridge_points:

            p1 = ridge_p[0]
            p2 = ridge_p[1]

            # Add p2 to the list for p1
            neighbors[p1,neighbors_counts[p1]] = p2
            neighbors_counts[p1] += 1

            # Add p1 to the list for p2
            neighbors[p2,neighbors_counts[p2]] = p1
            neighbors_counts[p2] += 1

        if False:
            n1 = neighbors[0,:]
            print(n1)
            print(neighbors[n1[0],:])

        # Mark connected points that meet the minimum density criteria.
        if verbose:
            print("Marking connected regions")
        min_density = density_factor * ave_density
        visited = numpy.zeros(n_locs, dtype = numpy.int32)

        def neighborsList(index):
            nlist = []
            for i in range(neighbors_counts[index]):
                loc_index = neighbors[index,i]
                if (visited[loc_index] == 0):
                    nlist.append(neighbors[index,i])
                    visited[loc_index] = 1
            return nlist

        cluster_id = 0
        for i in range(n_locs):
            if (visited[i] == 0):
                visited[i] = 1
                if (density[i] > min_density):
                    cluster_elt = [i]
                    c_size = 1
                    to_check = neighborsList(i)
                    while (len(to_check) > 0):

                        # Remove last localization from the list.
                        loc_index = to_check[-1]
                        to_check = to_check[:-1]

                        # If the localization has sufficient density add to cluster and
                        # check neighbors.
                        if (density[loc_index] > min_density):
                            to_check += neighborsList(loc_index)
                            cluster_elt.append(loc_index)
                            c_size += 1

                        # Mark as visited.
                        visited[loc_index] = 1

                    # Mark the cluster if there are enough localizations in the cluster.
                    if (c_size > min_size):
                        if verbose:
                            print("cluster", cluster_id, "size", c_size)
                        for elt in cluster_elt:
                            labels[elt] = cluster_id
                        cluster_id += 1

        if verbose:
            print(cluster_id, "clusters")

        # Save the clustering results.
        cl_dict["x"] = x
        cl_dict["y"] = y
        cl_dict["z"] = z
        cl_dict["density"] = density
        cl_dict["category"] = c
        cl_h5.addClusters(labels, cl_dict)

        # Save clustering info.
        info = "voronoi,df,{0:0.3f},ms,{1:d}".format(density_factor,min_size)
        cl_h5.setClusteringInfo(info)

    
def voronoiAnalysis(h5_name, density_factor, min_size = 30):
    """
    Runs voronoi() and clusterStats() on a storm-analysis format HDF5
    format file.

    h5_name - The name of the HDF5 file.
    density_factor - Voronoi clustering density factor.
    min_size - Minimum size cluster for cluster statistics.

    Note: This ignores z values and category.
    """
    findClusters(h5_name, density_factor, min_size)
    dbscanAnalysis.clusterStats(h5_name, min_size)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Voronoi based clustering following Levet, Nature Methods, 2015')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the HDF5 localizations file.")
    parser.add_argument('--density', dest='density', type=float, required=True,
                        help = "The density multiplier to be in a cluster. The median polygon size is multiplied by this value to give a threshold for polygon size to be in a cluster.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=30,
                        help = "The minimum cluster size to include when calculating cluster statistics.")

    args = parser.parse_args()

    voronoiAnalysis(args.hdf5, args.density, args.min_size)


