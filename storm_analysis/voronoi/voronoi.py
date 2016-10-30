#!/usr/bin/env python
#
# A Python implementation of some of the ideas in the SR-Tesseler paper.
# Basically this does is calculate the area (in pixels) of the Voroni
# region around a localization and stores that in the localizations fit
# area field.
#
# Note: This ignores the localization category.
#
# Note: This will handle up to on the order of 1M localizations. Analysis
#       of files with a lot more localizations than this will likely
#       take a long time to analyze.
#
# Hazen 09/16
#

import numpy
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import Polygon

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


def voronoi(mlist_name, clist_name, density_factor, min_size, verbose = True):

    i3_data_in = readinsight3.loadI3GoodOnly(mlist_name)
    n_locs = i3_data_in['xc'].size
    points = numpy.column_stack((i3_data_in['xc'], i3_data_in['yc']))

    print("Creating Voronoi object.")
    vor = Voronoi(points)

    print("Calculating 2D region sizes.")
    for i, region_index in enumerate(vor.point_region):
        if ((i%10000) == 0):
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
            i3_data_in['a'][i] = 1.0/area

    # Used median density based threshold.
    ave_density = numpy.median(i3_data_in['a'])
    if verbose:
        print("Min density", numpy.min(i3_data_in['a']))
        print("Max density", numpy.max(i3_data_in['a']))
        print("Median density", ave_density)

    # Record the neighbors of each point.
    max_neighbors = 40
    neighbors = numpy.zeros((n_locs, max_neighbors), dtype = numpy.int32) - 1
    neighbors_counts = numpy.zeros((n_locs), dtype = numpy.int32)

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
    print("Marking connected regions")
    i3_data_in['lk'] = -1
    min_density = density_factor * ave_density
    visited = numpy.zeros((n_locs), dtype = numpy.int32)

    def neighborsList(index):
        nlist = []
        for i in range(neighbors_counts[index]):
            loc_index = neighbors[index,i]
            if (visited[loc_index] == 0):
                nlist.append(neighbors[index,i])
                visited[loc_index] = 1
        return nlist

    cluster_id = 2
    for i in range(n_locs):
        if (visited[i] == 0):
            if (i3_data_in['a'][i] > min_density):
                cluster_elt = [i]
                c_size = 1
                to_check = neighborsList(i)
                while (len(to_check) > 0):

                    # Remove last localization from the list.
                    loc_index = to_check[-1]
                    to_check = to_check[:-1]

                    # If the localization has sufficient density add to cluster and check neighbors.
                    if (i3_data_in['a'][loc_index] > min_density):
                        to_check += neighborsList(loc_index)
                        cluster_elt.append(loc_index)
                        c_size += 1

                    # Mark as visited.
                    visited[loc_index] = 1

                # Mark the cluster if there are enough localizations in the cluster.
                if (c_size > min_size):
                    print(cluster_id, c_size)
                    for elt in cluster_elt:
                        i3_data_in['lk'][elt] = cluster_id
                cluster_id += 1
            visited[i] = 1

    print(cluster_id, "clusters")
    
    # Save the data.
    print("Saving results")
    i3_data_out = writeinsight3.I3Writer(clist_name)
    i3_data_out.addMolecules(i3_data_in)
    i3_data_out.close()
