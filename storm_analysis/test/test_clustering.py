#!/usr/bin/env python
"""
Tests for DBSCAN and Voronoi clustering.
"""
import numpy

import storm_analysis
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.dbscan.clusters_sa_h5py as clSAH5Py
import storm_analysis.dbscan.dbscan_analysis as dbscanAnalysis


def test_dbscan_clustering_1():
    numpy.random.seed(1)
    
    filename = "test_clusters_sa_h5py.hdf5"
    h5_name = storm_analysis.getPathOutputTest(filename)
    storm_analysis.removeFile(h5_name)

    # Write tracks data.
    category = numpy.zeros(10, dtype = numpy.int32)
    x = 10.0 * numpy.arange(10)
    y = 10.0 * numpy.arange(10)
    z = numpy.zeros(10)

    with saH5Py.SAH5Py(h5_name, is_existing = False) as h5:
        h5.setMovieInformation(1,1,2,"")
        h5.setPixelSize(1.0)
        
        for i in range(100):
            tracks = {"category" : category,
                      "x" : x + numpy.random.normal(scale = 0.1, size = 10),
                      "y" : y + numpy.random.normal(scale = 0.1, size = 10),
                      "z" : z + numpy.random.normal(scale = 0.1, size = 10)}
            
            h5.addTracks(tracks)

    # Cluster data with DBSCAN.
    dbscanAnalysis.findClusters(h5_name, 2.0, 10)

    # Check clustering results.
    with clSAH5Py.SAH5Clusters(h5_name) as cl_h5:
        assert(cl_h5.getNClusters() == 9)
        for index, cluster in cl_h5.clustersIterator(skip_unclustered = False):
            for elt in ['x', 'y', 'z']:
                assert(cluster[elt].size == 100)
                dev = numpy.std(cluster[elt])
                assert(dev > 0.085)
                assert(dev < 0.115)
        

if (__name__ == "__main__"):
    test_dbscan_clustering_1()

    
