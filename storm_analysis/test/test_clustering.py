#!/usr/bin/env python
"""
Tests for DBSCAN and Voronoi clustering.
"""
import numpy

import storm_analysis
import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.dbscan.clusters_sa_h5py as clSAH5Py
import storm_analysis.dbscan.dbscan_analysis as dbscanAnalysis
import storm_analysis.voronoi.voronoi_analysis as voronoiAnalysis


def test_dbscan_clustering_1():
    numpy.random.seed(1)
    
    filename = "test_clustering_sa_h5py.hdf5"
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
        assert(cl_h5.getNClusters() == 10)
        for index, cluster in cl_h5.clustersIterator(skip_unclustered = False, fields = ["x", "y", "z"]):
            for elt in ['x', 'y', 'z']:
                assert(cluster[elt].size == 100)
                dev = numpy.std(cluster[elt])
                assert(dev > 0.085)
                assert(dev < 0.115)

    # Calculate common cluster statistics.
    stats_name = dbscanAnalysis.clusterStats(h5_name, 100)

    # Check statistics.
    stats = numpy.loadtxt(stats_name, skiprows = 1)
    assert(stats.shape[0] == 10)
    assert(numpy.allclose(stats[:,0], numpy.arange(10) + 1))
    assert(numpy.allclose(stats[:,1], numpy.zeros(10)))
    assert(numpy.allclose(stats[:,2], numpy.zeros(10) + 100.0))
    assert(numpy.allclose(stats[:,3], x, rtol = 0.1, atol = 1.0))
    assert(numpy.allclose(stats[:,4], y, rtol = 0.1, atol = 1.0))
    assert(numpy.allclose(stats[:,5], z, rtol = 0.1, atol = 20.0))


def test_voronoi_clustering_1():
    numpy.random.seed(1)
    
    filename = "test_clustering_sa_h5py.hdf5"
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

    # Cluster data with voronoi.
    voronoiAnalysis.findClusters(h5_name, 0.1, 10, verbose = False)

    # Check clustering results.
    with clSAH5Py.SAH5Clusters(h5_name) as cl_h5:
        assert(cl_h5.getNClusters() == 10)
        for index, cluster in cl_h5.clustersIterator(skip_unclustered = True, fields = ["x", "y", "z"]):
            for elt in ['x', 'y', 'z']:
                dev = numpy.std(cluster[elt])
                assert(dev > 0.07)
                assert(dev < 0.12)

    # Calculate common cluster statistics.
    stats_name = dbscanAnalysis.clusterStats(h5_name, 50, verbose = False)

    # Check statistics.
    stats = numpy.loadtxt(stats_name, skiprows = 1)
    index = numpy.argsort(stats[:,3])
    assert(stats.shape[0] == 10)
    assert(numpy.allclose(stats[:,0], numpy.arange(10) + 1))
    assert(numpy.allclose(stats[:,1], numpy.zeros(10)))
    assert(numpy.count_nonzero(numpy.greater(stats[:,2], 80.0 * numpy.ones(10))) == 10)
    assert(numpy.allclose(stats[index,3], x, rtol = 0.2, atol = 2.0))
    assert(numpy.allclose(stats[index,4], y, rtol = 0.2, atol = 2.0))
    assert(numpy.allclose(stats[index,5], z, rtol = 0.2, atol = 20.0))

    
if (__name__ == "__main__"):
    #test_dbscan_clustering_1()
    test_voronoi_clustering_1()
    
