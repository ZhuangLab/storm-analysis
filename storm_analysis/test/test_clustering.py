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
        

def _test_dbscan_clustering():
    
    # Test dbscan
    import shutil

    # Copy alist to the output directory so that the DBSCAN results end up in the right place.
    alist_data = storm_analysis.getData("test/data/test_clustering_list.bin")
    alist_output = storm_analysis.getPathOutputTest("test_clustering_alist.bin")
    shutil.copyfile(alist_data, alist_output)
    
    from storm_analysis.dbscan.dbscan_analysis import dbscanAnalysis
    dbscanAnalysis(alist_output, 0)

    # Verify number of clusters found.
    stats_file = storm_analysis.getPathOutputTest("test_clustering_aclusters_stats.txt")
    with open(stats_file) as fp:
        n_clusters = len(fp.readlines())
    if (n_clusters != 99):
        raise Exception("DBSCAN did not identify the expected number of clusters.")

    # Make pictures.
    clist_name = storm_analysis.getPathOutputTest("test_clustering_aclusters_size_list.bin")
    image_name = storm_analysis.getPathOutputTest("test_clustering_db")

    from storm_analysis.dbscan.cluster_images import clusterImages
    clusterImages(clist_name, "DBSCAN Clustering", 50, 20, image_name, [256, 256])


def _test_voronoi_clustering():
    
    # Test voronoi
    alist_name = storm_analysis.getData("test/data/test_clustering_list.bin")
    output_dir = storm_analysis.getPathOutputTest("./")

    from storm_analysis.voronoi.voronoi_analysis import voronoiAnalysis
    voronoiAnalysis(alist_name, 0.1, output_dir)

    # Verify number of clusters found.
    stats_file = storm_analysis.getPathOutputTest("test_clustering_srt_stats.txt")
    with open(stats_file) as fp:
        n_clusters = len(fp.readlines())
    if (n_clusters != 100):
        raise Exception("Voronoi did not identify the expected number of clusters.")
    
    # Make pictures.
    clist_name = storm_analysis.getPathOutputTest("test_clustering_srt_size_list.bin")
    image_name = storm_analysis.getPathOutputTest("test_clustering_vr")

    from storm_analysis.dbscan.cluster_images import clusterImages
    clusterImages(clist_name, "Voronoi Clustering", 50, 20, image_name, [256, 256])
    

if (__name__ == "__main__"):
    test_dbscan_clustering_1()

    
