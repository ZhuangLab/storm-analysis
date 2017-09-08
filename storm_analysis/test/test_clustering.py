#!/usr/bin/env python

import storm_analysis


def test_dbscan_clustering():
    
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
        n_clusters = fp.readlines()
    if (n_clusters != 99):
        raise Exception("DBSCAN did not identify the expected number of clusters.")

    # Make pictures.
    clist_name = storm_analysis.getPathOutputTest("test_clustering_aclusters_size_list.bin")
    image_name = storm_analysis.getPathOutputTest("test_clustering_db")

    from storm_analysis.dbscan.cluster_images import clusterImages
    clusterImages(clist_name, "DBSCAN Clustering", 50, 20, image_name, [256, 256])


def test_voronoi_clustering():
    
    # Test voronoi
    alist_name = storm_analysis.getData("test/data/test_clustering_list.bin")
    output_dir = storm_analysis.getPathOutputTest("./")

    from storm_analysis.voronoi.voronoi_analysis import voronoiAnalysis
    voronoiAnalysis(alist_name, 0.1, output_dir)

    # Verify number of clusters found.
    stats_file = storm_analysis.getPathOutputTest("test_clustering_srt_stats.txt")
    with open(stats_file) as fp:
        n_clusters = fp.readlines()
    if (n_clusters != 100):
        raise Exception("Voronoi did not identify the expected number of clusters.")
    
    # Make pictures.
    clist_name = storm_analysis.getPathOutputTest("test_clustering_srt_size_list.bin")
    image_name = storm_analysis.getPathOutputTest("test_clustering_vr")

    from storm_analysis.dbscan.cluster_images import clusterImages
    clusterImages(clist_name, "Voronoi Clustering", 50, 20, image_name, [256, 256])
    

if (__name__ == "__main__"):
    test_dbscan_clustering()
    test_voronoi_clustering()
    
