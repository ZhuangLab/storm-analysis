#!/usr/bin/env python

import storm_analysis


def test_track_average_correct():

    mlist_name = storm_analysis.getData("test/data/test_drift_mlist.bin")
    settings = storm_analysis.getData("test/data/test_drift.xml")
    alist_name = storm_analysis.getPathOutputTest("test_drift_alist.bin")

    from storm_analysis.sa_utilities.track_average_correct import trackAverageCorrect
    trackAverageCorrect(mlist_name, alist_name, settings)


def test_dbscan():

    alist_name = storm_analysis.getPathOutputTest("test_drift_alist.bin")

    from storm_analysis.dbscan.dbscan_analysis import dbscanAnalysis
    dbscanAnalysis(alist_name, 0)

    clist_name = storm_analysis.getPathOutputTest("test_drift_aclusters_size_list.bin")
    image_name = storm_analysis.getPathOutputTest("test_drift_db")

    from storm_analysis.dbscan.cluster_images import clusterImages
    clusterImages(clist_name, "DBSCAN Clustering", 50, 20, image_name)


def test_voronoi():

    alist_name = storm_analysis.getPathOutputTest("test_drift_alist.bin")
    output_dir = storm_analysis.getPathOutputTest("./")

    from storm_analysis.voronoi.voronoi_analysis import voronoiAnalysis
    voronoiAnalysis(alist_name, 1.25, output_dir)

    clist_name = storm_analysis.getPathOutputTest("test_drift_asrt_size_list.bin")
    image_name = storm_analysis.getPathOutputTest("test_drift_vr")

    from storm_analysis.dbscan.cluster_images import clusterImages
    clusterImages(clist_name, "Voronoi Clustering", 50, 20, image_name)    
    

if (__name__ == "__main__"):
    test_track_average_correct()
    test_dbscan()
    test_voronoi()
