#!/usr/bin/env python
"""
Generate synthetic localization data for clustering.

Hazen 08/18
"""

import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py


def makeClusters(h5_name, n_clusters, cluster_size, n_background):
    """
    Make the fake localization data. The data is saved as tracks. The
    movie size is 300 x 200 with 100nm pixels.

    h5_name - The name of the HDF5 file.
    n_clusters - The number of clusters.
    cluster_size - The number of localizations per cluster.
    n_background - The number of background localizations.
    """

    margin = 20
    x_size = 300
    y_size = 200
    
    with saH5Py.SAH5Py(h5_name, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(x_size, y_size, 1, "")
        h5.setPixelSize(100.0)

        cx = numpy.random.uniform(low = margin, high = x_size - margin, size = n_clusters)
        cy = numpy.random.uniform(low = margin, high = y_size - margin, size = n_clusters)

        # Background tracks are uniformly distributed across tracks.
        bg_per_grp = int(round(n_background/cluster_size))

        # There are as many track groups as elements in a cluster.
        for i in range(cluster_size):
            x = cx + numpy.random.normal(scale = 1.0, size = n_clusters)
            y = cy + numpy.random.normal(scale = 1.0, size = n_clusters)

            # Add background
            x = numpy.concatenate((x, numpy.random.uniform(high = x_size, size = bg_per_grp)))
            y = numpy.concatenate((y, numpy.random.uniform(high = y_size, size = bg_per_grp)))

            tracks = {"category" : numpy.ones(x.size),
                      "x" : x,
                      "y" : y,
                      "z" : numpy.zeros(x.size)}

            h5.addTracks(tracks)
