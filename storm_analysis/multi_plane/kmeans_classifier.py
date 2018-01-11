#!/usr/bin/env python
"""
Given a k-means codebook and a (tracked) localizations file, assign
a color and a distance to the nearest cluster center. These fields
are added to tracks and are available as 'km_color' and 'km_distance'.
Tracks that are to far away from a cluster center are assigned a color
value of 9.

Hazen 01/18
"""
import numpy
import scipy
import scipy.cluster

import storm_analysis.sa_library.sa_h5py as saH5Py

def KMeansClassifier(codebook, h5_filename, max_distance = 80):
    """
    Note: The maximum distance is in percent, so '80' means that the 20%
          of the localizations that are most distant from a cluster center
          will put in category 9.
    """
    n_channels = codebook.shape[1]
    with saH5Py.SAH5Py(h5_filename) as h5:
        assert (n_channels == h5.getNChannels()), "Codebook size does not match data."

        # Create height names array.
        height_names = ["height"]
        for i in range(1,h5.getNChannels()):
            height_names.append(h5.getChannelPrefix(i) + "height")

        # Iterate over tracks.
        for i, tracks in enumerate(h5.tracksIterator(fields = height_names)):

            # Create features array.
            features = None
            for j, elt in enumerate(height_names):

                # Initialize features array if necessary.
                if features is None:
                    features = numpy.zeros((tracks[elt].size, h5.getNChannels()))
                features[:,j] = tracks[elt]

            # Normalize by total height.
            total = numpy.sum(features, axis = 1)
            for j in range(features.shape[0]):
                features[j,:] = features[j,:]/total[j]
    
            # Whiten the features as recommended by Scipy.
            features = scipy.cluster.vq.whiten(features)

            # Classify using codebook.
            [km_color, km_distance] = scipy.cluster.vq.vq(features, codebook)
            dist_max = numpy.percentile(km_distance, max_distance)

            # Put top XX% in distance in color value 9.
            mask = (km_distance > dist_max)
            km_color[mask] = 9

            # Save results.
            h5.addTrackData(km_color, i, "km_color")
            h5.addTrackData(km_distance, i, "km_distance")
        
    
if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Use K-Means codebook to classify localizations in a file.')

    parser.add_argument('--codebook', dest='codebook', type=str, required=True,
                        help = "A K-Means codebook.")
    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "Localization HDF5 file.")
    parser.add_argument('--max_dist', dest='max_dist', type=float, required=False, default = 80.0,
                        help = "The maximum distance from a cluster center to keep as a percentile (default is 80%).")
    
    args = parser.parse_args()

    codebook = numpy.load(args.codebook)
    KMeansClassifier(codebook, args.hdf5, max_distance = args.max_dist)
