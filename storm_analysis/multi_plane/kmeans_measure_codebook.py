#!/usr/bin/env python
"""
Measure vectors for kmeans classification. Typically this would
be run on data that comes from samples of sparse single color
dyes. For example, one slide would be a Alexa-647 and another
would be Biotium CF680.

Multiple HDF5 files can be merged into a single file with
sa_utilities.merge_hdf5.py

Hazen 01/18
"""
import numpy
import scipy
import scipy.cluster

import storm_analysis.sa_library.sa_h5py as saH5Py

def KMeansMeasureVectors(h5_filename, n_categories):

    with saH5Py.SAH5Py(h5_filename) as h5:
        assert h5.hasTracks(), "No tracking information."        
        
        # Load heights.
        #
        # Note: This loads all the tracks so this file is hopefully not
        #       too large.
        height_names = ["height"]
        for i in range(1,h5.getNChannels()):
            height_names.append(h5.getChannelPrefix(i) + "height")

        tracks = h5.getTracks(fields = height_names)

        # Create features array.
        features = None
        for i, elt in enumerate(height_names):

            # Initialize features array if necessary.
            if features is None:
                features = numpy.zeros((tracks[elt].size, h5.getNChannels()))
            features[:,i] = tracks[elt]
    
    # Normalize by total height.
    total = numpy.sum(features, axis = 1)
    for i in range(features.shape[0]):
        features[i,:] = features[i,:]/total[i]
    
    # Whiten the features as recommended by Scipy.
    features = scipy.cluster.vq.whiten(features)

    # K-Means clustering.
    [codebook, distortion] = scipy.cluster.vq.kmeans(features, n_categories)

    return codebook


if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Determine codebook for k-means classification.')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "Localization HDF5 file.")
    parser.add_argument('--ndyes', dest='ndyes', type=int, required=True,
                        help = "The number of different types of dyes to (attempt) to resolve.")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the file for the k-means codebook.")

    args = parser.parse_args()

    codebook = KMeansMeasureVectors(args.hdf5, args.ndyes)

    print("Codebook:")
    print(codebook)
    
    numpy.save(args.output, codebook)
    
