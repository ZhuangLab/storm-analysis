#!/usr/bin/env python
"""
Measure vectors for kmeans classification. Typically this would
be run on data that comes from samples of sparse single color
dyes. For example, one slide would be a Alexa-647 and another
would be Biotium CF680.

Hazen 09/17
"""
import numpy
import scipy
import scipy.cluster

def KMeansMeasureVectors(height_fnames, n_categories):
    
    #
    # Load height data from each file.
    #
    features = None
    for i, ht_fname in enumerate(height_fnames):
        h_data = numpy.load(ht_fname)

        # Initialize features array if necessary.
        if features is None:
            features = h_data
        else:
            features = numpy.concatenate((features, h_data), axis = 0)

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

    parser.add_argument('--heights', dest='heights', type=str, required=True, nargs = '*',
                        help = "One or more numpy height files as created by multiplane/merge_heights.py.")
    parser.add_argument('--ndyes', dest='ndyes', type=int, required=True,
                        help = "The number of different types of dyes to (attempt) to resolve.")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the file for the k-means codebook.")

    args = parser.parse_args()

    codebook = KMeansMeasureVectors(args.heights, args.ndyes)

    print("Codebook:")
    print(codebook)
    
    numpy.save(args.output, codebook)
    
