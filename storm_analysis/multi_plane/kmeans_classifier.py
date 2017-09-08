#!/usr/bin/env python
"""
Given a k-means codebook and a set of paired localization files 
create a new localization file with the k-mean category in the
'c' field and distance from the nearest cluster center in
the 'i' field.

Hazen 09/17
"""
import numpy
import scipy
import scipy.cluster

import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3

def KMeansClassifier(codebook, input_basename, output_name, extensions = [".bin", "_ch1.bin", "_ch2.bin", "_ch3.bin"], max_distance = 80):
    """
    Note: 
      1. The default is that there are 4 color channels / cameras.
      2. The maximum distance is in percent, so '80' means that the 20%
         of the localizations that most distant from a cluster center
         will put in category 9.
    """
    n_channels = codebook.shape[1]
    assert (n_channels == len(extensions)), "Codebook size does not match data."

    # Create a reader for each file.
    i3_readers = []
    for ext in extensions:
        i3_name = input_basename + ext
        print(i3_name)
        i3_readers.append(readinsight3.I3Reader(i3_name))

    # Create writer for the results.
    i3_out = writeinsight3.I3Writer(output_name)

    # Read first block of the first channel data.
    i3_data = [i3_readers[0].nextBlock()]
    while (i3_data[0] is not False):
        print("working..")

        # Read the data from the other channels.
        for i in range(1,len(i3_readers)):
            i3_data.append(i3_readers[i].nextBlock())

        # Load height data for each channel.
        features = numpy.zeros((i3_data[0].size, n_channels))
        for i in range(len(i3_readers)):
            features[:,i] = i3_data[i]['h']

        # Normalize by total height.
        total = numpy.sum(features, axis = 1)
        for i in range(features.shape[0]):
            features[i,:] = features[i,:]/total[i]
    
        # Whiten the features as recommended by Scipy.
        features = scipy.cluster.vq.whiten(features)

        # Classify using codebook.
        [category, distance] = scipy.cluster.vq.vq(features, codebook)
        dist_max = numpy.percentile(distance, max_distance)

        # Put top XX% in distance in category 9 (the discard category).
        mask = (distance > dist_max)
        category[mask] = 9
            
        #
        # Store category and distance in the 'c' and 'i' field respectively.
        #
        i3_data[0]['c'] = category
        i3_data[0]['i'] = distance

        i3_out.addMolecules(i3_data[0])

        # Load the next block of data.
        i3_data = [i3_readers[0].nextBlock()]

    # Close output file
    i3_out.close()

    
if (__name__ == "__main__"):
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Use K-Means codebook to classify localizations in a file.')

    parser.add_argument('--codebook', dest='codebook', type=str, required=True,
                        help = "A K-Means codebook.")
    parser.add_argument('--basename', dest='basename', type=str, required=True,
                        help = "The basename for the localization files.")
    parser.add_argument('--output', dest='output', type=str, required=True,
                        help = "The name of the file for the categorized localizations.")
    parser.add_argument('--max_dist', dest='max_dist', type=float, required=False, default = 80.0,
                        help = "The maximum distance from a cluster center to keep as a percentile (default is 80%).")
    
    args = parser.parse_args()

    codebook = numpy.load(args.codebook)
    KMeansClassifier(codebook, args.basename, args.output, max_distance = args.max_dist)
