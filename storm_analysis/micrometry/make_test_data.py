#!/usr/bin/env python
"""
Generate test data sets for micrometry. This will create
a files called 'locs1.bin' and 'locs2.bin' in the current
directory.

Hazen 07/17
"""
import argparse
import numpy

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.sa_library.writeinsight3 as writeinsight3

parser = argparse.ArgumentParser(description = 'Generate test data for Micrometry')

parser.add_argument('--total', dest='total', type=int, required=True,
                    help = "The total number of localizations.")
parser.add_argument('--match', dest='match', type=float, required=True,
                    help = "The fraction of localizations that match.")

args = parser.parse_args()

im_size = 512

# Create matching points.
total_match = int(args.total * args.match)
m_data = i3dtype.createDefaultI3Data(total_match)

i3dtype.posSet(m_data, "x", numpy.random.uniform(high = im_size, size = total_match))
i3dtype.posSet(m_data, "y", numpy.random.uniform(high = im_size, size = total_match))

# Create noise 1.
total_noise = args.total - total_match
n1_data = i3dtype.createDefaultI3Data(total_noise)

i3dtype.posSet(n1_data, "x", numpy.random.uniform(high = im_size, size = total_noise))
i3dtype.posSet(n1_data, "y", numpy.random.uniform(high = im_size, size = total_noise))

# Create noise 2.
n2_data = i3dtype.createDefaultI3Data(total_noise)

i3dtype.posSet(n2_data, "x", numpy.random.uniform(high = im_size, size = total_noise))
i3dtype.posSet(n2_data, "y", numpy.random.uniform(high = im_size, size = total_noise))

# Save data sets.
with writeinsight3.I3Writer("locs1.bin") as i3w:
    i3w.addMolecules(m_data)
    i3w.addMolecules(n1_data)
    
with writeinsight3.I3Writer("locs2.bin") as i3w:
    i3w.addMolecules(m_data)
    i3w.addMolecules(n2_data)

