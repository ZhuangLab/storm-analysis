#!/usr/bin/env python
"""
Print out a mapping file. A debugging tool.

Hazen 05/17
"""

import argparse
import pickle

parser = argparse.ArgumentParser(description = 'Print out a mapping file.')

parser.add_argument('--mapping', dest='mapping', type=str, required=True,
                    help = "The name of the mapping file.")

args = parser.parse_args()

with open(args.mapping, 'rb') as fp:
    mappings = pickle.load(fp)
    
for key in sorted(mappings):
    print(key, mappings[key])
