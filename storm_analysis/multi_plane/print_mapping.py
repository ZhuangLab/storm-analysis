#!/usr/bin/env python
"""
Print out a mapping file. Primarily a debugging tool.

Hazen 05/17
"""

import pickle

def printMapping(filename):
    with open(filename, 'rb') as fp:
        mappings = pickle.load(fp)
    
    for key in sorted(mappings):
        print(key, mappings[key])


if (__name__ == "__main__"):
    import argparse
    parser = argparse.ArgumentParser(description = 'Print out a mapping file.')

    parser.add_argument('--mapping', dest='mapping', type=str, required=True,
                        help = "The name of the mapping file.")

    args = parser.parse_args()

    printMapping(args.mapping)
