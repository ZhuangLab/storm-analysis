#!/usr/bin/env python
"""
Print out the values in a PSF file. Primarily a debugging tool.

Hazen 02/18
"""

import pickle

def printPSF(filename):
    with open(filename, 'rb') as fp:
        mappings = pickle.load(fp)
    
    for key in sorted(mappings):
        if not (key == "psf"):
            print(key, mappings[key])
        else:
            print("psf shape", mappings["psf"].shape)


if (__name__ == "__main__"):
    import argparse
    parser = argparse.ArgumentParser(description = 'Print out the values in a PSF file.')

    parser.add_argument('--psf', dest='psf', type=str, required=True,
                        help = "The name of the PSF file.")

    args = parser.parse_args()

    printPSF(args.psf)
