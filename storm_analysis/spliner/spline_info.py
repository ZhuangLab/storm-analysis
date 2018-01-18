#!/usr/bin/env python
"""
Prints out the details of a spline.

Hazen 01/18
"""
import argparse
import pickle


parser = argparse.ArgumentParser(description = 'Print spline information.')

parser.add_argument('--spline', dest='spline', type=str, required=True,
                    help = "The name of the spline file.")

args = parser.parse_args()

with open(args.spline, "rb") as fp:
    data = pickle.load(fp)

print("Size is", data["spline"].shape[0]/2, "pixels.")
print("Z range is", data["zmin"], "to", data["zmax"], "nanometers.")
