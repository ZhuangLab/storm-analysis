#!/usr/bin/env python
"""
Collate FRC results.

Hazen 01/18
"""
import numpy

import storm_analysis.diagnostics.frc.settings as settings


def collate():

    for i, reps in enumerate(settings.n_reps):
        data = numpy.loadtxt("test_{0:02d}/frc.txt".format(i+1))

        # Find point where FRC first drops below 0.17.
        for j in range(data.shape[0]):
            if (data[j,1] < 0.17):
                print("{0:d} reps, resolution is {1:0.1f} nm".format(reps, 1.0/data[j,0]))
                break


if (__name__ == "__main__"):
    collate()
    
        
        
