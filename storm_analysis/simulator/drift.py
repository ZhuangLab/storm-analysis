#!/usr/bin/env python
"""
Classes for creating different kinds of drifts.

Hazen 06/17
"""

import numpy
import random

import storm_analysis.sa_library.i3dtype as i3dtype
import storm_analysis.simulator.simbase as simbase


class Drift(simbase.SimBase):
    """
    Apply frame dependent drift to the localizations.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data):
        super(Drift, self).__init__(sim_fp, x_size, y_size, h5_data)
        

class DriftFromFile(Drift):
    """
    Add drift from a file. X and Y are in units of pixels, Z is in
    nanometers.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, drift_file):
        super(DriftFromFile, self).__init__(sim_fp, x_size, y_size, h5_data)
        self.saveJSON({"drift" : {"class" : "DriftFromFile",
                                  "drift_file" : drift_file}})
        self.drift_data = numpy.loadtxt(drift_file)
    
    def drift(self, frame_number, h5_data):

        # Check that frame number is in range.
        if (frame_number >= self.drift_data.shape[0]):
            frame_number = self.drift_data.shape[0] - 1
            
        dx = self.drift_data[frame_number, 0]
        dy = self.drift_data[frame_number, 1]
        dz = self.drift_data[frame_number, 2]

        h5_data["x"] += dx
        h5_data["y"] += dy
        h5_data["z"] += dz

        
