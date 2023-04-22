#!/usr/bin/env python
"""
Base class for the various simulation classes.

Hazen 12/16
"""

import json


class SimException(Exception):
    pass

class SimBase(object):
    """
    Base class for all the simulation classes.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data):
        self.h5_data = h5_data
        self.sim_fp = sim_fp
        self.x_size = x_size
        self.y_size = y_size

    def checkSize(self, x_size, y_size):
        if (self.x_size != x_size) or (self.y_size != y_size):
            raise SimException("Size mismatch (X) " + str(self.x_size) + " != " + x_size + " or (Y) " + self.y_size + " != " + y_size)

    def cleanup(self):
        pass
        
    def saveJSON(self, elt):
        self.sim_fp.write(json.dumps(elt) + "\n")

    def setSize(self, x_size, y_size):
        self.x_size = x_size
        self.y_size = y_size


#
# The MIT License
#
# Copyright (c) 2016 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
