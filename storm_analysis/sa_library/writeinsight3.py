#!/usr/local/bin/python
#
# Writes Insight3 format binary molecule lists.
#
# Hazen 6/09
#

import numpy
import struct

import i3dtype

def _putV(fp, format, data):
    fp.write(struct.pack(format, data))

class I3Writer():

    def __init__(self, filename, frames = 1):
        self.molecules = 0
        self.fp = open(filename, "wb")
        _putV(self.fp, "4s", "M425")
        _putV(self.fp, "i", frames)
        _putV(self.fp, "i", 6)
        _putV(self.fp, "i", 0)

    def __enter__(self):
        return self

    def __exit__(self, etype, value, traceback):
        if self.fp:
            self.close()
            
    #
    # This is for localizations identified by the original DAOSTORM
    # algorithm, not the 3D-DAOSTORM algorithm.
    #
    def addDAOSTORMMolecules(self, frame, xc, yc, br, be, msky, niter, sharp, chi, err):
        #
        # DAOSTORM -> Insight3 format mapping.
        #
        # xc - xcenter
        # yc - ycenter
        # br - brightness -> peak height
        # be - brightness error (?) -> peak area
        # msky - background -> peak background
        # niter - fit iterations
        # sharp - sharpness (?) -> peak angle
        # chi - fit quality -> peak width
        # err - error flag -> link
        #

        i3data = i3dtype.createDefaultI3Data(xc.size)
        i3dtype.posSet(i3data, 'x', xc)
        i3dtype.posSet(i3data, 'y', yc)
        i3dtype.setI3Field(i3data, 'h', br)
        i3dtype.setI3Field(i3data, 'a', be)
        i3dtype.setI3Field(i3data, 'bg', msky)
        i3dtype.setI3Field(i3data, 'fi', niter)
        i3dtype.setI3Field(i3data, 'phi', sharp)
        i3dtype.setI3Field(i3data, 'w', chi)
        i3dtype.setI3Field(i3data, 'lk', err)

        self.addMolecules(i3data)
        
    def addMolecules(self, i3data):
        i3data.tofile(self.fp)
        self.molecules += i3data['x'].size
        #self.fp.flush()

    # Various Convenience functions
    def addMoleculesWithXY(self, x, y):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        self.addMolecules(i3data)

    def addMoleculesWithXYAFrame(self, x, y, pa, frame):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'a', pa)
        i3dtype.setI3Field(i3data, 'fr', frame)
        self.addMolecules(i3data)

    def addMoleculesWithXYAItersFrame(self, x, y, pa, iters, frame):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'a', pa)
        i3dtype.setI3Field(i3data, 'fi', iters)
        i3dtype.setI3Field(i3data, 'fr', frame)
        self.addMolecules(i3data)

    def addMoleculesWithXYCat(self, x, y, cat):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'c', cat)
        self.addMolecules(i3data)

    def addMoleculesWithXYCatF(self, x, y, cat,f):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'c', cat)
        i3dtype.setI3Field(i3data, 'fr', f)
        self.addMolecules(i3data)

    def addMoleculesWithXYF(self, x, y, f):
        self.addMoleculesWithXYFrame(x, y, f)

    def addMoleculesWithXYFrame(self, x, y, frame):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'fr', frame)
        self.addMolecules(i3data)

    def addMoleculesWithXYI(self, x, y, pi):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'i', pi)
        self.addMolecules(i3data)

    def addMoleculesWithXYICat(self, x, y, pi, cat):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'i', pi)
        i3dtype.setI3Field(i3data, 'c', cat)
        self.addMolecules(i3data)

    def addMoleculesWithXYIFrame(self, x, y, pi, frame):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'i', pi)
        i3dtype.setI3Field(i3data, 'fr', frame)
        self.addMolecules(i3data)

    def addMoleculesWithXYIWFrame(self, x, y, pi, width, frame):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.setI3Field(i3data, 'i', pi)
        i3dtype.setI3Field(i3data, 'w', width)
        i3dtype.setI3Field(i3data, 'fr', frame)
        self.addMolecules(i3data)

    def addMoleculesWithXYZ(self, x, y, z):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.posSet(i3data, 'z', z)
        self.addMolecules(i3data)

    def addMoleculesWithXYZF(self, x, y, z, f):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.posSet(i3data, 'z', z)
        i3dtype.setI3Field(i3data, 'fr', f)
        self.addMolecules(i3data)

    def addMoleculesWithXYZI(self, x, y, z, pi):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.posSet(i3data, 'z', z)
        i3dtype.setI3Field(i3data, 'i', pi)
        self.addMolecules(i3data)

    def addMoleculesWithXYZIFrame(self, x, y, z, pi, f):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.posSet(i3data, 'z', z)
        i3dtype.setI3Field(i3data, 'i', pi)
        i3dtype.setI3Field(i3data, 'fr', f)
        self.addMolecules(i3data)

    def addMoleculesWithXYZCat(self, x, y, z, cat):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.posSet(i3data, 'z', z)
        i3dtype.setI3Field(i3data, 'c', cat)
        self.addMolecules(i3data)

    def addMoleculesWithXYZICat(self, x, y, z, pi, cat):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.posSet(i3data, 'z', z)
        i3dtype.setI3Field(i3data, 'i', pi)
        i3dtype.setI3Field(i3data, 'c', cat)
        self.addMolecules(i3data)

    def addMoleculesWithXYZICatFrame(self, x, y, z, pi, cat, f):
        i3data = i3dtype.createDefaultI3Data(x.size)
        i3dtype.posSet(i3data, 'x', x)
        i3dtype.posSet(i3data, 'y', y)
        i3dtype.posSet(i3data, 'z', z)
        i3dtype.setI3Field(i3data, 'i', pi)
        i3dtype.setI3Field(i3data, 'c', cat)
        i3dtype.setI3Field(i3data, 'fr', f)
        self.addMolecules(i3data)

    #
    # This is for localization identified by 3D-DAOSTORM.
    #
    def addMultiFitMolecules(self, molecules, x_size, y_size, frame, nm_per_pixel, inverted=False):
        n_molecules = molecules.shape[0]
        
        h = molecules[:,0]
        if inverted:
            xc = y_size - molecules[:,1]
            yc = x_size - molecules[:,3]
            wx = 2.0*molecules[:,2]*nm_per_pixel
            wy = 2.0*molecules[:,4]*nm_per_pixel
        else:
            xc = molecules[:,3] + 1
            yc = molecules[:,1] + 1
            wx = 2.0*molecules[:,4]*nm_per_pixel
            wy = 2.0*molecules[:,2]*nm_per_pixel

        bg = molecules[:,5]
        zc = molecules[:,6] * 1000.0  # fitting is done in um, insight works in nm
        st = numpy.round(molecules[:,7])
        err = molecules[:,8]
        
        # calculate peak area, which is saved in the "a" field.
        parea = 2.0*3.14159*h*molecules[:,2]*molecules[:,4]

        ax = wy/wx
        ww = numpy.sqrt(wx*wy)
        
        i3data = i3dtype.createDefaultI3Data(xc.size)
        i3dtype.posSet(i3data, 'x', xc)
        i3dtype.posSet(i3data, 'y', yc)
        i3dtype.posSet(i3data, 'z', zc)
        i3dtype.setI3Field(i3data, 'h', h)
        i3dtype.setI3Field(i3data, 'bg', bg)
        i3dtype.setI3Field(i3data, 'fi', st)
        i3dtype.setI3Field(i3data, 'a', parea)
        i3dtype.setI3Field(i3data, 'w', ww)
        i3dtype.setI3Field(i3data, 'ax', ax)
        i3dtype.setI3Field(i3data, 'fr', frame)
        i3dtype.setI3Field(i3data, 'i', err)
        self.addMolecules(i3data)

    def close(self):
        print "Added", self.molecules
        _putV(self.fp, "i", 0)
        self.fp.seek(12)
        _putV(self.fp, "i", self.molecules)
        self.fp.close()
    

#
# The MIT License
#
# Copyright (c) 2012 Zhuang Lab, Harvard University
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
