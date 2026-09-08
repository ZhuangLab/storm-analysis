#!/usr/bin/env python
"""
Classes for simulating different dye photophysics.

Hazen 11/16
"""

import numpy
import random

import storm_analysis.sa_library.sa_h5py as saH5Py

import storm_analysis.simulator.simbase as simbase


class PhotoPhysics(simbase.SimBase):
    """
    Returns location and intensity (peak height in photons) of
    the emitters that are on in the current frame.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data):
        super(PhotoPhysics, self).__init__(sim_fp, x_size, y_size, h5_data)

        
class AlwaysOn(PhotoPhysics):
    """
    All the emitters are on all the time.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, photons = 2000):
        super(AlwaysOn, self).__init__(sim_fp, x_size, y_size, h5_data)
        self.saveJSON({"photophysics" : {"class" : "AlwaysOn",
                                         "photons" : str(photons)}})
        self.h5_data['sum'] = photons * numpy.ones(self.h5_data['x'].size)

    def getEmitters(self, frame):
        temp = {}
        for key in self.h5_data:
            temp[key] = self.h5_data[key].copy()
            
        return temp

    
class AlwaysOnMC(PhotoPhysics):
    """
    All the emitters are on all the time, adapted for multicolor simulations.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, color = 0, photons = 2000):
        super(AlwaysOnMC, self).__init__(sim_fp, x_size, y_size, h5_data)
        self.saveJSON({"photophysics" : {"class" : "AlwaysOnMC",
                                         "color" : str(color),
                                         "photons" : str(photons)}})
        mask = (self.h5_data['color'] != color)
        self.h5_data['sum'] = photons * numpy.ones(self.h5_data['x'].size)
        self.h5_data['sum'][mask] = 0.1 * self.h5_data['sum'][mask]

    def getEmitters(self, frame):
        temp = {}
        for key in self.h5_data:
            temp[key] = self.h5_data[key].copy()
            
        return temp


class Duplicate(PhotoPhysics):
    """
    This just returns the localizations in an existing HDF5 file. It is perhaps
    most useful for simulating multiplane data.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, h5_name, cx = None, cy = None, z_offset = None):
        super(Duplicate, self).__init__(sim_fp, x_size, y_size, h5_data)
        self.saveJSON({"photophysics" : {"class" : "Duplicate",
                                         "h5_name" : h5_name,
                                         "cx" : str(cx),
                                         "cy" : str(cy),
                                         "z_offset" : str(z_offset)}})

        self.h5_data = saH5Py.SAH5Py(h5_name)
        self.cx = cx
        self.cy = cy
        self.z_offset = z_offset

    def getEmitters(self, frame):
        locs = self.h5_data.getLocalizationsInFrame(frame)

        temp = {}
        
        # Check for no localization data for this frame. If there are none add the minimum
        # needed fields.
        if not locs:
            for elt in ['x', 'y', 'z', 'sum']:
                temp[elt] = numpy.array([])
                    
        else:
            for elt in locs:
                temp[elt] = locs[elt]

            # Apply transforms if requested.
            if self.cx is not None:
                xf = self.cx[0] + self.cx[1] * temp['x'] + self.cx[2] * temp['y']
                yf = self.cy[0] + self.cy[1] * temp['x'] + self.cy[2] * temp['y']
                temp['x'] = xf
                temp['y'] = yf

            if self.z_offset is not None:
                temp['z'] += self.z_offset
                
        return temp
    
    
class STORM(PhotoPhysics):
    """
    Emitters go on and off with an exponentially distributed waiting times.
    Emitters can have a different brightness when on.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data):
        super(STORM, self).__init__(sim_fp, x_size, y_size, h5_data)
        self.n_emitters = self.h5_data['x'].size
        
        # Sub-classes need to set these as appropriate.
        self.photons = None
        self.off_time = None
        self.on_time = None

    def getEmitters(self, frame):
        integrated_on = numpy.zeros(self.n_emitters)

        #
        # This is a little complicated because we are trying to include accurate
        # modeling of emitters that turned on/off more than once in a single frame.
        #
        for i in range(self.n_emitters):
                
            # The easy case, no change in state in the current frame.
            if (self.next_transistion[i] >= (frame + 1.0)):
                if self.am_on[i]:
                    integrated_on[i] = 1.0

            else:
                last_transistion = frame
                while (self.next_transistion[i] < (frame + 1.0)):
                    if self.am_on[i]:
                        integrated_on[i] += self.next_transistion[i] - last_transistion

                    self.am_on[i] = not self.am_on[i]
                    last_transistion = self.next_transistion[i]
                    if self.am_on[i]:
                        self.next_transistion[i] += numpy.random.exponential(self.on_time[i])
                    else:
                        self.next_transistion[i] += numpy.random.exponential(self.off_time[i])
                    
                # Turned on and not off again in the current frame.
                if self.am_on[i]:
                    integrated_on[i] += (frame + 1.0) - last_transistion

        # Set sum and return only those emitters that have sum > 0.
        self.h5_data['sum'] = integrated_on * self.photons
        mask = (self.h5_data['sum'] > 0.0)

        temp = {}
        for key in self.h5_data:
            temp[key] = self.h5_data[key][mask].copy()

        return temp    


class ComplexSTORM(STORM):
    """
    Each emitter has different on & off kinetics.
    """
    def __init__(self, sim_fp, x_size, y_size, h5_data):
        super(ComplexSTORM, self).__init__(sim_fp, x_size, y_size, h5_data)    

        # The on/off time and intensity come from the H5 file. We'll crash
        # here if these fields are not available.
        self.photons = self.h5_data['photons']
        self.off_time = self.h5_data['off_time']
        self.on_time = self.h5_data['on_time']

        self.n_emitters = self.h5_data['x'].size
        self.saveJSON({"photophysics" : {"class" : "ComplexSTORM",
                                         "photons" : str(numpy.mean(self.photons)),
                                         "on_time" : str(numpy.mean(self.on_time)),
                                         "off_time" : str(numpy.mean(self.off_time))}})

        # Initially all the emitters are off.
        self.am_on = numpy.zeros(self.n_emitters, dtype = numpy.bool_)
        self.next_transistion = numpy.random.exponential(self.off_time, self.n_emitters)
        
    
class SimpleSTORM(STORM):
    """
    Each emitter on for 1 frame out of every 1000 frames 
    on average, both are exponentially distributed.

    Args:
        on_time : Average on time in frames.
        off_time : Average off time in frames.

    """
    def __init__(self, sim_fp, x_size, y_size, h5_data, photons = 2000, on_time = 1.0, off_time = 1000.0):
        super(SimpleSTORM, self).__init__(sim_fp, x_size, y_size, h5_data)
        
        self.photons = numpy.ones(self.n_emitters)*photons
        self.off_time = numpy.ones(self.n_emitters)*off_time
        self.on_time = numpy.ones(self.n_emitters)*on_time

        self.n_emitters = self.h5_data['x'].size
        self.saveJSON({"photophysics" : {"class" : "SimpleSTORM",
                                         "photons" : str(photons),
                                         "on_time" : str(on_time),
                                         "off_time" : str(off_time)}})

        # Initially all the emitters are off.
        self.am_on = numpy.zeros(self.n_emitters, dtype = numpy.bool_)
        self.next_transistion = numpy.random.exponential(self.off_time, self.n_emitters)


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
