#!/usr/bin/env python
"""
Generate simulated. The basic idea is that you provide a list of
localizations in Insight3 .bin format and these are used to
generate a series of images using the following steps:

Initialization:
  1. locs = h5.getLocalizations()
  2. bg = background.Background()
  3. camera = camera.Camera()
  4. drift = drift.Drift()
  5. pp = photophysics.Photophysics()
  6. psf = psf.PSF()


Generation:
  1. image = numpy.zeros()
  2. image += bg.getBackground()
  3. cur_locs = pp.getEmitters()
  4. drift.drift(cur_locs, frame_number)
  5. image += psf.getPSFs(cur_locs)
  6. image = camera.readImage(image)
  7. saveimage()
  8. savelocs()

Note: This is expected to set the 'height','sum' and 'background' fields 
    in the output to the correct values. The values 'x', 'y', 'z', 'xsigma' 
    and 'ysigma' are just passed through.

Hazen 01/18
"""

import json
import numpy

import storm_analysis.sa_library.datawriter as datawriter
import storm_analysis.sa_library.sa_h5py as saH5Py


class Simulate(object):

    def __init__(self,
                 background_factory = None,
                 camera_factory = None,
                 drift_factory = None,
                 photophysics_factory = None,
                 psf_factory = None,
                 dither = False,
                 x_size = 256,
                 y_size = 256,
                 **kwds):
        """
        The factory variables should be functions that return the correct class
        to run a simulation with the following signature:

        factory_fn(sim_settings, x_size, y_size, h5_data_in)

        """
        super(Simulate, self).__init__(**kwds)
        
        self.bg_factory = background_factory
        self.cam_factory = camera_factory
        self.drift_factory = drift_factory
        self.pphys_factory = photophysics_factory
        self.psf_factory = psf_factory

        self.dither = dither
        self.x_size = x_size
        self.y_size = y_size

    def setBackgroundFactory(self, new_factory):
        self.bg_factory = new_factory

    def setCameraFactory(self, new_factory):
        self.cam_factory = new_factory

    def setDriftFactory(self, new_drift_factory):
        self.drift_factory = new_drift_factory

    def setPhotoPhysicsFactory(self, new_factory):
        self.pphys_factory = new_factory

    def setPSFFactory(self, new_factory):
        self.psf_factory = new_factory

    def simulate(self, dax_file, bin_file, n_frames):

        #
        # Initialization.
        #
        movie_data = datawriter.inferWriter(dax_file,
                                            width = self.x_size,
                                            height = self.y_size)
        with saH5Py.SAH5Py(bin_file) as h5:
            h5_data_in = h5.getLocalizations()

        out_fname_base = dax_file[:-4]
        h5_data_out = saH5Py.SAH5Py(filename = out_fname_base + "_ref.hdf5",
                                    is_existing = False,
                                    overwrite = True)
        h5_data_out.setMovieInformation(self.x_size, self.y_size, n_frames, "")
        
        sim_settings = open(out_fname_base + "_sim_params.txt", "w")

        #
        # Create the user-specified class instances that will do
        # most of the actual work of the simulation.
        #
        bg = self.bg_factory(sim_settings, self.x_size, self.y_size, h5_data_in)
        cam = self.cam_factory(sim_settings, self.x_size, self.y_size, h5_data_in)
        drift = None
        if self.drift_factory is not None:
            drift = self.drift_factory(sim_settings, self.x_size, self.y_size, h5_data_in)
        pp = self.pphys_factory(sim_settings, self.x_size, self.y_size, h5_data_in)
        psf = self.psf_factory(sim_settings, self.x_size, self.y_size, h5_data_in)

        sim_settings.write(json.dumps({"simulation" : {"bin_file" : bin_file,
                                                       "x_size" : str(self.x_size),
                                                       "y_size" : str(self.y_size)}}) + "\n")

        #
        # Generate the simulated movie.
        #
        for i in range(n_frames):

            # Generate the new image.
            image = numpy.zeros((self.x_size, self.y_size))

            # Get the emitters that are on in the current frame.
            cur_h5 = pp.getEmitters(i)

            print("Frame", i, cur_h5['x'].size, "emitters")

            # Dither points x,y values if requested. This is useful for things
            # like looking for pixel level biases in simulated data with gridded
            # localizations.
            #
            if self.dither:
                cur_h5['x'] += numpy.random.uniform(size = cur_h5['x'].size) - 0.5
                cur_h5['y'] += numpy.random.uniform(size = cur_h5['y'].size) - 0.5

            # Add background to image.
            image += bg.getBackground(i)

            # Set 'bg' parameter of the emitters.
            cur_h5 = bg.getEmitterBackground(cur_h5)

            # Apply drift to the localizations.
            if drift is not None:
                drift.drift(i, cur_h5)
            
            # Foreground
            image += psf.getPSFs(cur_h5)

            # Camera
            image = cam.readImage(image)

            # Save the image.
            movie_data.addFrame(numpy.transpose(image))

            # Save the molecule locations.
            h5_data_out.addLocalizations(cur_h5, i)

        movie_data.close()
        h5_data_out.close()
        sim_settings.close()


# The simplest simulation.

if (__name__ == "__main__"):

    import argparse

    import storm_analysis.simulator.background as background
    import storm_analysis.simulator.camera as camera
    import storm_analysis.simulator.photophysics as photophysics
    import storm_analysis.simulator.psf as psf

    parser = argparse.ArgumentParser(description = "A simple simulation example.")

    parser.add_argument('--dax', dest='dax_file', type=str, required=True,
                        help = "The name of the dax file to save the simulated STORM movie.")
    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the HDF5 file containing the emitter locations.")
    parser.add_argument('--frames', dest='frames', type=int, required=True,
                        help = "The length of the movie in frames.")
    parser.add_argument('--photons', dest='photons', type=float, required=True,
                        help = "The integral of a single emitter in photons.")

    args = parser.parse_args()
 
    sim = Simulate(lambda settings, xs, ys, h5data : background.UniformBackground(settings, xs, ys, h5data),
                   lambda settings, xs, ys, h5data : camera.Ideal(settings, xs, ys, h5data, 100.0),
                   lambda settings, xs, ys, h5data : photophysics.AlwaysOn(settings, xs, ys, h5data, args.photons),
                   lambda settings, xs, ys, h5data : psf.GaussianPSF(settings, xs, ys, h5data, 160.0))

    sim.simulate(args.dax_file, args.hdf5, args.frames)


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
