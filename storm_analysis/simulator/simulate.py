#!/usr/bin/python
#
# Generate simulated. The basic idea is that you provide a list of
# localizations in Insight3 .bin format and these are used to
# generate a series of images using the following steps:
#
# Initialization:
#   1. locs = readinsight3.loadI3File(mlist_file)
#   2. bg = background.Background(locs)
#   3. camera = camera.Camera()
#   4. pp = photophysics.Photophysics(locs, intensity)
#   5. psf = psf.PSF()
#
#
# Generation:
#   1. image = numpy.zeros((256, 256))
#   2. image += bg.getBackground()
#   3. cur_locs = pp.getEmitters()
#   4. image += psf.getPSFs(cur_locs)
#   5. image = camera.readImage(image)
#   6. saveimage(image)
#   7. saveloc(cur_locs, frame_number)
#
# Note: This expected to set the 'h' and 'bg' fields in the
#       output list to the correct values. The values 'x', 'y',
#       'z', 'ax' and 'w' are just passed through. Other values
#       such as 'a' and 'i' are not set and are will likely
#       be incorrect.
#
#
# Hazen 11/16
#

import json
import numpy

import storm_analysis.sa_library.daxwriter as daxwriter
import storm_analysis.sa_library.readinsight3 as readinsight3
import storm_analysis.sa_library.writeinsight3 as writeinsight3


class Simulate(object):

    def __init__(self, background_factory, camera_factory, photophysics_factory, psf_factory, x_size = 256, y_size = 256):
        """
        The factor variables should be functions that return the correct class
        to run a simulation with the following signature:

        factory_fn(sim_settings, x_size, y_size, i3_data_in)

        """
        self.bg_factory = background_factory
        self.cam_factory = camera_factory
        self.pphys_factory = photophysics_factory
        self.psf_factory = psf_factory
        self.x_size = x_size
        self.y_size = y_size

    def setBackgroundFactory(self, new_factory):
        self.bg_factory = new_factory

    def setCameraFactory(self, new_factory):
        self.cam_factory = new_factory

    def setPhotoPhysicsFactory(self, new_factory):
        self.pphys_factory = new_factory

    def setPSFFactory(self, new_factory):
        self.psf_factory = new_factory

    def simulate(self, dax_file, bin_file, n_frames):

        #
        # Initialization.
        #
        dax_data = daxwriter.DaxWriter(dax_file, self.x_size, self.y_size)
        i3_data_in = readinsight3.loadI3File(bin_file)

        out_fname_base = dax_file[:-4]
        i3_data_out = writeinsight3.I3Writer(out_fname_base + "_olist.bin")
        sim_settings = open(out_fname_base + "_sim_params.txt", "w")

        #
        # Create the user-specified class instances that will do
        # most of the actual work of the simulation.
        #
        bg = self.bg_factory(sim_settings, self.x_size, self.y_size, i3_data_in)
        cam = self.cam_factory(sim_settings, self.x_size, self.y_size, i3_data_in)
        pp = self.pphys_factory(sim_settings, self.x_size, self.y_size, i3_data_in)
        psf = self.psf_factory(sim_settings, self.x_size, self.y_size, i3_data_in)

        sim_settings.write(json.dumps({"simulation" : {"bin_file" : bin_file,
                                                       "x_size" : str(self.x_size),
                                                       "y_size" : str(self.y_size)}}) + "\n")

        #
        # Generate the simulated movie.
        #
        for i in range(n_frames):
            print("Generating frame:", i)

            # Generate the new image.
            image = numpy.zeros((self.x_size, self.y_size))
            cur_i3 = pp.getEmitters(i)

            # Background
            image += bg.getBackground(i)
            cur_i3 = bg.getEmitterBackground(cur_i3)

            # Foreground
            image += psf.getPSFs(cur_i3)

            # Camera
            image = cam.readImage(image)

            # Save the image.
            dax_data.addFrame(image)

            # Save the molecule locations.
            cur_i3['fr'] = i + 1
            i3_data_out.addMolecules(cur_i3)

        dax_data.close()
        i3_data_out.close()
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
    parser.add_argument('--bin', dest='i3bin', type=str, required=True,
                        help = "The name of the Insight3 file containing the emitter locations.")
    parser.add_argument('--frames', dest='frames', type=int, required=True,
                        help = "The length of the movie in frames.")
    parser.add_argument('--height', dest='height', type=float, required=True,
                        help = "The peak height of a single emitter in photons.")

    args = parser.parse_args()
 
    sim = Simulate(lambda settings, xs, ys, i3data : background.UniformBackground(settings, xs, ys, i3data),
                   lambda settings, xs, ys, i3data : camera.Ideal(settings, xs, ys, i3data, 100.0),
                   lambda settings, xs, ys, i3data : photophysics.AlwaysOn(settings, xs, ys, i3data, args.height),
                   lambda settings, xs, ys, i3data : psf.GaussianPSF(settings, xs, ys, i3data, 160.0),
                   256,
                   256)

    sim.simulate(args.dax_file, args.i3bin, args.frames)


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
