#!/usr/bin/env python
"""
Tracks fiducials across frames in a SMLM movie.

This adds the "fiducial_id" to field to each localization. Numbers
less than zero in this field mean that the localization was not
assigned to a fiducial.

Hazen 09/18
"""
import math
import numpy

import storm_analysis
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysys.sa_utilities.tracker as tracker


class FiducialException(storm_analysis.SAException):
    pass


class Fiducial(tracker.Track):
    """
    Class to store & manipulate a single fiducial.
    """
    def addLocalization(self, loc, index):
        """
        Add the localization to the fiducial.
        """
        self.length += 1

        self.tx = loc["x"][index]
        self.ty = loc["y"][index]
        self.last_added = 0
        
    def getCenter(self):
        """
        Return current track center.
        """
        return [self.tx, self.ty]

    def getFiducialId(self):
        return self.track_id


def trackFiducials(sa_hdf5_filename, max_gap = 0, radius = 0.0, reference_frame = 0):
    """
    max_gap - The maximum number of frames with no objects before a track is 
              considered to be terminated.
    radius - Maximum distance for an object to be in a track in pixels.
    reference_frame - The localizations in this frame will be used as the
                      fiducial reference locations.
    """

    with saH5Py.SAH5Py(sa_hdf5_filename) as h5:
        ref_locs = h5.getLocalizationsInFrame(reference_frame, fields = ["x", "y"])

        if not ref_locs:
            raise FiducialException("No localizations in frame " + str(reference_frame))

        # Create fiducials list.
        fiducials = []
        for i in range(ref_locs["x"].size):
            fd = Fiducial(track_id = i)
            fd.addLocalization(ref_locs, i)
            fiducials.append(df)

        # Iterate over localizations.
        for fnum, locs in h5.localizationsIterator(skip_empty = False, fields = ["x", "y"]):

            # User feedback.
            if ((fnum%500)==0):
                print(" processing frame {0:0d}, {1:0d} fiducials".format(fnum, len(fiducials)))

            # Check if we still have any 'live' fiducials.
            if (len(fiducials) == 0):

                # Mark all localizations as non-fiducial, i.e. -1.
                if bool(locs):
                    locs_fd_id = numpy.zeros(locs["x"].size, dtype = numpy.int32) - 1
                    h5.addLocalizationData(locs_fd_id, fnum, "fiducial_id")

                continue
                    
            # Check that the frame had localizations, assign them to fiducials if it did.
            if bool(locs):

                # Create numpy array for storage of the fiducial id for each localization.
                locs_fd_id = numpy.zeros(locs["x"].size, dtype = numpy.int32) - 1
                
                # Create arrays with current fiducial centers. This is also increments
                # the fiducials last added counter.
                tx = numpy.zeros(len(fiducials))
                ty = numpy.zeros(len(fiducials))
                for i, elt in enumerate(fiducials):
                    [tx[i], ty[i]] = elt.getCenter()
                    elt.incLastAdded()
                
                kd_locs = iaUtilsC.KDTree(locs["x"], locs["y"])
                
                # Query for nearest localization to each fiducial.
                index_fd = kd_locs.nearest(tx, ty, radius)[1]

                # Add fiducial information to localizations.
                for i in range(index_fd.size):
                    if (index_fd[i] > -1):
                        locs_fd_id[index_fd[i]] = fiducials[index_fd[i]].getFiducialId()

                h5.addLocalizationData(locs_fd_id, fnum, "fiducial_id")
                                          
                # Clean up KD trees.
                kd_locs.cleanup()
            
            # Otherwise just increment the fiducials last added counter.
            else:
                for elt in fiducials:
                    elt.incLastAdded()

            # Remove fiducials that have not had any localizations added for max_gap frames.
            temp = fiducials
            fiducials = []
            for elt in temp:
                if(elt.getLastAdded() < max_gap):
                    fiducials.append(elt)
