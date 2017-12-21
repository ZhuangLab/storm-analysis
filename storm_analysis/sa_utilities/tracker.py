#!/usr/bin/env python
"""
Tracks (and averages) localizations across frames in a SMLM movie.

Notes:
1. Drift correction must be done before tracking.

2. This is different from the original in that localizations are 
   assigned to the nearest track, not to all of the tracks that 
   they are close enough to.

3. This does the averaging at the same time that it does the
   tracking.

Hazen 12/17
"""
import math
import numpy

import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py


class Track(object):
    """
    Class to store & manipulate a single track.
    """

    def __init__(self, category = None, **kwds):
        super(Track, self).__init__(**kwds)

        self.category = category
        self.last_added = 0
        self.length = 0
        self.tw = 0.0
        self.tx = None
        self.ty = None
        
    def addLocalization(self, loc):
        """
        Add localization to the track.
        """
        if "sum" in loc:
            self.tw += math.sqrt(loc["sum"])
        else:
            self.tw += 1
        self.tx += loc["x"]
        self.ty += loc["y"]
        self.last_added = 0

    def incLastAdded(self):
        self.last_added += 1
        
    def getCenter(self):
        """
        Return current track center.
        """
        return [self.tx/self.tw, self.ty/self.tw]

    def getLastAdded(self):
        return self.last_added


def tracker(sa_hdf5_filename, descriptor = "", radius = 0.0):
    """
    descriptor - a string containing the frame designation.
    radius - maximum distance for an object to be in a track in pixels.
    """

    fnum = 0
    current_tracks = []
    with saH5Py.SAH5Py(sa_hdf5_filename) as h5:
        tw = TrackWriter(h5)
        while(fnum < h5.getMovieLength()):

            # Determine current frame description.
            fdesc = 1
            if(len(descriptor) > 0):
                fdesc = int(descriptor[(fnum%len(descriptor))])
                
            # Go to the next frame if this is an activation frame.
            if(fdesc == 0):
                fnum += 1
                continue

            # The category is the descriptor minus 1.
            category = fdesc - 1

            # Load the localization x and y positions.
            locs_xy = h5.getLocalizationsInFrame(fnum,
                                                 drift_corrected = True,
                                                 fields = ["x", "y"])

            # If there are no localizations, terminate all the current tracks.
            if not bool(locs_xy):
                for elt in current_tracks:
                    tw.writeTrack(elt)
                current_tracks = []

            # Create KD tree from current tracks. This is also increments
            # the tracks last added counter.
            tx = numpy.zeros(len(current_tracks))
            ty = numpy.zeros(len(current_tracks))
            for i, elt in enumerate(current_tracks):
                [tx[i], ty[i]] = elt.getCenter()
                elt.incLastAdded()
                
            kd = iaUtilsC.KDTree(tx, ty)

            # Query with localizations.
            index = kd.nearest(locs_xy["x"], locs_xy["y"], radius)

            # Add localizations to tracks.
            for i in range(locs_xy["x"].size):
                if (index[i] > -1):
                    current_tracks[index[i]].addLocalization(...)

            # Remove tracks that did not have any localizations added.
            temp = current_tracks
            current_tracks = []
            for elt in temp:
                if(elt.getLastAdded != 0):
                    tw.writeTrack(elt)
                else:
                    current_tracks.append(elt)

            # Start new tracks from the localizations that were not in a track.
            for i in range(locs_xy["x"].size):
                if (index[i] == -1):
                    tr = Track(category = category)
                    tr.addLocalization(...)
                    current_tracks.append(tr)

        # Write the remaining tracks & close the track writer.
        for elt in current_tracks:
            tw.writeTrack(elt)
            
                
