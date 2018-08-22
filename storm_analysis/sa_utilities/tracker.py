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

    def __init__(self, category = None, frame_number = None, track_id = None, **kwds):
        super(Track, self).__init__(**kwds)

        self.category = category
        self.frame_number = frame_number
        self.last_added = 0
        self.length = 0
        self.props = {}
        self.track_id = track_id
        self.tw = 0.0
        self.tx = 0.0
        self.ty = 0.0
        self.tz = 0.0

    def addLocalization(self, loc, index):
        """
        Add the localization to the track.
        """
        self.length += 1

        # Update track center.
        #
        # FIXME: Weight values based on Mortensen error estimate? The current
        #        weighting is just a simple sqrt().
        #
        w = 1
        if "sum" in loc:
            w += math.sqrt(loc["sum"][index])

        self.tw += w
        self.tx += w*loc["x"][index]
        self.ty += w*loc["y"][index]
        if "z" in loc:
            self.tz += w*loc["z"][index]
        self.last_added = 0

        # Record the localization properties.
        for key in loc:
            if key in self.props:
                self.props[key] += loc[key][index]
            else:
                self.props[key] = loc[key][index]

    def incLastAdded(self):
        self.last_added += 1
        
    def getCenter(self):
        """
        Return current track center.
        """
        return [self.tx/self.tw, self.ty/self.tw]

    def getLastAdded(self):
        return self.last_added

    def getProperties(self):
        return self.props


class TrackWriter(object):
    """
    This handles saving the tracks in the HDF5 file in blocks.

    FIXME: Automatically normalize those fields for which this is appropriate?
           Skip values that it makes no sense to include in the tracks?
    """
    def __init__(self, sa_h5 = None, **kwds):
        super(TrackWriter, self).__init__(**kwds)

        self.data = {}
        self.h5 = sa_h5
        self.n_data = 0
        self.n_tracks = 0

    def finish(self, verbose = True):
        if(self.n_data > 0):
            temp = {}
            for field in self.data:
                temp[field] = self.data[field][:self.n_data]
            self.h5.addTracks(temp)

        print("Added", self.n_tracks, "tracks")

    def writeTrack(self, track):

        # Get the track's properties.
        t_props = track.getProperties()

        # Check if we need to initialize self.data.
        if not bool(self.data):
            for field in t_props:
                if (t_props[field].dtype == numpy.int32):
                    self.data[field] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.int32)
                else:
                    self.data[field] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.float64)
            self.data["category"] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.int32)
            self.data["frame_number"] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.int32)
            self.data["track_id"] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.int64)
            self.data["track_length"] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.int32)
            self.data["x"] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.float64)
            self.data["y"] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.float64)
            self.data["z"] = numpy.zeros(saH5Py.track_block_size, dtype = numpy.float64)

        # Add track to data.
        for field in t_props:
            self.data[field][self.n_data] = t_props[field]
        self.data["category"][self.n_data] = track.category
        self.data["frame_number"][self.n_data] = track.frame_number
        self.data["track_id"][self.n_data] = track.track_id
        self.data["track_length"][self.n_data] = track.length
        self.data["x"][self.n_data] = track.tx/track.tw
        self.data["y"][self.n_data] = track.ty/track.tw
        self.data["z"][self.n_data] = track.tz/track.tw

        # Increment counters.
        self.n_data += 1
        self.n_tracks += 1

        # Add tracks to the HDF5 file if we have a full block.
        if (self.n_data == saH5Py.track_block_size):
            self.h5.addTracks(self.data)
            self.n_data = 0
    

def tracker(sa_hdf5_filename, descriptor = "", max_gap = 0, radius = 0.0):
    """
    descriptor - A string containing the frame designation.
    max_gap - The maximum number of frames with no objects before a track is 
              considered to be terminated.
    radius - Maximum distance for an object to be in a track in pixels.
    """
    # Just set localization category if radius is zero or negative.
    if (radius <= 0.0):
        with saH5Py.SAH5Py(sa_hdf5_filename) as h5:
            for fnum, locs in h5.localizationsIterator():

                # Determine current frame description.
                fdesc = 1
                if(len(descriptor) > 0):
                    fdesc = int(descriptor[(fnum%len(descriptor))])
                
                # The category is the descriptor minus 1.
                category = fdesc - 1

                h5.addCategory(category,fnum)

    # Otherwise do the tracking.
    else:
        track_id = 0
        current_tracks = []
        with saH5Py.SAH5Py(sa_hdf5_filename) as h5:
            tw = TrackWriter(h5)
            for fnum, locs in h5.localizationsIterator(skip_empty = False):

                # User feedback.
                if ((fnum%500)==0):
                    print(" processing frame {0:0d}, {1:0d} tracks".format(fnum, track_id))

                # Determine current frame description.
                fdesc = 1
                if(len(descriptor) > 0):
                    fdesc = int(descriptor[(fnum%len(descriptor))])
                
                # The category is the descriptor minus 1.
                category = fdesc - 1

                # Add/update localization category.
                if bool(locs):
                    h5.addCategory(category,fnum)

                # Go to the next frame if this is an activation frame.
                if(fdesc == 0):
                    continue

                # Check that the frame had localizations, assign them if it did.
                index_locs = None
                locs_track_id = None
                if bool(locs):

                    # Create numpy array for storage of the track id for each localization.
                    locs_track_id = numpy.zeros(locs["x"].size, dtype = numpy.int64)
                
                    # Create arrays with current track centers. This is also increments
                    # the tracks last added counter.
                    tx = numpy.zeros(len(current_tracks))
                    ty = numpy.zeros(len(current_tracks))
                    for i, elt in enumerate(current_tracks):
                        [tx[i], ty[i]] = elt.getCenter()
                        elt.incLastAdded()
                
                    kd_locs = iaUtilsC.KDTree(locs["x"], locs["y"])
                    kd_tracks = iaUtilsC.KDTree(tx, ty)
                
                    # Query KD trees.
                    index_locs = kd_tracks.nearest(locs["x"], locs["y"], radius)[1]
                    index_tracks = kd_locs.nearest(tx, ty, radius)[1]
                
                    # Add localizations to tracks. The localization must be the closest
                    # one to the track and vice-versa. We're trying to avoid multiple
                    # localizations in a single frame in the track, and one localization
                    # in multiple tracks.
                    #
                    for i in range(locs["x"].size):
                        if (index_locs[i] > -1):
                            
                            # Check that the track and the localization agree that each
                            # is closest to the other.
                            #
                            tr = None
                            if (index_tracks[index_locs[i]] == i):
                                tr = current_tracks[index_locs[i]]
                                tr.addLocalization(locs, i)
                            else:
                                tr = Track(category = category,
                                           frame_number = fnum,
                                           track_id = track_id)
                                tr.addLocalization(locs, i)
                                current_tracks.append(tr)
                                track_id += 1

                            locs_track_id[i] = tr.track_id

                    # Clean up KD trees.
                    kd_locs.cleanup()
                    kd_tracks.cleanup()

                # Otherwise just increment the current tracks last added counter.
                else:
                    for elt in current_tracks:
                        elt.incLastAdded()

                # Remove tracks that have not had any localizations added for
                # max_gap frames.
                temp = current_tracks
                current_tracks = []
                for elt in temp:
                    if(elt.getLastAdded() > max_gap):
                        tw.writeTrack(elt)
                    else:
                        current_tracks.append(elt)

                # Start new tracks from the localizations that were not in
                # a track.
                if index_locs is not None:
                    for i in range(locs["x"].size):
                        if (index_locs[i] < 0):
                            tr = Track(category = category,
                                       frame_number = fnum,
                                       track_id = track_id)
                            tr.addLocalization(locs, i)
                            current_tracks.append(tr)
                            locs_track_id[i] = tr.track_id
                            track_id += 1

                # Add track information for localizations.
                if locs_track_id is not None:
                    h5.addTrackID(locs_track_id, fnum)

            # Write the remaining tracks & close the track writer.
            for elt in current_tracks:
                tw.writeTrack(elt)

            tw.finish()
