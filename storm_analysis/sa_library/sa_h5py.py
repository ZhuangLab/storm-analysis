#!/usr/bin/env python
"""
Wraps h5py (http://www.h5py.org/).

Hazen 12/17
"""
import h5py
import numpy
import os
import sys
import time

import storm_analysis.sa_library.grid_c as gridC

# This is the maximum size of a dataset in a group of tracks.
track_block_size = 10000


class SAH5PyException(Exception):
    pass


def isSAHDF5(filename):
    """
    Queries if 'filename' is a storm-analysis HDF5 file.
    """

    # Make sure the file exists.
    if not os.path.exists(filename):
        raise SAH5PyException(filename + " does not exist.")

    # Try and open it and check that it has the expected attributes.
    try:
        with SAH5Py(filename) as h5:
            if ('version' in h5.hdf5.attrs) and ('sa_type' in h5.hdf5.attrs):
                return True
    except OSError:
        pass
    except IOError:
        pass

    return False

def loadLocalizations(filename, fields = None):
    """
    This is convenience function for loading all the localizations in a
    HDF5 file. Only recommended for use on relatively small files.
    """
    with SAH5Reader(filename) as h5:
        locs = h5.getLocalizations(fields = fields)
    return locs

def loadTracks(filename, fields = None):
    """
    This is convenience function for loading all the tracks in a
    HDF5 file. Only recommended for use on relatively small files.
    """
    with SAH5Reader(filename) as h5:
        tracks = h5.getTracks(fields = fields)
    return tracks

def saveLocalizations(filename, locs, frame_number = 0):
    """
    This is convenience function for creating a (simple) HDF5 file from
    a dictionary of localization data.
    """
    with SAH5Py(filename, is_existing = False, overwrite = True) as h5:
        h5.setMovieInformation(1,1,1,"")
        h5.addLocalizations(locs, frame_number)
        

class SAH5Py(object):
    """
    HDF5 file reader/writer.

    The default is behavior is that this will be used to read/modify
    existing HDF5 files. To create new HDF5 files the 'is_existing'
    keyword argument must be set to False.

    Important differences between this format and the old Insight3 format: 
    1. We don't swap the x/y axises on saving. Note that this also applies
       to widths in x/y.
    2. We dropped the single pixel offset in x/y.
    3. We use 0 based frame indexing like the movie reader.
    4. Z is in microns, not nanometers.
    5. These files are required to have pixel size.
    6. Metadata (XML string) is required.
    7. Movie width, height, length are required.

    The internal structure is one group per frame analyzed, with
    each localization property saved as a separate dataset.

    Localizations that have been tracked and averaged together are
    stored as tracks in groups with a maximum dataset size of 
    'track_block_size'.

    The metadata is stored in the 'metadata.xml' root dataset as a 
    variable length unicode string.

    The 'sa_type' attribute records what generated this file, one
    of the SMLM analysis programs for example, or another program
    that for example was used to merge one or more of these files.
    """
    def __init__(self, filename = None, is_existing = True, overwrite = False, sa_type = 'unknown', **kwds):
        super(SAH5Py, self).__init__(**kwds)

        self.last_write_time = time.time()
        self.n_track_groups = 0
        self.total_added = 0

        if is_existing:
            if os.path.exists(filename):
                self.hdf5 = h5py.File(filename, "r+")
                self.existing = True
            else:
                raise SAH5PyException("file '" + filename + "' not found.")

        else:
            if os.path.exists(filename):
                if overwrite:
                    os.remove(filename)
                else:
                    raise SAH5PyException("file '" + filename + "' already exists.")
            self.hdf5 = h5py.File(filename, "w")
            self.hdf5.attrs['analysis_finished'] = 0
            self.hdf5.attrs['n_channels'] = 1
            self.hdf5.attrs['sa_type'] = sa_type
            self.hdf5.attrs['version'] = 0.1
            self.existing = False
            
    def __enter__(self):
        return self

    def __exit__(self, etype, value, traceback):
        if self.hdf5:
            # If this file already existed don't print the number of
            # localizations added.
            if self.existing:
                self.close(verbose = False)
            else:
                self.close(verbose = True)

    def addCategory(self, category, frame_number):
        """
        Add/set the category field of each localization.
        """
        assert isinstance(category, int)
        grp = self.getGroup(frame_number)
        cat = category * numpy.ones(grp.attrs['n_locs'], dtype = numpy.int32)
        self.addLocalizationData(cat, frame_number, "category")

    def addLocalizationData(self, np_data, frame_number, field_name):
        """
        Add/set localization data in an existing group.
        """
        grp = self.getGroup(frame_number)
        assert(np_data.size == grp.attrs['n_locs'])
        
        if not field_name in grp:
            grp.create_dataset(field_name, data = np_data)
        else:
            grp[field_name][()] = np_data
            
    def addLocalizations(self, localizations, frame_number, channel = None):
        """
        Add localization data to the HDF5 file. Each element of localizations
        is stored in it's own dataset. In the case of multi-channel data, the
        data from the other channels is stored in datasets with the prefix 
        'cX_' where X is the channel number. The channel 0 data must be added
        first to create the group.
        """
        if (channel is not None):
            assert(isinstance(channel, int))

            if (channel >= self.hdf5.attrs['n_channels']):
                self.hdf5.attrs['n_channels'] = channel + 1
            
        grp_name = self.getGroupName(frame_number)
        if (channel is None) or (channel == 0):
            grp = self.hdf5.create_group(grp_name)

            # Add initial values for drift correction. These only apply to
            # channel 0.
            grp.attrs['dx'] = 0.0
            grp.attrs['dy'] = 0.0
            grp.attrs['dz'] = 0.0

            # Update counter. Notes:
            #
            # 1. This assumes the existance of the "x" field.
            # 2. This is not necessarily the same as the total number of
            #    localizations in a file as there could for example have
            #    been an analysis restart.
            #
            grp.attrs['n_locs'] = localizations["x"].size
            self.total_added += localizations["x"].size
        else:
            grp = self.getGroup(frame_number)

        for key in localizations:
            d_name = key
            if (channel is not None) and (channel > 0):
                d_name = self.getChannelPrefix(channel) + key
            grp.create_dataset(d_name, data = localizations[key])

        # Flush the file once a minute.
        #
        # FIXME: Not sure if this a bad idea, as for example this might
        #        already be handled in some way by HDF5.
        # 
        current_time = time.time()
        if (current_time > (self.last_write_time + 60.0)):
            self.last_write_time = current_time
            self.hdf5.flush()
            
    def addMetadata(self, metadata):
        """
        Add metadata to the HDF5 file, this is the contents of the XML file
        that was used to analyze the data along with some information about
        the size of the movie.
        """
        assert(isinstance(metadata, str))

        #dt = h5py.special_dtype(vlen = unicode)
        dt = h5py.special_dtype(vlen = str)

        # The +10 was choosen arbitrarily.
        dset_size = (int(len(metadata)+10),)
        dset = self.hdf5.create_dataset("metadata.xml", dset_size, dtype = dt)
        dset[:len(metadata)] = metadata

    def addMovieInformation(self, movie_reader):
        """
        Store some properties of the movie as attributes.
        """
        self.setMovieInformation(movie_reader.getMovieX(),
                                 movie_reader.getMovieY(),
                                 movie_reader.getMovieL(),
                                 movie_reader.hashID())

    def addTrackID(self, track_id, frame_number):
        """
        Add/set the track id field of each localization.
        """
        self.addLocalizationData(track_id, frame_number, "track_id")

    def addTracks(self, tracks):
        """
        Add tracks to the HDF5 file. Tracks are one or more localizations 
        that have been averaged together.

        Note that all the tracks have to be added in a single instantiation.
        If you close this object then the new object will start over and 
        overwrite any existing tracking information.
        """
        # Create tracks group, if necessary.
        if(self.n_track_groups == 0):

            # Delete old tracking information, if any.
            if("tracks" in self.hdf5):
                del self.hdf5["tracks"]
            self.hdf5.create_group("tracks")

        track_grp = self.getTrackGroup()
        grp = track_grp.create_group(self.getTrackGroupName(self.n_track_groups))

        # Add the tracks.
        for field in tracks:
            grp.create_dataset(field, data = tracks[field])
        grp.attrs['n_tracks'] = tracks["x"].size

        self.n_track_groups += 1
        track_grp.attrs['n_groups'] = self.n_track_groups

    def addTrackData(self, np_data, index, field_name):
        """
        Add an additional field to a tracks group.
        """
        track_grp = self.getTrackGroup()
        grp = track_grp[self.getTrackGroupName(index)]
        
        assert(np_data.size == grp.attrs['n_tracks'])

        if not field_name in grp:
            grp.create_dataset(field_name, data = np_data)
        else:
            grp[field_name][()] = np_data
        
    def addLocalizationZ(self, z_vals, frame_number):
        """
        Add/set the z field of each localization.
        """
        self.addLocalizationData(z_vals, frame_number, "z")

    def close(self, verbose = False):
        if verbose:
            print("Added", self.total_added)
        self.hdf5.close()

    def getChannelPrefix(self, channel_number):
        return "c" + str(channel_number) + "_"
        
    def getDatasets(self, group, fields):
        datasets = {}
        
        # Return all the datasets in the group.
        if fields is None:
            for field in group:
                datasets[field] = group[field][()]

        # Return only the fields that the user requested.
        else:
            for field in fields:
                datasets[field] = group[field][()]

        return datasets
        
    def getDriftCorrection(self, frame_number):
        grp = self.getGroup(frame_number)
        if grp is not None:
            return [grp.attrs['dx'], grp.attrs['dy'], grp.attrs['dz']]
        else:
            raise SAH5PyException("getDriftCorrection(), no such frame " + str(frame_number))
        
    def getFileType(self):
        return self.hdf5.attrs['sa_type']
    
    def getFileVersion(self):
        return self.hdf5.attrs['version']

    def getGroup(self, frame_number):
        grp_name = self.getGroupName(frame_number)
        if grp_name in self.hdf5:
            return self.hdf5[grp_name]

    def getGroupName(self, frame_number):
        assert(isinstance(frame_number, int))
        return "/fr_" + str(frame_number)

    def getLocalizations(self, drift_corrected = False, fields = None):
        return self.getLocalizationsInFrameRange(0,
                                                 self.hdf5.attrs['movie_l'],
                                                 drift_corrected = drift_corrected,
                                                 fields = fields)
    
    def getLocalizationsInFrame(self, frame_number, drift_corrected = False, fields = None):

        locs = {}
        grp = self.getGroup(frame_number)
        
        if (grp is not None) and (grp.attrs['n_locs'] > 0):
            locs = self.getDatasets(grp, fields)

        if drift_corrected and bool(locs):
            if "x" in locs:
                locs["x"] += grp.attrs['dx']
            if "y" in locs:
                locs["y"] += grp.attrs['dy']
            if "z" in locs:
                locs["z"] += grp.attrs['dz']

        return locs

    def getLocalizationsInFrameRange(self, start, stop, drift_corrected = False, fields = None):
        """
        Return the localizations in the range start <= frame number < stop.
        """
        assert(stop > start)
        locs = {}
        for i in range(start, stop):
            temp = self.getLocalizationsInFrame(i,
                                                fields = fields,
                                                drift_corrected = drift_corrected)
            if(not bool(temp)):
                continue
            
            for field in temp:
                if field in locs:
                    locs[field] = numpy.concatenate((locs[field], temp[field]))
                else:
                    locs[field] = temp[field]

        return locs
        
    def getMetadata(self):
        if "metadata.xml" in self.hdf5:
            return self.hdf5["metadata.xml"][0]
        raise SAH5PyException("No metadata!")

    def getMovieInformation(self):
        """
        Return the dimensions of the movie and it's hash value ID.
        """
        if 'movie_l' in self.hdf5.attrs:
            return [int(self.hdf5.attrs['movie_x']),
                    int(self.hdf5.attrs['movie_y']),
                    int(self.hdf5.attrs['movie_l']),
                    self.hdf5.attrs['movie_hash_value']]

        raise SAH5PyException("Movie information is not available! Maybe no frames were analyzed?")

    def getMovieLength(self):
        """
        Return the length of the movie.
        """
        if 'movie_l' in self.hdf5.attrs:
            return int(self.hdf5.attrs['movie_l'])

        raise SAH5PyException("Movie length is not available! Maybe no frames were analyzed?")

    def getNChannels(self):
        """
        Return the number of channels.
        """
        return(int(self.hdf5.attrs['n_channels']))
        
    def getNLocalizations(self):
        n_locs = 0
        for i in range(self.getMovieLength()):
            grp = self.getGroup(i)
            if(grp is not None):
                n_locs += grp.attrs['n_locs']
        return n_locs

    def getNTracks(self):
        if(not self.hasTracks()):
            return 0
        track_grp = self.getTrackGroup()
        n_tracks = 0
        for i in range(track_grp.attrs['n_groups']):
            n_tracks += track_grp[self.getTrackGroupName(i)].attrs['n_tracks']
        return n_tracks

    def getPixelSize(self):
        if 'pixel_size' in self.hdf5.attrs:
            return float(self.hdf5.attrs['pixel_size'])
        raise SAH5PyException("No pixel size information!")    

    def getTracksByIndex(self, index, fields):
        """
        Mostly for internal use. See trackIterator() for the recommended way
        to access the tracks.
        """

        # Check that there are tracks.
        if (not self.hasTracks()):
            return {}
        track_grp = self.getTrackGroup()
        t_grp_name = self.getTrackGroupName(index)

        # Check that the requested group of tracks exists.
        if (not t_grp_name in track_grp):
            return {}

        # Return the tracks.
        return self.getDatasets(track_grp[t_grp_name], fields)
    
    def getTrackGroup(self):
        return self.hdf5["tracks"]
        
    def getTrackGroupName(self, index):
        return "tracks_" + str(index)

    def getTracks(self, fields = None):
        """
        This will return all the tracks in the file in a single dictionary. Not
        recommended for use with large files..
        """
        all_tracks = {}
        for tracks in self.tracksIterator(fields = fields):
            for field in tracks:
                if field in all_tracks:
                    all_tracks[field] = numpy.concatenate((all_tracks[field], tracks[field]))
                else:
                    all_tracks[field] = tracks[field]

        return all_tracks
                        
    def hasLocalizationsField(self, field_name):
        """
        Return True if the localizations have the dataset 'field_name'.
        """
        for fnum, locs in self.localizationsIterator():
            if field_name in locs:
                return True
            else:
                return False

    def hasTracks(self):
        return ("tracks" in self.hdf5)

    def hasTracksField(self, field_name):
        """
        Return True if the tracks have the dataset 'field_name'.

        FIXME: We probably don't need to actually load 10k tracks just to
               determine what fields are available.
        """
        for tracks in self.tracksIterator():
            if field_name in tracks:
                return True
            else:
                return False

    def isAnalysisFinished(self):
        return (self.hdf5.attrs['analysis_finished'] != 0)
            
    def isAnalyzed(self, frame_number):
        return self.getGroup(frame_number) is not None
        
    def isExisting(self):
        """
        Return TRUE if the underlying HDF5 file already existed, FALSE if
        we just created it.
        """
        return self.existing

    def localizationsIterator(self, drift_corrected = True, fields = None, skip_empty = True):
        """
        An iterator for getting all the localizations in a for loop. This is
        probably the easiest way to process all the localizations in a file.
        It returns a two-element list, [frame number, localizations], where
        localizations is dictionary containing the localization properities
        as numpy arrays.

        for fnum, locs in h5.localizationsIterator():
            ..

        If skip_empty is false you will get empty locs dictionaries for frames
        that have no localizations.
        """
        # This should be a zero length generator if there are no localizations.
        for i in range(self.getMovieLength()):
            locs = self.getLocalizationsInFrame(i,
                                                drift_corrected = drift_corrected,
                                                fields = fields)
            
            # Check if locs is empty and we should skip.
            if not bool(locs) and skip_empty:
                continue
            else:
                yield [i, locs]

    def setAnalysisFinished(self, finished):
        if finished:
            self.hdf5.attrs['analysis_finished'] = 1
        else:
            self.hdf5.attrs['analysis_finished'] = 0            
    
    def setDriftCorrection(self, frame_number, dx = 0.0, dy = 0.0, dz = 0.0):
        grp = self.getGroup(frame_number)
        if grp is not None:
            grp.attrs['dx'] = dx
            grp.attrs['dy'] = dy
            grp.attrs['dz'] = dz
        else:
            raise SAH5PyException("setDriftCorrection(), no such frame " + str(frame_number))

    def setMovieInformation(self, movie_x, movie_y, movie_l, hash_value):
        """
        Store some properties of the movie as attributes.
        """
        self.hdf5.attrs['movie_hash_value'] = hash_value
        self.hdf5.attrs['movie_l'] = movie_l
        self.hdf5.attrs['movie_x'] = movie_x
        self.hdf5.attrs['movie_y'] = movie_y

    def setPixelSize(self, pixel_size):
        """
        Add pixel size in information (in nanometers).
        """
        self.hdf5.attrs['pixel_size'] = pixel_size

    def splitByChannel(self, locs):
        """
        Split a dictionary of localizations by channel.
        """
        assert bool(locs)

        # Create a list of dictionaries.
        split_locs = [{} for i in range(self.getNChannels())]

        for elt in locs:
            used = False
            for i in range(1, self.getNChannels()):
                prefix = self.getChannelPrefix(i)
                if elt.startswith(prefix):
                    split_locs[i][elt[len(prefix):]] = locs[elt]
                    used = True
                    break

            if not used:
                split_locs[0][elt] = locs[elt]
        
        return split_locs
        
    def tracksIterator(self, fields = None):
        """
        An iterator for getting all the tracks in a for loop. This approach
        is preferred over loading all the tracks into memory at once. You
        can optionally specify that it only return certain fields.

        for tracks in h5.tracksIterator():
            ..

        The track center in x,y,z are the 'x', 'y' and 'z' fields. The other
        fields may need to be normalized by the track length.
        """
        # This should be a zero length generator if there are no tracks.
        if (not self.hasTracks()):
            for i in range(0):
                yield {}

        else:
            track_grp = self.getTrackGroup()
            for i in range(track_grp.attrs['n_groups']):
                yield self.getTracksByIndex(i, fields)


class SAH5Grid(SAH5Py):
    """
    A sub-class of SAH5Py that has the ability to grid the tracks. This
    is similar to sa_library.drift_utilities.SAH5DriftCorrection.
    """
    def __init__(self, scale = None, z_bins = 1, **kwds):
        super(SAH5Grid, self).__init__(**kwds)
        
        self.im_shape_2D = (self.hdf5.attrs['movie_x']*scale,
                            self.hdf5.attrs['movie_y']*scale)
        self.im_shape_3D = (self.hdf5.attrs['movie_x']*scale,
                            self.hdf5.attrs['movie_y']*scale,
                            z_bins)
        self.scale = scale
        self.z_bins = z_bins

    def gridTracks2D(self, dx = 0.0, dy = 0.0, verbose = False):
        image = numpy.zeros(self.im_shape_2D, dtype = numpy.int32)
        for locs in self.tracksIterator(fields = ["x", "y"]):
            if verbose:
                sys.stdout.write(".")
                sys.stdout.flush()

            f_x = locs["x"] + dx
            f_y = locs["y"] + dy
            
            i_x = numpy.floor(f_x*self.scale).astype(numpy.int32)
            i_y = numpy.floor(f_y*self.scale).astype(numpy.int32)
            gridC.grid2D(i_x, i_y, image)
            
        if verbose:
            sys.stdout.write("\n")
            
        return image

    def gridTracks3D(self, z_min, z_max, dx = 0.0, dy = 0.0, verbose = False):
        z_scale = float(self.z_bins)/(z_max - z_min)
        image = numpy.zeros(self.im_shape_3D, dtype = numpy.int32)
        for locs in self.tracksIterator(fields = ["x", "y", "z"]):
            if verbose:
                sys.stdout.write(".")
                sys.stdout.flush()
                
            # Add to image.
            f_x = locs["x"] + dx
            f_y = locs["y"] + dy
            f_z = (locs["z"] - z_min)*z_scale
            
            i_x = numpy.floor(f_x*self.scale).astype(numpy.int32)
            i_y = numpy.floor(f_y*self.scale).astype(numpy.int32)
            i_z = numpy.floor(f_z.astype(numpy.int32))
            gridC.grid3D(i_x, i_y, i_z, image)

        if verbose:
            sys.stdout.write("\n")
            
        return image


class SAH5Reader(SAH5Py):
    """
    HDF5 file read only access.

    Use this if you only want to read the localization file. It won't
    lock the file restricting it to a single process.
    """
    def __init__(self, filename = None):

        # This used be close()
        self.total_added = 0
        
        if os.path.exists(filename):
            self.hdf5 = h5py.File(filename, "r")
            self.existing = True
        else:
            raise SAH5PyException("file '" + filename + "' not found.")        
        

        
if (__name__ == "__main__"):

    import os
    import sys

    from xml.etree import ElementTree

    if (len(sys.argv) != 2):
        print("usage: <hdf5_file>")
        exit()

    if not os.path.exists(sys.argv[1]):
        print("File", sys.argv[1], "not found.")
        exit()
        
    with SAH5Reader(sys.argv[1]) as h5:
        try:
            metadata = h5.getMetadata()
        except SAH5PyException as ex:
            print(ex)
            metadata = None
            
        if metadata is not None:
            metadata = ElementTree.fromstring(metadata)

            print(" meta data:")
            for node in sorted(metadata, key = lambda node: node.tag):
                print("    " + node.tag.strip() + " - " + node.text.strip())

        print()
        [mx, my, ml, mh] = h5.getMovieInformation()
        print("Size:", mx, "x", my, "pixels")
        print("Frames:", h5.getMovieLength())
        print("Localizations:", h5.getNLocalizations())
        print("Tracks:", h5.getNTracks())

        if(h5.getNLocalizations() > 0):
            print()
            print("Localization statistics:")

            locs = h5.getLocalizations()
            print(locs["x"].size)
            for field in locs:
                print("  {0:15} {1:.3f} {2:.3f} {3:.3f} {4:.3f}".format(field,
                                                                        numpy.mean(locs[field]),
                                                                        numpy.std(locs[field]),
                                                                        numpy.min(locs[field]),
                                                                        numpy.max(locs[field])))
