#!/usr/bin/env python
"""
Wraps h5py (http://www.h5py.org/).

Hazen 12/17
"""
import h5py
import numpy
import os


class SAH5Py(object):
    """
    HDF5 file reader/writer.

    Unlike the Insight3 format I believe we can do both at the same
    time with HDF5.

    The internal structure is one group per frame analyzed, with
    each localization property saved as a separate dataset.
   
    The metadata is stored in the 'metadata.xml' root dataset as a 
    variable length unicode string.
    """
    def __init__(self, filename = None, **kwds):
        super(SAH5Py, self).__init__(**kwds)

        if os.path.exists(filename):
            self.hdf5 = h5py.File(filename, "r+")
        else:
            self.hdf5 = h5py.File(filename, "w")
            self.hdf5.attrs['version'] = 0.1
            
    def __enter__(self):
        return self

    def __exit__(self, etype, value, traceback):
        if self.hdf5:
            self.hdf5.close()

    def addLocalizations(self, localizations, frame_number):
        pass
    
    def addMetadata(self, metadata):
        """
        Add metadata to the HDF5 file, this is the contents of the XML file
        that was used to analyze the data along with some information about
        the size of the movie.
        """
        assert(isinstance(metadata, str))
        assert(not "metadata.xml" in self.hdf5)

        #dt = h5py.special_dtype(vlen = unicode)
        dt = h5py.special_dtype(vlen = str)

        # The +10 was choosen arbitrarily.
        dset_size = (int(len(metadata)+10),)
        dset = self.hdf5.create_dataset("metadata.xml", dset_size, dtype = dt)
        dset[:len(metadata)] = metadata
        self.hdf5.flush()

    def close(self):
        self.hdf5.close()

    def getFileVersion(self):
        return self.hdf5.attrs['version']
    
    def getMetadata(self):
        return self.hdf5["metadata.xml"][0]
    
