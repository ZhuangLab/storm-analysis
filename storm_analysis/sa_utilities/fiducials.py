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
import scipy.interpolate

import storm_analysis
import storm_analysis.sa_library.ia_utilities_c as iaUtilsC
import storm_analysis.sa_library.sa_h5py as saH5Py
import storm_analysis.sa_utilities.tracker as tracker


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


class SAH5Fiducials(saH5Py.SAH5Py):
    """
    A sub-class of SAH5Py for working with fiducials.
    """
    def addFiducialID(self, fiducial_id, frame_number):
        """
        Add/set the fiducial id field of each localization.
        """
        self.addLocalizationData(fiducial_id, frame_number, "fiducial_id")

    def averageFiducials(self, drift_corrected = False, fields = None, min_frac_occupancy = 0.9, preload_all = True):
        """
        Returns a dictionary with the requested fields average across all
        the frames. Missing values are handled by linear interpolation.

        drift_corrected - Use drift corrected values for fiducial locations.
        fields - (Localization) fields to include in average, None = all of them.
        min_frac_occupancy - The minimum fraction of the frames that a 
                             fiducial needs to be present to be included
                             in the average.
        preload_all - Load all the fiducials at initialization.

        return [averaged_fields, number_averaged]
        """
        ml = self.getMovieLength()
        
        ave_fdcl = {}
        fr = numpy.arange(ml)
        min_occ = int(min_frac_occupancy * ml)
        n_fiducials = 0
        for fdcl in self.fiducialsIterator(fields = fields, drift_corrected = drift_corrected):
            
            # Skip fiducials with no data. Can this actually happen? I
            # think there would always be at least one element.
            if not fdcl:
                continue
            
            # Verify this fiducial was present in most of the frames.
            if (fdcl["frame"].size < min_occ):
                continue

            n_fiducials += 1
            for elt in fdcl:
                if (elt != "frame"):
        
                    # Use linear interpolation to fill in missing values.
                    lin_int = scipy.interpolate.interp1d(fdcl["frame"],
                                                         fdcl[elt],
                                                         kind = 'linear',
                                                         fill_value = 'extrapolate')

                    if not elt in ave_fdcl:
                        ave_fdcl[elt] = lin_int(fr)
                    else:
                        ave_fdcl[elt] += lin_int(fr)

        if ave_fdcl:
            for elt in ave_fdcl:
                ave_fdcl[elt] = ave_fdcl[elt]/float(n_fiducials)

        return [ave_fdcl, n_fiducials]

    def getFiducial(self, f_id, drift_corrected = False, fields = None):
        """
        Get the data for a single fiducial.
        """
        # Add 'fiducial_id' as we'll also need this field.
        if fields is not None:
            fields.append("fiducial_id")
            
        return self.getFiducialData(f_id, drift_corrected, fields)
        
    def getFiducialData(self, f_id, drift_corrected, fields):
        """
        This is designed for internal use. 

        This returns the data for a single fiducial.
        """
        index = 0
        fdcl = {}
        ml = self.getMovieLength()
        for fnum, locs in self.localizationsIterator(drift_corrected = drift_corrected, fields = fields):
                
            if not ("fiducial_id" in locs):
                raise FiducialException("File contains no fiducial information.")

            # Create fdcl dictionary if it is empty.
            if not fdcl:
                for elt in locs:
                    if not (elt == "fiducial_id"):
                        fdcl[elt] = numpy.zeros(ml, dtype = locs[elt].dtype)
                fdcl["frame"] = numpy.zeros(ml, dtype = numpy.int32)

            # Check for current fiducial in this frame.
            temp = numpy.where(locs["fiducial_id"] == f_id)[0]

            # There should only be 0 or 1 matches.
            assert (temp.size < 2)
                
            if (temp.size == 1):
                fd_index = temp[0]
                for elt in fdcl:
                    if (elt != "frame"):
                        fdcl[elt][index] = locs[elt][fd_index]
                fdcl["frame"][index] = fnum                            
                index += 1

        # Chop off extra data in arrays (if any).
        if fdcl:
            for elt in fdcl:
                fdcl[elt] = fdcl[elt][:index]

        return fdcl

    def getFiducialDataAll(self, drift_corrected, fields):
        """
        This is designed for internal use. 

        This returns the data for all of the fiducials.
        """

        # Create storage for fiducial dictionaries.
        indices =[]
        fdcls = []
        for i in range(self.getNFiducials()):
            indices.append(0)
            fdcls.append({})

        ml = self.getMovieLength()
        for fnum, locs in self.localizationsIterator(drift_corrected = drift_corrected, fields = fields):
                
            if not ("fiducial_id" in locs):
                raise FiducialException("File contains no fiducial information.")

            for i, fdcl in enumerate(fdcls):

                # Create fdcl dictionary if it is empty.
                if not fdcl:
                    for elt in locs:
                        if not (elt == "fiducial_id"):
                            fdcl[elt] = numpy.zeros(ml, dtype = locs[elt].dtype)
                    fdcl["frame"] = numpy.zeros(ml, dtype = numpy.int32)

                # Check for current fiducial in this frame.
                temp = numpy.where(locs["fiducial_id"] == i)[0]

                # There should only be 0 or 1 matches.
                assert (temp.size < 2)
                
                if (temp.size == 1):
                    index = indices[i]
                    fd_index = temp[0]
                    for elt in fdcl:
                        if (elt != "frame"):
                            fdcl[elt][index] = locs[elt][fd_index]
                        fdcl["frame"][index] = fnum                            
                    indices[i] += 1

        # Chop off extra data in arrays (if any).
        for i, fdcl in enumerate(fdcls):
            for elt in fdcl:
                fdcl[elt] = fdcl[elt][:indices[i]]

        return fdcls
        
    def getNFiducials(self):
        return (int(self.hdf5.attrs['n_fiducials']))

    def fiducialsIterator(self, fields = None, drift_corrected = False, preload_all = True):
        """
        An iterator for working with fiducials.

        for fdcl in h5.fiducialsIterator():
           ..

        fdcl is dictionary with all the fields for each localization in
        each frame that is a fiducial. Note that this never includes
        the "fiducial_id" field.

        fields - (Localization) fields (i.e. "x", "y", etc..), None = all of them.
        drift_corrected - Use drift corrected values for fiducial locations.
        preload_all - Load all the fiducials at initialization.
        
        Note: The default (preload_all = True) is much faster since it avoids multiple
              iterations across the HDF5 file, but it assumes that there is enough
              memory to store all the fiducial data.
        """
        # Add 'fiducial_id' as we'll also need this field.
        if fields is not None:
            fields.append("fiducial_id")

        if preload_all:
            fdcls = self.getFiducialDataAll(drift_corrected, fields)
            for i in range(len(fdcls)):
                yield fdcls[i]
        else:
            for i in range(self.getNFiducials()):
                yield self.getFiducialData(i, drift_corrected, fields)

    def hasFiducials(self):
        return "n_fiducials" in self.hdf5.attrs

    def setNFiducials(self, n_fiducials):
        self.hdf5.attrs['n_fiducials'] = n_fiducials
        

def trackFiducials(sa_hdf5_filename, max_gap = 0, radius = 0.0, reference_frame = 0):
    """
    max_gap - The maximum number of frames with no objects before a track is 
              considered to be terminated.
    radius - Maximum distance for an object to be in a track in pixels.
    reference_frame - The localizations in this frame will be used as the
                      fiducial reference locations.
    """

    with SAH5Fiducials(sa_hdf5_filename) as h5:
        ref_locs = h5.getLocalizationsInFrame(reference_frame, fields = ["x", "y"])

        if not ref_locs:
            raise FiducialException("No localizations in frame " + str(reference_frame))

        # Create fiducials list.
        fiducials = []
        for i in range(ref_locs["x"].size):
            fd = Fiducial(track_id = i)
            fd.addLocalization(ref_locs, i)
            fiducials.append(fd)

        h5.setNFiducials(len(fiducials))
        
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
                    h5.addFiducialID(locs_fd_id, fnum)

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
                        fiducials[i].addLocalization(locs, index_fd[i])
                        locs_fd_id[index_fd[i]] = fiducials[i].getFiducialId()

                h5.addFiducialID(locs_fd_id, fnum)
                                          
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
                if(elt.getLastAdded() <= max_gap):
                    fiducials.append(elt)
                else:
                    print("Lost fiducial", elt.getFiducialId(), "at frame", fnum)


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Fiducial tracking')

    parser.add_argument('--bin', dest='hdf5', type=str, required=True,
                        help = "The name of the localizations file.")
    parser.add_argument('--radius', dest='max_gap', type=float, required=True,
                        help = "Tracking radius in pixels.")
    parser.add_argument('--max_gap', dest='max_gap', type=int, required=False, default = 0,
                        help = "Maximum allowed gap in frames for fiducials. Default is 0.")
    parser.add_argument('--ref_frame', dest='ref_frame', type=int, required=False, default = 0,
                        help = "Frame to use for fiducial reference positions. Default is 0.")

    args = parser.parse_args()

    trackFiducials(args.hdf5, max_gap = args.max_gap, radius = args.radius, reference_frame = args.ref_frame)
