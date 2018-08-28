#!/usr/bin/env python
"""
A sub-class of SA5Py for working with clusters.

Hazen 08/18
"""
import numpy

import storm_analysis.sa_library.sa_h5py as saH5Py


class SAH5Clusters(saH5Py.SAH5Py):
    """
    A sub-class of SAH5Py designed for use in clustering.
    """
    def addClusters(self, cluster_id, cluster_data):
        """
        Add clustering information to the H5 file.

        cluster_id is a numpy array containing the cluster number
                   assigned to each localization.

        cluster_data is a dictionary of numpy arrays that includes
                     localization information such as track group
                     and track element as well as other data that
                     should be associated with each cluster. Each
                     of the numpy arrays in cluster_data must be
                     the same length as the cluster_id array.

        It is assumed that negative values in cluster_id are
        used for the localizations that were not assigned to a
        cluster.
        """
        # Delete off clustering information, if any.
        if self.hasClusters():
            del self.hdf5["clusters"]

        clusters_grp = self.hdf5.create_group("clusters")

        # Check that arrays are the correct size.
        for data in cluster_data:
            assert (cluster_data[data].size == cluster_id.size), "Incorrect size for data " + data

        # Add clustering information.
        start = numpy.amin(cluster_id)
        end = numpy.amax(cluster_id) + 1
        n_clusters = 0

        # Add empty cluster zero if all the localizations / tracks
        # were assigned to a cluster.
        if (start >= 0):
            cl_grp = clusters_grp.create_group(self.getClusterName(n_clusters))
            cl_grp.attrs['cl_size'] = 0
            n_clusters += 1

        for i in range(start, end):
            cl_mask = (cluster_id == i)
            cl_size = numpy.count_nonzero(cl_mask)
            if (cl_size > 0):
                cl_grp = clusters_grp.create_group(self.getClusterName(n_clusters))
                for field in cluster_data:
                    cl_grp.create_dataset(field, data = cluster_data[field][cl_mask])
                cl_grp.attrs['cl_size'] = cl_size
                
                n_clusters += 1

        # The first 'cluster' is not actually a cluster, it contains all the
        # track / localizations that were not assigned to a cluster.
        clusters_grp.attrs['n_clusters'] = n_clusters - 1

        # Record fields that were saved with the cluster.
        clusters_grp.attrs['fields'] = ",".join(list(cluster_data.keys()))
    
    def clustersIterator(self, fields = None, min_size = 1, skip_unclustered = True):
        """
        An iterator for getting all the clusters in a for loop.

        for cl in h5.clustersIterator():
            ..

        'cl' will be a dictionary containing all the fields, or
        all the requested fields for each track / localization in
        the cluster.

        You can set the fields parameter to a list if you only
        want some of the fields, fields = ['x', 'y', 'z'] for 
        example.

        You can use the min_size parameter to only return those
        clusters above a certain size. This should be at least
        as large as minimum cluster size used when clustering.

        Cluster 0 contains all the tracks / localizations that
        were not assigned to a cluster. The default behavior is
        to not return this cluster.
        """
        if (not self.hasClusters()):
            return
                
        else:
            start = 1
            if not skip_unclustered:
                start = 0
                
            for i in range(start, self.getNClusters() + 1):
                cl_grp = self.getClusterGroup(i)
                if (cl_grp.attrs['cl_size'] >= min_size):
                    yield [i, self.getClusterData(i, fields = fields)]

    def getCluster(self, group):
        """
        Get the data in a single cluster. The dictionary this returns
        only includes information about how to look up the tracks or 
        localizations in the group.
        """
        cl_dict = {}
        for field in group:
            cl_dict[field] = group[field][()]
        return cl_dict

    def getClusterData(self, index, fields = None):
        """
        Return a dictionary containing all the fields, or all 
        the requested fields for each track / localization in 
        the cluster.

        Notes: 
        1. The recommended approach is to use clustersIterator().
        2. Clustering is always done on drift corrected data.
        3. This is not very efficient as it loads all the data
           for each track group / frame for every cluster element,
           then only uses a single element of the data.
        """
        cl = {}
        
        # Get the cluster data. These are the arrays that were
        # saved with the cluster. They include information that
        # is used to find the localization / track.
        #
        cl_dict = self.getCluster(self.getClusterGroup(index))

        if not(cl_dict):
            return cl

        # Check if we have all the data in the cluster so we can
        # shortcut having to pull data from the tracks / localizations.
        #
        # Note: If the user does not specify particular fields this
        #       is automatically false.
        #
        if fields is not None:
            cl_fields = self.getClustersFields()

            have_all_fields = True
            for elt in fields:
                if not elt in cl_fields:
                    have_all_fields = False
                    break

            if have_all_fields:
                for elt in fields:
                    cl[elt] = cl_dict[elt]
                return cl

        # Check for fields that localizations / tracks won't have.
        cl_fields = None
        sa_fields = None
        if fields is not None:
            cl_fields = []
            sa_fields = []
            for elt in fields:

                # These are available in the cluster.
                #
                # Note that if a field is also available with the
                # localization/track we'll still get the cluster's
                # version. The assumption is that they are the same.
                #
                if elt in self.getClustersFields():
                    cl_fields.append(elt)

                # Everything else is assumed to available in the
                # localizations / tracks.
                else:
                    sa_fields.append(elt)

        # Is the data for localizations or tracks?
        if "frame" in cl_dict:
            cl_size = cl_dict["frame"].size
            for i in range(cl_size):
                locs = self.getLocalizationsInFrame(int(cl_dict["frame"][i]),
                                                    drift_corrected = True,
                                                    fields = sa_fields)
                if not cl:
                    for field in locs:
                        cl[field] = numpy.zeros(cl_size, locs[field].dtype)
                        
                for field in locs:
                    cl[field][i] = locs[field][cl_dict["loc_id"][i]]

        else:
            cl_size = cl_dict["track_id"].size
            for i in range(cl_size):
                tracks = self.getTracksByIndex(int(cl_dict["track_id"][i]),
                                               fields = sa_fields)
                if not cl:
                    for field in tracks:
                        cl[field] = numpy.zeros(cl_size, tracks[field].dtype)
                        
                for field in tracks:
                    cl[field][i] = tracks[field][cl_dict["loc_id"][i]]

        # Also add the cluster data arrays.
        if cl_fields is None:
            for field in cl_dict:

                # Only add if we don't already have data for this field.
                if not field in cl:
                    cl[field] = cl_dict[field]
                    
        else:
            for field in cl_fields:
                cl[field] = cl_dict[field]
           
        return cl
    
    def getClusterGroup(self, index):
        """
        Gets a cluster group, this is group contains all the data
        in a single cluster.
        """
        cl_name = self.getClusterName(index)
        if cl_name in self.getClusters():
            return self.getClusters()[cl_name]

    def getClusteringInfo(self):
        """
        Returns the (short) string describing how the clustering was done.
        """
        return self.getClusters().attrs["info"]

    def getClusterName(self, index):
        return "cl_" + str(index)
        
    def getClusters(self):
        return self.hdf5["clusters"]

    def getClustersFields(self):
        return self.getClusters().attrs["fields"].split(",")

    def getDataForClustering(self):
        """
        This return the X/Y/Z locations of all the tracks or localizations
        in a clustering friendly format.

        Returns [x, y, z, c, cluster_data] where x, y are in pixels, z is
                in microns and c is the localization category.
        """
        cluster_data = {}

        # Load tracks.
        if self.hasTracks():
            total_tracks = self.getNTracks()

            loc_id = numpy.zeros(total_tracks, dtype = numpy.int)
            track_id = numpy.zeros(total_tracks, dtype = numpy.int)
            x = numpy.zeros(total_tracks)
            y = numpy.zeros(total_tracks)
            z = numpy.zeros(total_tracks)
            c = numpy.zeros(total_tracks, dtype = numpy.int32)
            
            start = 0
            for i, tracks in enumerate(self.tracksIterator(fields = ['x', 'y', 'z', 'category'])):
                n_tracks = tracks['x'].size
                end = start + n_tracks
                
                loc_id[start:end] = numpy.arange(n_tracks)
                track_id[start:end] = i
                x[start:end] = tracks['x']
                y[start:end] = tracks['y']
                z[start:end] = tracks['z']
                c[start:end] = tracks['category']

                start += n_tracks

            cluster_data['loc_id'] = loc_id
            cluster_data['track_id'] = track_id

        # Load localizations:
        else:
            total_locs = self.getNLocalizations()

            frame = numpy.zeros(total_locs, dtype = numpy.int)            
            loc_id = numpy.zeros(total_locs, dtype = numpy.int)
            x = numpy.zeros(total_locs)
            y = numpy.zeros(total_locs)
            z = numpy.zeros(total_locs)
            c = numpy.zeros(total_locs, dtype = numpy.int32)

            fields = ['x', 'y', 'category']
            
            # Check if the localizations have the 'z' field.
            for f_num, locs in self.localizationsIterator(fields = fields):
                if "z" in locs:
                    fields.append('z')
                break

            start = 0
            for f_num, locs in self.localizationsIterator(fields = fields):
                n_locs = locs['x'].size
                end = start + n_locs

                frame[start:end] = f_num
                loc_id[start:end] = numpy.arange(n_locs)
                x[start:end] = locs['x']
                y[start:end] = locs['y']
                if 'z' in fields:
                    z[start:end] = locs['z']
                c[start:end] = locs['category']

                start += n_locs

            cluster_data['frame'] = frame
            cluster_data['loc_id'] = loc_id
            
        return [x, y, z, c, cluster_data]
                             
    def getNClusters(self):
        """
        Return the number of clusters.
        """
        return self.getClusters().attrs['n_clusters']
        
    def hasClusters(self):
        return ("clusters" in self.hdf5)

    def setClusteringInfo(self, info):
        """
        info is a (short) string describing how the clustering was done.
        """
        self.getClusters().attrs["info"] = info
        
