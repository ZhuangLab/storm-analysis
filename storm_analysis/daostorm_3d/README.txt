
A sample run (execute this command in the 3d_daostorm/sample_data directory):

(1) .dax format data

$ python
>>> from storm_analysis.daostorm_3d.mufit_analysis import analyze
>>> analyze("comp.dax", "comp_mlist.bin", "3d_zfit.xml")
Peak finding
('Frame:', 0, 430, 430)
('Frame:', 1, 426, 856)

('Added', 856)

Tracking
Molecules: 856 (comp_mlist.bin)
Descriptor: 1
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 856 tracks
Analysis complete

(2) .tif format data

$ python
>>> from storm_analysis.daostorm_3d.mufit_analysis import analyze
>>> analyze("comp.tif", "comp_mlist.bin", "3d_zfit.xml")
Peak finding
('Frame:', 0, 430, 430)
('Frame:', 1, 426, 856)

('Added', 856)

Tracking
Molecules: 856 (comp_mlist.bin)
Descriptor: 1
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 856 tracks
Analysis complete


comp.dax is the STORM movie in .dax format. This is a raw 16 bit unsigned integer 
format. Each .dax file has .inf file associated with it that provides some key details
that are necessary for an application to be able to read the .dax file. The .inf file
must have the same name as the .dax file (e.g. movie.dax <-> movie.inf). An example
.dax file and it's associated .inf file can be found in the sample_data folder.

comp_mlist.bin is a binary file that will contain the localizations returned
   by the analysis in Insight3 format. Note that the _mlist.bin extension is
   an important part of the name and it is best not to substitute this for
   something else. Also, if the analysis has been run previously an attempt
   will be made to restart the analysis from where the last analysis finished/
   crashed. If you want "fresh" analyis you'll need to delete this file.

3d_zfit.xml describes how to do the analysis. Other example xml file can be
   found in storm_analysis/test/data directory.


In addition, when the tracking radius is greater than zero a file called
comp_alist.bin will be created. In this file, localizations in subsequent
frames that are within radius of a localization(s) in previous frame will be
averaged together to create a single localization.

If drift correction is performed, a file called comp_drift.txt will be created
containing the XYZ and offsets by frame of the localizations. This correction
is automatically applied to _alist.bin file if it is created, otherwise it
is only applied to _mlist.bin file.

The file contains the molecule localizations in Insight3 format. A
seperate program for visualizing the localizations is available by
request from the Zhuang lab. Alternatively, you can use the
bin_to_tagged_spot_file.py program in the sa_utilities directory
to convert the .bin file to a tagged spot file format (.tsf) file.
These files can be visualized using Micro-Manager.
