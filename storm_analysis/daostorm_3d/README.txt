
Python Programs:

batch_analysis - Used to analyze an entire directory of movies with the
   same parameters.
   
find_peaks - Configuration function for 3D-DAOSTORM peak finding.

mufit_analysis - Use this to perform peak finding and fitting using 3D-DAOSTORM.

z_calibration - Fit wx/wy versus z curves for use in 3D astigmatism imaging.


A sample run (execute this command in the 3d_daostorm/sample_data directory):

(1) .dax format data

$ python ../mufit_analysis.py --movie comp.dax --bin comp.hdf5 --xml 3d_zfit.xml 
Peak finding
Frame: 1 426 426
Frame: 2 423 849

Added 849
   0 fits lost due to Cholesky failure
   0 fits lost to image margin
   0 fits lost to negative value fit function
   2 fits lost to negative height
   0 fits lost to negative width

Tracking
Molecules: 849 (comp_mlist.bin)
Descriptor: 1
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 849 tracks
Analysis complete


(2) .tif format data

$ python ../mufit_analysis.py --movie comp.tif --bin comp.hdf5 --xml 3d_zfit.xml 
Peak finding
Frame: 1 426 426
Frame: 2 423 849

Added 849
   0 fits lost due to Cholesky failure
   0 fits lost to image margin
   0 fits lost to negative value fit function
   2 fits lost to negative height
   0 fits lost to negative width

Tracking
Molecules: 849 (comp_mlist.bin)
Descriptor: 1
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 849 tracks
Analysis complete


comp.dax is the STORM movie in .dax format. This is a raw 16 bit unsigned integer 
   format. Each .dax file has .inf file associated with it that provides some key details
   that are necessary for an application to be able to read the .dax file. The .inf file
   must have the same name as the .dax file (e.g. movie.dax <-> movie.inf). An example
   .dax file and it's associated .inf file can be found in the sample_data folder.

comp.hdf5 is a binary file that will contain the localizations returned
   by the analysis in HDF5 format. If the analysis has been run previously an
   attempt will be made to restart the analysis from where the last analysis
   finished/crashed. If you want "fresh" analyis you'll need to delete this file.

3d_zfit.xml describes how to do the analysis. Other example xml file can be
   found in storm_analysis/test/data directory.

If drift correction is performed, a file called comp_drift.txt will be created
containing the XYZ and offsets by frame of the localizations.
