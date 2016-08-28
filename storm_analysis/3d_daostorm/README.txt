
The 3D-DAOSTORM code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
Python - PIL image library - http://www.pythonware.com/products/pil/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)

LAPACK - http://www.netlib.org/lapack/


In order to use the 3D-DAOSTORM code you will first need to compile several
C programs and libraries:

Libraries (These are all in the sa_library folder):
 grid - C functions for gridding 2D and 3D data.

 ia_utilities - A collection of C image analysis utility functions used by 3D-DAOSTORM.

 multi_fit - The core C localization fitting routines.


Programs:
 (These programs can all be found in the sa_utilities folder)
 avemlist - This program averages all the objects in a track into a single object.

 apply-drift-correction - Applies a previously determined drift correction to each
    of the localizations.

 fitz - This program is used to determine the z position from the localization x and y
    widths and a previously determined calibration curve.

 tracker - This program tracks objects across multiple frames & assigns the
    appropriate category to each object (i.e. specific or non-specific activation, etc.)


Linux compilation examples:
gcc -fPIC -g -c -Wall ia_utilities.c
gcc -shared -Wl,-soname,ia_utilities.so.1 -o ia_utilities.so.1.0.1 ia_utilities.o -lc

gcc -fPIC -g -c -Wall multi_fit.c
gcc -shared -Wl,-soname,multi_fit.so.1 -o multi_fit.so.1.0.1 multi_fit.o -lc -llapack

gcc tracker.c -o tracker -lm
gcc avemlist.c -o avemlist -lm
gcc fitz.c -o fitz -lm
gcc apply-drift-correction.c -o apply-drift-correction

Once these are compiled you will need to create symlinks in the same directories:
ln -s ia_utilities.so.1.0.1 ia_utilities.so
ln -s multi_fit.so.1.0.1 multi_fit.so


Windows compilation examples (using a 64 bit MinGW compiler):
C:\MinGW64\bin\x86_64-w64-mingw32-gcc -c ia_utilities.c
C:\MinGW64\bin\x86_64-w64-mingw32-gcc -shared -o ia_utilities.dll ia_utilities.o

C:\MinGW64\bin\x86_64-w64-mingw32-gcc -c multi_fit.c
C:\MinGW64\bin\x86_64-w64-mingw32-gcc -shared -o multi_fit.dll multi_fit.o -llapack -Lc:\path\to\lapack

C:\MinGW64\bin\x86_64-w64-mingw32-gcc tracker.c -o tracker
C:\MinGW64\bin\x86_64-w64-mingw32-gcc avemlist.c -o avemlist
C:\MinGW64\bin\x86_64-w64-mingw32-gcc fitz.c -o fitz
C:\MinGW64\bin\x86_64-w64-mingw32-gcc apply-drift-correction.c -o apply-drift-correction


Once you have compiled all the C libraries you will need to add the project root 
directory (i.e. the directory one up from the location of this file) to your Python 
path (Also, at least on windows, the llapack library needs to in a directory on PATH). 
Then you should be able to go ahead and run the 3D-DAOSTORM analysis. You can run it 
on a .dax format STORM movie using the following command:

python /path/to/mufit_analysis.py movie.dax movie_mlist.bin analysis_params.xml

Or on .tif format STORM movie with this command:
python /path/to/mufit_analysis.py movie.tif movie_mlist.bin analysis_params.xml

movie.dax is the STORM movie in .dax format. This is a raw 16 bit unsigned integer 
format. Each .dax file has .inf file associated with it that provides some key details
that are necessary for an application to be able to read the .dax file. The .inf file
must have the same name as the .dax file (e.g. movie.dax <-> movie.inf). An example
.dax file and it's associated .inf file can be found in the sample_data folder.

movie_mlist.bin is a binary file that will contain the localizations returned
   by the analysis in Insight3 format. Note that the _mlist.bin extension is
   an important part of the name and it is best not to substitute this for
   something else. Also, if the analysis has been run previously an attempt
   will be made to restart the analysis from where the last analysis finished/
   crashed. If you want "fresh" analyis you'll need to delete this file.

analysis_params.xml describes how to do the analysis. See the example.xml for
   an example.

In addition, when the tracking radius is greater than zero a file called
movie_alist.bin will be created. In this file, localizations in subsequent
frames that are within radius of a localization(s) in previous frame will be
averaged together to create a single localization.

If drift correction is performed, a file called movie_drift.txt will be created
containing the XYZ and offsets by frame of the localizations. This correction
is automatically applied to _alist.bin file if it is created, otherwise it
is applied to _mlist.bin file.


A sample run: (execute this command in the 3d_daostorm/sample_data directory)
Analyze the .dax format sample data:
python ../mufit_analysis.py comp.dax comp_mlist.bin 3d_zfit.xml

Or analyze the .tif format sample data:
python ../mufit_analysis.py comp.tif comp_mlist.bin 3d_zfit.xml

If the program executes correctly you should see the following output:

Peak finding
Frame: 0 430 430
Frame: 1 426 856

Added 856

Tracking
Molecules: 856 (comp_mlist.bin)
Descriptor: 1
Processing molecule 0 in frame 0 (tracker)
Finished processing
Found 856 tracks
Analysis complete


It will create the following files:
comp_mlist.bin

The comp_alist.bin file is not created as tracking radius was set to zero
in the 3d_zfit.xml file.

This file contains the molecule localizations in Insight3 format. A
seperate program for visualizing the localizations is available by
request from the Zhuang lab. Alternatively, you can use the
bin_to_tagged_spot_file.py program in the sa_utilities directory
to convert the .bin file to a tagged spot file format (.tsf) file.
These files can be visualized using Micro-Manager.
