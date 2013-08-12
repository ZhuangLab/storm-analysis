
Python Programs:

bin_to_tagged_spot_file.py - Convert a .bin file to the Micro-Manager tagged spot file
   (tsf) format. Tagged spot file format files can then be rendered with programs
   such as Micro-Manager and ViSP (1).

read_tagged_spot_file.py - Read .tsf format file. This is useful mostly as a debugging
   aid to make sure that the .tsf file gotten written properly (1).

track_average_correct.py - This does tracking, averaging and drift correction on a 
   molecule list file. It can be useful if you want to redo these parts without
   redoing the localization step, which takes most of the time.

xyz-drift-correction.py - This will determine (but not apply) the drift correction
   for a STORM movie.


C Programs:

tracker - This program tracks objects across multiple frames & assigns the
   appropriate category to each object (i.e. specific or non-specific activation, etc.)

avemlist - This program averages all the objects in a track into a single object.

fitz - This program is used to determine the z position from the localization x and y
   widths and a previously determined calibration curve.

apply-drift-correction - Applies a previously determined drift correction to each
   of the localizations.


Linux compilation examples:
gcc tracker.c -o tracker -lm
gcc avemlist.c -o avemlist -lm
gcc fitz.c -o fitz -lm
gcc apply-drift-correction.c -o apply-drift-correction


Windows compilation examples (using a 64 bit MinGW compiler):
C:\MinGW64\bin\x86_64-w64-mingw32-gcc tracker.c -o tracker
C:\MinGW64\bin\x86_64-w64-mingw32-gcc avemlist.c -o avemlist
C:\MinGW64\bin\x86_64-w64-mingw32-gcc fitz.c -o fitz
C:\MinGW64\bin\x86_64-w64-mingw32-gcc apply-drift-correction.c -o apply-drift-correction


(1) These programs require Google protocol buffers. I have also included the 
   (generated) TSFProto_pb2.py file that describes the tsf format. However, is the
   tsf format gets changed in Micro-Manager this may fall out of sync.

