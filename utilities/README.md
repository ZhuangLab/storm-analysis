
Python Programs:




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

