
The RCC code has the following dependencies:

Python
Python - numpy library - http://numpy.scipy.org/
Python - scipy library - http://numpy.scipy.org/
Python - matplotlib - http://matplotlib.org/
(64 bit windows installers for these libraries are available here http://www.lfd.uci.edu/~gohlke/pythonlibs/)

You will need to add the project root directory (i.e. the directory one
up from the location of this file) to your Python path.

Files:
 create_drift_data.py - Creates .bin files for testing drift correction.

 rcc-drift-correction.py - The program for calculating drift during STORM
   movies using the RCC approach.

Usage:
 > python rcc-drift-correction.py movie_mlist.bin movie_drift.txt 500 2

The file movie_mlist.bin is the input molecule list file. The results are
saved in the movie_drift.txt file which specifies and the calculated drift
in x,y and z for every frame. 500 is the number of frames to bin into a
single STORM image for the purposes of calculating correlations between
frames. 2 is the up-scaling factor to use for the STORM image. For example,
if the original data is 256x256 then the STORM image will be 512x512.

Also, if you add an additional parameter on the end then the z drift
calculation will not be performed:

 > python rcc-drift-correction.py movie_mlist.bin movie_drift.txt 500 2 1


A sample run:
(execute this command in the rcc directory)
python rcc-drift-correction.py test.bin test_drift.txt 50 1 1

If this works correctly you will see the following output:

Version: M425
Frames: 1
Status: 6
Molecules: 50000

Could not find movie file for test.bin assuming 256x256x300
Performing XY correction.
offset between frame ranges  0 - 50  and  50 - 100
 ->  0.013969514764 0.00769523363573 good

offset between frame ranges  0 - 50  and  100 - 150
 ->  0.00760025915068 0.000176801037412 good

offset between frame ranges  0 - 50  and  150 - 200
 ->  0.00366980012723 0.00904294388613 good

offset between frame ranges  0 - 50  and  200 - 250
 ->  0.00411372203374 -0.00398676872976 good

offset between frame ranges  0 - 50  and  250 - 300
 ->  0.0127954954698 -0.00165640268452 good

offset between frame ranges  50 - 100  and  100 - 150
 ->  25.9297811985 -13.5312941363 good

offset between frame ranges  50 - 100  and  150 - 200
 ->  -0.0207305933497 -0.00263870125823 good

offset between frame ranges  50 - 100  and  200 - 250
 ->  -0.00438144438073 -0.00733064476375 good

offset between frame ranges  50 - 100  and  250 - 300
 ->  -0.00419918100417 -0.00551800740828 good

offset between frame ranges  100 - 150  and  150 - 200
 ->  0.00555719703408 0.0129130288662 good

offset between frame ranges  100 - 150  and  200 - 250
 ->  -0.0107275414252 -0.00794143199144 good

offset between frame ranges  100 - 150  and  250 - 300
 ->  0.00871247592679 -0.00517362171666 good

offset between frame ranges  150 - 200  and  200 - 250
 ->  -3.13044346854e-05 -0.013009016387 good

offset between frame ranges  150 - 200  and  250 - 300
 ->  0.00923952455031 -0.0109188719213 good

offset between frame ranges  200 - 250  and  250 - 300
 ->  0.0102302285098 0.0025337630199 good

--
14 removing 5 with error 19.5028
13 removing 5 with error 4.88084
not removing 7 with error 4.8806
not removing 9 with error 4.87833
not removing 6 with error 4.87827
not removing 1 with error 4.87395
not removing 0 with error 4.87354
not removing 5 with error 4.87011
not removing 8 with error 4.86987
0 25 0.0 0.0
1 75 0.0130622312427 0.0050151203759
2 125 0.00662893755361 0.00132710300386
3 175 0.00582106690854 0.0104281092063
4 225 0.00348199810833 -0.00392163638026
5 275 0.0131545541808 -0.00157688884065
