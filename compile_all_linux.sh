#!/bin/bash

cd sa_library
gcc -fPIC -g -c -Wall -O3 multi_fit.c
gcc -shared -Wl,-soname,multi_fit.so.1 -o multi_fit.so.1.0.1 multi_fit.o -llapack
ln -s multi_fit.so.1.0.1 multi_fit.so

gcc -fPIC -g -c -Wall -O3 ia_utilities.c
gcc -shared -Wl,-soname,ia_utilities.so.1 -o ia_utilities.so.1.0.1 ia_utilities.o
ln -s ia_utilities.so.1.0.1 ia_utilities.so

gcc -fPIC -g -c -Wall -O3 grid.c
gcc -shared -Wl,-soname,grid.so.1 -o grid.so.1.0.1 grid.o
ln -s grid.so.1.0.1 grid.so

cd ../sa_utilities
gcc tracker.c -o tracker -lm
gcc apply-drift-correction.c -o apply-drift-correction
gcc fitz.c -o fitz -lm
gcc avemlist.c -o avemlist -lm

cd ../frc
gcc -fPIC -g -c -Wall frc.c -O3
gcc -shared -Wl,-soname,frc.so.1 -o frc.so.1.0.1 frc.o
ln -s frc.so.1.0.1 frc.so

cd ../sCMOS
gcc -fPIC -g -c -Wall -O3 scmos_utilities.c
gcc -shared -Wl,-soname,scmos_utilities.so.1 -o scmos_utilities.so.1.0.1 scmos_utilities.o
ln -s scmos_utilities.so.1.0.1 scmos_utilities.so


#./L1H/homotopy_common.c
#./L1H/homotopy_imagea.c
#./L1H/fista_lib.c
#./L1H/homotopy_imagea_common.c
#./L1H/homotopy_gpu.c
#./L1H/homotopy_general.c
#./L1H/homotopy_sse.c
#./L1H/homotopy_storm.c
