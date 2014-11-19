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

cd ../L1H
gcc -fPIC -g -c -Wall -O3 homotopy_common.c
gcc -fPIC -g -c -Wall -O3 homotopy_general.c
gcc -fPIC -g -c -Wall -O3 homotopy_storm.c

#
# FIXME: Do I have the right flag here for SSE2?
#    Does this happen automatically anyway for gcc on linux?
#
gcc -fPIC -g -c -Wall -O3 -msse2 homotopy_sse.c
gcc -fPIC -g -c -Wall -O3 homotopy_imagea_common.c
gcc -fPIC -g -c -Wall -O3 homotopy_imagea.c

gcc -shared -Wl,-soname,homotopy_general.so.1 -o homotopy_general.so.1.0.1 homotopy_general.o homotopy_common.o -llapack
ln -s homotopy_general.so.1.0.1 homotopy_general.so
gcc -shared -Wl,-soname,homotopy_storm.so.1 -o homotopy_storm.so.1.0.1 homotopy_storm.o homotopy_common.o -llapack
ln -s homotopy_storm.so.1.0.1 homotopy_storm.so
gcc -shared -Wl,-soname,homotopy_sse.so.1 -o homotopy_sse.so.1.0.1 homotopy_sse.o homotopy_common.o -llapack
ln -s homotopy_sse.so.1.0.1 homotopy_sse.so

gcc -shared -Wl,-soname,homotopy_ia_storm.so.1 -o homotopy_ia_storm.so.1.0.1 homotopy_imagea.o homotopy_storm.o homotopy_imagea_common.o homotopy_common.o -llapack
ln -s homotopy_ia_storm.so.1.0.1 homotopy_ia_storm.so
gcc -shared -Wl,-soname,homotopy_ia_sse.so.1 -o homotopy_ia_sse.so.1.0.1 homotopy_imagea.o homotopy_sse.o homotopy_imagea_common.o homotopy_common.o -llapack
ln -s homotopy_ia_sse.so.1.0.1 homotopy_ia_sse.so

