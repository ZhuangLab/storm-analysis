#!/bin/bash

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

