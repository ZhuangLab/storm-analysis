#!/bin/bash

gcc -fPIC -g -c -Wall -O3 tracker.c
gcc -shared -Wl,-soname,tracker.so.1 -o tracker.so.1.0.1 tracker.o
ln -s tracker.so.1.0.1 tracker.so

gcc -fPIC -g -c -Wall -O3 apply-drift-correction.c
gcc -shared -Wl,-soname,apply-drift-correction.so.1 -o apply-drift-correction.so.1.0.1 apply-drift-correction.o
ln -s apply-drift-correction.so.1.0.1 apply-drift-correction.so

gcc -fPIC -g -c -Wall -O3 fitz.c
gcc -shared -Wl,-soname,fitz.so.1 -o fitz.so.1.0.1 fitz.o
ln -s fitz.so.1.0.1 fitz.so

gcc -fPIC -g -c -Wall -O3 avemlist.c
gcc -shared -Wl,-soname,avemlist.so.1 -o avemlist.so.1.0.1 avemlist.o
ln -s avemlist.so.1.0.1 avemlist.so

