#!/bin/bash

gcc -fPIC -g -c -Wall -O3 kdtree.c
gcc -shared -Wl,-soname,kdtree.so.1 -o kdtree.so.1.0.1 kdtree.o
ln -s kdtree.so.1.0.1 kdtree.so

gcc -fPIC -g -c -Wall -O3 dbscan.c
gcc -shared -Wl,-soname,dbscan.so.1 -o dbscan.so.1.0.1 dbscan.o kdtree.o
ln -s dbscan.so.1.0.1 dbscan.so
