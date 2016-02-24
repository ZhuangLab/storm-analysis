#!/bin/bash

gcc -fPIC -g -c -Wall -O3 scmos_utilities.c
gcc -shared -Wl,-soname,scmos_utilities.so.1 -o scmos_utilities.so.1.0.1 scmos_utilities.o
ln -s scmos_utilities.so.1.0.1 scmos_utilities.so

