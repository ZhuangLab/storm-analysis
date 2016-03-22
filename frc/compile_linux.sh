#!/bin/bash

gcc -fPIC -g -c -Wall frc.c -O3
gcc -shared -Wl,-soname,frc.so.1 -o frc.so.1.0.1 frc.o
ln -s frc.so.1.0.1 frc.so
