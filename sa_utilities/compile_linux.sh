#!/bin/bash

gcc tracker.c -o tracker -lm
gcc apply-drift-correction.c -o apply-drift-correction
gcc fitz.c -o fitz -lm
gcc avemlist.c -o avemlist -lm

