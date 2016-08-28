#!/bin/bash

gcc tracker.c -o tracker -lc -lm
gcc apply-drift-correction.c -o apply-drift-correction -lc
gcc fitz.c -o fitz -lc -lm
gcc avemlist.c -o avemlist -lc -lm
