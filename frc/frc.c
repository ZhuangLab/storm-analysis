/*
 *
 * C utilities for calculating FRC following Nieuwenhuizen, 
 * Nature Methods, 2013.
 *
 * Hazen 10/14
 *
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall frc.c
 *  gcc -shared -Wl,-soname,frc.so.1 -o frc.so.1.0.1 frc.o
 *  (and you need to create the symlink frc.so frc.so.1.0.1)
 * 
 * Windows (64bit, MinGW):
 *  c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -c frc.c -O3
 *  c:\MinGW64\bin\x86_64-w64-mingw32-gcc.exe -shared -o frc.dll frc.o
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Function Declarations */
void calc_frc(double *, double *, double *, int *, int, int);

/*
 * calc_frc()
 *
 * Calculate the FRC given the image FFTs and somewhere
 * to store the results. Note that the FFTs are actually 
 * complex numbers with the real part first, followed by 
 * the imaginary part. Note also that it is assumed that 
 * the zero frequency is located in the center of the 
 * array, not in the corners.
 *
 * fft1 - Fourier transform of the first image (complex).
 * fft2 - Fourier transform of the second image (complex).
 * frc - Storage for the FRC results (double).
 * frc_counts - Number of values in each bin of frc.
 * y_size - Size of the FFT array in y (arr.shape[0]).
 * x_size - Size of the FFT array in x (arr.shape[1]).
 */
void calc_frc(double *fft1, double *fft2, double *frc, int *frc_counts, int y_size, int x_size)
{
  int cx,cy,dx,dy,i,j,q;
  int real,ima;
  double g1,g2,g1g2;

  cx = x_size/2;
  cy = y_size/2;

  for(i=0;i<y_size;i++){
    dy = i - cy;
    for(j=0;j<x_size;j++){
      dx = j - cx;
      real = i*2*x_size+2*j;
      ima = i*2*x_size+2*j+1;
      q = (int)(sqrt(dx*dx+dy*dy)+0.5);
      g1 = fft1[real]*fft1[real] + fft1[ima]*fft1[ima];
      g2 = fft2[real]*fft2[real] + fft2[ima]*fft2[ima];
      g1g2 = fft1[real]*fft2[real] + fft1[ima]*fft2[ima];

      /*
	printf("%d %d %f %f\n", i, j, fft1[real], fft1[ima]);
      */

      frc[q] += g1g2/sqrt(g1*g2);
      frc_counts[q] += 1;
    }
  }
}

/*
 * The MIT License
 *
 * Copyright (c) 2014 Zhuang Lab, Harvard University
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
