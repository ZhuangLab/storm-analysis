/*
 * C library for sCMOS image processing. This is based on the algorithms described in:
 *
 * "Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms"
 * F. Huang et al. Nature Methods, 10, p653-658.
 *
 * Hazen 10/13
 *
 * Compilation instructions:
 *
 * Windows:
 *  gcc -c smooth.c -O3
 *  gcc -shared -o smooth.dll smooth.o
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Function Declarations */
void deregularize(double *, double *, double *, double *, double *, int);
void regularize(double *, double *, double *, double *, double *, int);
void smooth(double *, double *, double *, double *, double *, int, int, int, int);

/*
 * deregularize()
 *
 * Given the initial image & camera offset, variance and gain parameters
 * this "deregularizes" the image, basically it is the inverse of regularize().
 *
 * output - The "deregularized" image.
 * image - The regularized image.
 * offset - The camera offset values.
 * variance - The camera variance values.
 * gain - The camera gain values.
 * size - The total number of pixels in the image.
 */
void deregularize(double *output, double *image, double *offset, double *variance, double *gain, int size)
{
  int i;
  double inv_gain;

  for(i=0;i<size;i++){
    inv_gain = 1.0/gain[i];
    output[i] = (image[i] - variance[i]*inv_gain*inv_gain)*gain[i] + offset[i];
  }
}

/*
 * regularize()
 *
 * Given the initial image & camera offset, variance and gain parameters
 * this "regularizes" (if that is the right word) the image. See
 * section 3.2 in the paper supplement.
 *
 * output - The "regularized" image.
 * image - The original image.
 * offset - The camera offset values.
 * variance - The camera variance values.
 * gain - The camera gain values.
 * size - The total number of pixels in the image.
 */
void regularize(double *output, double *image, double *offset, double *variance, double *gain, int size)
{
  int i;
  double inv_gain;

  for(i=0;i<size;i++){
    inv_gain = 1.0/gain[i];
    output[i] = (image[i] - offset[i])*inv_gain + variance[i]*inv_gain*inv_gain;
  }
}

/*
 * smooth()
 *
 * Given the initial image & camera offset, variance and gain parameters
 * this returns a smoothed image for peak finding. See section 3.1 in the
 * paper supplement.
 *
 * smoothed - The smoothed image (output).
 * image - The original image.
 * offset - Camera offset values.
 * inv_variance - 1.0/Camera variance values.
 * inv_gain - Camera gain values.
 * x_size - Size of the image in x.
 * y_size - Size of the image in y.
 * hw1 - Half width of the smaller convolution kernel.
 * hw2 - Half width of the larger convolution kernel.
 */
void smooth(double *smoothed, double *image, double *offset, double *inv_variance, double *inv_gain, int x_size, int y_size, int hw1, int hw2)
{
  int i,j,k,l,m,n,o;
  double t1, t2;
  double *work1;

  /* Calculate (image - offset)/(gain * variance) */
  work1 = (double *)malloc(sizeof(double)*x_size*y_size);

  for(i=0;i<x_size;i++){
    k = i*y_size;
    for(j=0;j<y_size;j++){
      l = k+j;
      work1[l] = (image[l] - offset[l])*inv_variance[l]*inv_gain[l];
    }
  }

  /* Smooth image */
  for(i=hw2;i<(x_size-hw2);i++){
    for(j=hw2;j<(y_size-hw2);j++){
      k = i*y_size+j;

      /* Calculate smaller convolution kernel. */
      t1 = 0.0;
      t2 = 0.0; 
      for(l=-hw1;l<=hw1;l++){
	n = k + l*y_size;
	for(m=-hw1;m<=hw1;m++){
	  o = n + m;
	  t1 += work1[o];
	  t2 += inv_variance[o];
	}
      }
      smoothed[k] = t1/t2;

      /* Calculate larger convolution kernel. */
      t1 = 0.0;
      t2 = 0.0; 
      for(l=-hw2;l<=hw2;l++){
	n = k + l*y_size;
	for(m=-hw2;m<=hw2;m++){
	  o = n + m;
	  t1 += work1[o];
	  t2 += inv_variance[o];
	}
      }
      
      smoothed[k] -= t1/t2;
    }
  }

  free(work1);
}

/*
 * The MIT License
 *
 * Copyright (c) 2013 Zhuang Lab, Harvard University
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
