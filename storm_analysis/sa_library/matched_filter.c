/*
 * C library for the matched filter approach. Basically 
 * this is just convolving the image with the PSF.
 *
 * This uses an FFT to do the convolution so there will
 * be edge effects.
 *
 * Hazen 3/16
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall matched_filter.c
 *  gcc -shared -Wl,-soname,matched_filter.so.1 -o matched_filter.so.1.0.1 matched_filter.o -lc -lfftw3
 *  ln -s matched_filter.so.1.0.1 matched_filter.so
 *
 * Windows:
 *  gcc -c -O3 matched_filter.c
 *  gcc -shared -o matched_filter.dll matched_filter.o -lfftw3-3 c:\path\to\libfftw3-3.dll
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>

/* Structures & Types */
struct filter_struct {
  int fft_size;
  int image_size;
  int x_size;
  int y_size;
  double normalization;

  double *fft_vector;
  
  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *fft_vector_fft;
  fftw_complex *psf_fft;
};
typedef struct filter_struct filter;

/* Function Declarations */
void cleanup(filter *);
void convolve(filter *, double *, double *);
filter *initialize(double *, int, int, int);

/* Functions */

/*
 * cleanup()
 *
 * filter - A pointer to a filter structure.
 */
void cleanup(filter *flt)
{
  free(flt->fft_vector);
  
  fftw_destroy_plan(flt->fft_backward);
  fftw_destroy_plan(flt->fft_forward);

  fftw_free(flt->fft_vector_fft);
  fftw_free(flt->psf_fft);
}

/*
 * convolve()
 *
 * Convolve image with psf.
 *
 * flt - A pointer to a filter structure.
 * image - The image (must be the same size as the original psf image).
 * result - Pre-allocated storage for the result of the convolution.
 */
void convolve(filter *flt, double *image, double *result)
{
  int i;
  double c,r;

  /* Compute FFT of the image. */
  for(i=0;i<flt->image_size;i++){
    flt->fft_vector[i] = image[i];
  }
  fftw_execute(flt->fft_forward);

  /* Multiple by FFT of the PSF and compute inverse FFT. */
  for(i=0;i<flt->fft_size;i++){
    r = flt->fft_vector_fft[i][0] * flt->psf_fft[i][0] - flt->fft_vector_fft[i][1] * flt->psf_fft[i][1];
    c = flt->fft_vector_fft[i][0] * flt->psf_fft[i][1] + flt->fft_vector_fft[i][1] * flt->psf_fft[i][0];
    flt->fft_vector_fft[i][0] = r;
    flt->fft_vector_fft[i][1] = c;
  }
  fftw_execute(flt->fft_backward);

  /* Copy into result, */
  for(i=0;i<flt->image_size;i++){
    result[i] = flt->fft_vector[i] * flt->normalization;
  }  
}

/*
 * initialize()
 *
 * Set things up for FFT convolution.
 *
 * psf - the psf (x_size, y_size).
 * x_size - the size of the psf in x (slow dimension).
 * y_size - the size of the psf in y (fast dimension).
 * estimate - 0/1 to just use an estimated FFT plan. If you are only going to
 *            to use the FFT a few times this can be much faster.
 */
filter *initialize(double *psf, int x_size, int y_size, int estimate)
{
  int i;
  filter *flt;

  flt = (filter *)malloc(sizeof(filter));
  
  /* Initialize some variables. */
  flt->fft_size = x_size * (y_size/2 + 1);
  flt->image_size = x_size * y_size;
  
  flt->x_size = x_size;
  flt->y_size = y_size;
  flt->normalization = 1.0/((double)(x_size * y_size));

  /* Allocate storage. */
  flt->fft_vector = (double *)fftw_malloc(sizeof(double)*flt->image_size);
  flt->fft_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*flt->fft_size);
  flt->psf_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*flt->fft_size);

  /* Create FFT plans. */
  if (estimate){
    flt->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, flt->fft_vector, flt->fft_vector_fft, FFTW_ESTIMATE);
    flt->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, flt->fft_vector_fft, flt->fft_vector, FFTW_ESTIMATE);
  }
  else {
    flt->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, flt->fft_vector, flt->fft_vector_fft, FFTW_MEASURE);
    flt->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, flt->fft_vector_fft, flt->fft_vector, FFTW_MEASURE);    
  }

  /* Compute FFT of psf and save. */
  for(i=0;i<flt->image_size;i++){
    flt->fft_vector[i] = psf[i];
  }
  
  fftw_execute(flt->fft_forward);
  
  for(i=0;i<flt->fft_size;i++){
    flt->psf_fft[i][0] = flt->fft_vector_fft[i][0];
    flt->psf_fft[i][1] = flt->fft_vector_fft[i][1];
  }

  return flt;
}

/*
 * The MIT License
 *
 * Copyright (c) 2016 Zhuang Lab, Harvard University
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
