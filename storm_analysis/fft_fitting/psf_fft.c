/*
 * C library for using a FFT to shift a PSF in x,y and z
 * as well as taking it's derivatives along these axises.
 *
 * Hazen 2/16
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall psf_fft.c
 *  gcc -shared -Wl,-soname,psf_fft.so.1 -o psf_fft.so.1.0.1 psf_fft.o -lc -lfftw3
 *  ln -s psf_fft.so.1.0.1 psf_fft.so
 *
 * Windows:
 *  gcc -c -O3 psf_fft.c
 *  gcc -shared -o psf_fft.dll psf_fft.o -lfftw3-3 c:\path\to\libfftw3-3.dll
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>

/* Structures & Types */
struct psf_fft_struct {
  int fft_size;
  int psf_size;
  
  int x_size;
  int y_size;
  int z_size;
  
  double normalization;

  double *forward_vector;
  
  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *backward_vector;
  fftw_complex *psf_fft;
};
typedef struct psf_fft_struct psf_fft;

/* Function Declarations */
void cleanup(psf_fft *);
psf_fft *initialize(double *, int, int, int);

/* Functions */

/*
 * cleanup()
 *
 * filter - A pointer to a filter structure.
 */
void cleanup(psf_fft *pfft)
{
  free(pfft->forward_vector);

  fftw_destroy_plan(pfft->fft_backward);
  fftw_destroy_plan(pfft->fft_forward);

  fftw_free(pfft->backward_vector);
  fftw_free(pfft->psf_fft);
}

/*
 * initialize()
 *
 * Set things up for FFT based PSF calculations.
 *
 * Note: In the psf array, the z axis is the slowest followed
 *       by the y axis and then the x axis.
 *
 * psf - The psf (z_size, y_size, x_size).
 * z_size - The size of the psf in z (slowest dimension).
 * y_size - The size of the psf in y.
 * x_size - The size of the psf in x (slow dimension).
 */
psf_fft *initialize(double *psf, int z_size, int y_size, int x_size)
{
  int i;
  psf_fft *pfft;

  pfft = (psf_fft *)malloc(sizeof(psf_fft));
  
  /* Initialize some variables. */
  pfft->fft_size = z_size * y_size * (x_size/2 + 1);
  pfft->psf_size = z_size * y_size * x_size;

  pfft->x_size = x_size;
  pfft->y_size = y_size;
  pfft->z_size = z_size;
  
  pfft->normalization = 1.0/((double)(pfft->psf_size));

  /* Allocate storage. */
  pfft->forward_vector = (double *)fftw_malloc(sizeof(double)*pfft->psf_size);
  pfft->backward_vector = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pfft->fft_size);
  pfft->psf_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pfft->fft_size);

  /* Create FFT plans. */
  pfft->fft_forward = fftw_plan_dft_r2c_3d(z_size, y_size, x_size, pfft->forward_vector, pfft->backward_vector, FFTW_MEASURE);
  pfft->fft_backward = fftw_plan_dft_c2r_3d(z_size, y_size, x_size, pfft->backward_vector, pfft->forward_vector, FFTW_MEASURE);

  /* Compute FFT of psf and save. */
  for(i=0;i<pfft->psf_size;i++){
    pfft->forward_vector[i] = psf[i];
  }
  
  fftw_execute(pfft->fft_forward);
  
  for(i=0;i<pfft->fft_size;i++){
    pfft->psf_fft[i][0] = pfft->backward_vector[i][0];
    pfft->psf_fft[i][1] = pfft->backward_vector[i][1];
  }

  return pfft;
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
