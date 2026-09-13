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

#include "ft_math.h"


/* Structures & Types */
struct filter_struct {
  int fft_size;
  int image_size;
  int x_size;
  int y_size;

  double max_diff;

  double *fft_vector;
  double *old_image;
  
  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *fft_vector_fft;
  fftw_complex *psf_fft;
};
typedef struct filter_struct filter;

/*
 * A bank of filters that all get applied to the same image.
 *
 * The point of this is the forward FFT. Applying N filters with N separate
 * 'filter' structures transforms the image N times, and every one of those
 * transforms is identical. Here the image is transformed once and the result
 * is reused, so the cost goes from 2N transforms to N+1.
 *
 * image_fft has to be kept separate from work_fft because convolve() does its
 * complex multiply in place, which would destroy the transform we are trying
 * to reuse.
 */
struct filter_bank_struct {
  int fft_size;
  int image_size;
  int n_filters;
  int x_size;
  int y_size;

  double *fft_vector;

  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *image_fft;
  fftw_complex *work_fft;
  fftw_complex **psf_ffts;
};
typedef struct filter_bank_struct filterBank;

/* Function Declarations */
void cleanup(filter *);
void cleanupBank(filterBank *);
void convolve(filter *, double *, double *);
void convolveBank(filterBank *, double *, double *);
void convolveMemo(filter *, double *, double *);
filter *initialize(double *, double, int, int, int);
filterBank *initializeBank(double *, int, int, int, int);

/* Functions */


/*
 * cleanup()
 *
 * flt - A pointer to a filter structure.
 */
void cleanup(filter *flt)
{
  if (flt->old_image != NULL){
    free(flt->old_image);
  }
  
  fftw_free(flt->fft_vector);
  
  fftw_destroy_plan(flt->fft_backward);
  fftw_destroy_plan(flt->fft_forward);

  fftw_free(flt->fft_vector_fft);
  fftw_free(flt->psf_fft);

  free(flt);
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
  /* Compute FFT of the image. */
  ftmDoubleCopy(image, flt->fft_vector, flt->image_size);
  fftw_execute(flt->fft_forward);

  /* Multiple by FFT of the PSF and compute inverse FFT. */
  ftmComplexMultiply(flt->fft_vector_fft, flt->fft_vector_fft, flt->psf_fft, flt->fft_size, 0);
  fftw_execute(flt->fft_backward);

  /* Copy into result, */
  ftmDoubleCopy(flt->fft_vector, result, flt->image_size);
}


/*
 * convolveMemo()
 *
 * Convolve image with psf, but only if image is different enough from
 * the previous image, otherwise just return the previous result.
 *
 * flt - A pointer to a filter structure.
 * image - The image (must be the same size as the original psf image).
 * result - Pre-allocated storage for the result of the convolution.
 */
void convolveMemo(filter *flt, double *image, double *result)
{
  int i,different;

  different = 0;
  
  /* Just check, don't copy so that we don't slowly drift.. */
  for(i=0;i<flt->image_size;i++){
    if(fabs(image[i] - flt->old_image[i]) > flt->max_diff){
      different = 1;
      break;
    }
  }

  if (different){

    /* Copy into old image. */
    for(i=0;i<flt->image_size;i++){
      flt->old_image[i] = image[i];
    }
    //ftmDoubleCopy(image, flt->old_image, flt->image_size);
    
    /* Do the convolution. */
    convolve(flt, image, result);
  }
  else {
    
    /* Otherwise just return the previous result. */
    ftmDoubleCopy(flt->fft_vector, result, flt->image_size);
  }
}


/*
 * initialize()
 *
 * Set things up for FFT convolution.
 *
 * psf - the psf (x_size, y_size).
 * max_diff - if this it not zero configure for memoization of results.
 * x_size - the size of the psf in x (slow dimension).
 * y_size - the size of the psf in y (fast dimension).
 * estimate - 0/1 to just use an estimated FFT plan. If you are only going to
 *            to use the FFT a few times this can be much faster.
 */
filter *initialize(double *psf, double max_diff, int x_size, int y_size, int estimate)
{
  int i;
  double normalization;
  filter *flt;

  flt = (filter *)malloc(sizeof(filter));
  
  /* Initialize some variables. */
  flt->fft_size = x_size * (y_size/2 + 1);
  flt->image_size = x_size * y_size;
  
  flt->x_size = x_size;
  flt->y_size = y_size;
  normalization = 1.0/((double)(x_size * y_size));

  /* Check whether we are memoizing. */
  if(max_diff > 0.0){
    flt->max_diff = max_diff;
    flt->old_image = (double *)malloc(sizeof(double)*flt->image_size);

    for(i=0;i<flt->image_size;i++){
      flt->old_image[i] = 0;
    }
  }
  else{
    flt->max_diff = 0.0;
    flt->old_image = NULL;
  }

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

  ftmDoubleCopy(psf, flt->fft_vector, flt->image_size);
  fftw_execute(flt->fft_forward);
  ftmComplexCopyNormalize(flt->fft_vector_fft, flt->psf_fft, normalization, flt->fft_size);

  return flt;
}


/*
 * cleanupBank()
 *
 * bank - A pointer to a filterBank structure.
 */
void cleanupBank(filterBank *bank)
{
  int i;

  for(i=0;i<bank->n_filters;i++){
    fftw_free(bank->psf_ffts[i]);
  }
  free(bank->psf_ffts);

  fftw_free(bank->fft_vector);
  fftw_free(bank->image_fft);
  fftw_free(bank->work_fft);

  fftw_destroy_plan(bank->fft_backward);
  fftw_destroy_plan(bank->fft_forward);

  free(bank);
}


/*
 * convolveBank()
 *
 * Convolve image with every psf in the bank. This is the same arithmetic
 * that convolve() does, in the same order, just with the forward transform
 * hoisted out of the loop.
 *
 * bank - A pointer to a filterBank structure.
 * image - The image (must be the same size as the psfs the bank was made with).
 * results - Pre-allocated storage for n_filters images, in psf order.
 */
void convolveBank(filterBank *bank, double *image, double *results)
{
  int i;

  /* Compute FFT of the image, once. */
  ftmDoubleCopy(image, bank->fft_vector, bank->image_size);
  fftw_execute(bank->fft_forward);

  for(i=0;i<bank->n_filters;i++){

    /* Multiply by FFT of this psf and compute inverse FFT. */
    ftmComplexMultiply(bank->work_fft, bank->image_fft, bank->psf_ffts[i], bank->fft_size, 0);
    fftw_execute(bank->fft_backward);

    /* Copy into the matching slice of result. */
    ftmDoubleCopy(bank->fft_vector, results + i*bank->image_size, bank->image_size);
  }
}


/*
 * initializeBank()
 *
 * Set things up for FFT convolution with several psfs at once.
 *
 * psfs - n_filters psfs, contiguous, each (x_size, y_size).
 * n_filters - The number of psfs.
 * x_size - the size of the psfs in x (slow dimension).
 * y_size - the size of the psfs in y (fast dimension).
 * estimate - 0/1 to just use an estimated FFT plan.
 */
filterBank *initializeBank(double *psfs, int n_filters, int x_size, int y_size, int estimate)
{
  int i;
  double normalization;
  filterBank *bank;

  bank = (filterBank *)malloc(sizeof(filterBank));

  bank->fft_size = x_size * (y_size/2 + 1);
  bank->image_size = x_size * y_size;
  bank->n_filters = n_filters;
  bank->x_size = x_size;
  bank->y_size = y_size;

  normalization = 1.0/((double)(x_size * y_size));

  /* Allocate storage. */
  bank->fft_vector = (double *)fftw_malloc(sizeof(double)*bank->image_size);
  bank->image_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*bank->fft_size);
  bank->work_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*bank->fft_size);

  bank->psf_ffts = (fftw_complex **)malloc(sizeof(fftw_complex *)*n_filters);
  for(i=0;i<n_filters;i++){
    bank->psf_ffts[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*bank->fft_size);
  }

  /*
   * Create FFT plans. This has to happen before anything is written into
   * fft_vector, as FFTW_MEASURE overwrites its arrays while planning.
   *
   * One pair of plans serves the whole bank, rather than a pair per psf.
   */
  if (estimate){
    bank->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, bank->fft_vector, bank->image_fft, FFTW_ESTIMATE);
    bank->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, bank->work_fft, bank->fft_vector, FFTW_ESTIMATE);
  }
  else {
    bank->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, bank->fft_vector, bank->image_fft, FFTW_MEASURE);
    bank->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, bank->work_fft, bank->fft_vector, FFTW_MEASURE);
  }

  /* Compute FFT of each psf and save. */
  for(i=0;i<n_filters;i++){
    ftmDoubleCopy(psfs + i*bank->image_size, bank->fft_vector, bank->image_size);
    fftw_execute(bank->fft_forward);
    ftmComplexCopyNormalize(bank->image_fft, bank->psf_ffts[i], normalization, bank->fft_size);
  }

  return bank;
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
