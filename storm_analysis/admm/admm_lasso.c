/*
 * Use the ADMM approach to solve the lasso problem:
 *
 * minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1
 *
 * As explained described in this paper:
 *  http://www.stanford.edu/~boyd/papers/admm_distr_stats.html
 *
 * And this Matlab program:
 *  http://www.stanford.edu/~boyd/papers/admm/lasso/lasso.html
 *
 * Note that this has been specialized for the problem of
 * image deconvolution and makes heavy use of the FFT. This
 * program cannot be used to solve the generic lasso problem.
 *
 * The FFT package that is used is available here:
 *  http://www.fftw.org/
 *
 * The debian package is: fftw3-dev
 *
 * Hazen 04/14
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>

/* Functions. */
void cleanup(void);
void getXVector(double *);
void initialize(double *, double, int);
void iterate(double, int);
void newImage(double *);

/* Global variables. */
static int image_size;
static double normalization;
static double rho;

static int *mask;
static double *Atb;
static double *x_vector;
static double *z_vector;
static double *u_vector;

static fftw_plan fft_backward;
static fftw_plan fft_forward;

static double *a_vector;
static fftw_complex *a_vector_fft;
static fftw_complex *m1_fft;
static fftw_complex *psf_fft;


/*
 * cleanup
 */
void cleanup(void)
{
  free(mask);
  free(Atb);
  free(a_vector);
  free(x_vector);
  free(z_vector);
  free(u_vector);

  fftw_destroy_plan(fft_backward);
  fftw_destroy_plan(fft_forward);

  fftw_free(a_vector_fft);
  fftw_free(m1_fft);
  fftw_free(psf_fft);
}

/*
 * getXVector()
 *
 * Returns the current x vector in user supplied storage.
 *
 * xv_copy - Pre-allocated storage to copy the x vector into.
 */
void getXVector(double *xv_copy)
{
  int i;

  for(i=0;i<image_size;i++){
    xv_copy[i] = x_vector[i];
  }
}

/*
 * initialize()
 *
 * Set things up for analysis.
 *
 * psf - The psf, a 1D array with a total size of size*size.
 * rho_param - The rho parameter.
 * size - The X and Y dimensions of psf.
 */
void initialize(double *psf, double rho_param, int size)
{
  int i;

  image_size = size*size;
  normalization = 1.0/((double)image_size);
  rho = rho_param;

  /* Allocate storage. */
  mask = (int *)malloc(sizeof(int)*image_size);
  Atb = (double *)malloc(sizeof(double)*image_size);
  x_vector = (double *)malloc(sizeof(double)*image_size);
  z_vector = (double *)malloc(sizeof(double)*image_size);
  u_vector = (double *)malloc(sizeof(double)*image_size);

  a_vector = (double *)fftw_malloc(sizeof(double)*image_size);
  a_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*image_size);
  m1_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*image_size);
  psf_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*image_size);

  /* Create FFT plans. */
  fft_forward = fftw_plan_dft_r2c_2d(size, size, a_vector, a_vector_fft, FFTW_MEASURE);
  fft_backward = fftw_plan_dft_c2r_2d(size, size, a_vector_fft, a_vector, FFTW_MEASURE);  

  /* 
   * Compute FFT of the PSF & save in psf_fft. 
   */
  for(i=0;i<image_size;i++){
    a_vector[i] = psf[i];
  }
  fftw_execute(fft_forward);
  for(i=0;i<image_size;i++){
    psf_fft[i][0] = a_vector_fft[i][0];
    psf_fft[i][1] = a_vector_fft[i][1];
  }

  /* 
   * Compute m1_fft vector & valid mask. This relies on the psf computation have
   * been performed above so that a_vector_fft contains the fft of the psf.
   */
  for(i=0;i<image_size;i++){
    a_vector_fft[i][0] = a_vector_fft[i][0]*a_vector_fft[i][0] + a_vector_fft[i][1]*a_vector_fft[i][1];
    a_vector_fft[i][1] = 0.0;
  }
  fftw_execute(fft_backward);
  for(i=0;i<image_size;i++){
    a_vector[i] = a_vector[i]*normalization;
  }
  a_vector[0] += rho;
  fftw_execute(fft_forward);
  for(i=0;i<image_size;i++){
    m1_fft[i][0] = a_vector_fft[i][0];
    m1_fft[i][1] = a_vector_fft[i][1];
    if ((fabs(m1_fft[i][0]) > 0.0)||(fabs(m1_fft[i][1]) > 0.0)){
      mask[i] = 1;
    }
    else{
      mask[i] = 0;
    }
  }
}

/*
 * iterate()
 *
 * Performs one cycle of improvement.
 *
 * lambda - The lambda value to use.
 * pos_only - If 1, all values of z < lambda/rho are set to 0.0.
 */
void iterate(double lambda, int pos_only)
{
  int i;
  double imag, lr, mag, real, t1, t2;

  /* Calculate fft of (Atb + rho * (z - u)). */
  for(i=0;i<image_size;i++){
    a_vector[i] = Atb[i] + rho * (z_vector[i] - u_vector[i]);
  }
  fftw_execute(fft_forward);
  
  /* Divide by m1_fft (where valid). */
  for(i=0;i<image_size;i++){
    if(mask[i]){
      mag = 1.0/(m1_fft[i][0]*m1_fft[i][0] + m1_fft[i][1]*m1_fft[i][1]);
      real = a_vector_fft[i][0]*m1_fft[i][0] + a_vector_fft[i][1]*m1_fft[i][1];
      imag = a_vector_fft[i][1]*m1_fft[i][0] - a_vector_fft[i][0]*m1_fft[i][1];
      a_vector_fft[i][0] = mag * real;
      a_vector_fft[i][1] = mag * imag;
    }
  }

  /* Calculate x. */
  fftw_execute(fft_backward);
  for(i=0;i<image_size;i++){
    x_vector[i] = a_vector[i]*normalization;
  }

  /* Calculate z. */
  lr = lambda/rho;
  if (pos_only){
    for(i=0;i<image_size;i++){
      t1 = x_vector[i] + u_vector[i] - lr;
      z_vector[i] = (t1 > 0.0) ? t1 : 0.0;
    }
  }
  else{
    for(i=0;i<image_size;i++){
      t1 = x_vector[i] + u_vector[i];
      t2 = fabs(t1) - lr;
      if (t2 <= 0.0){
	z_vector[i] = 0.0;
      }
      else {
	z_vector[i] = (t1 > 0.0) ? t2 : -t2;
      }
    }
  }

  /* Update u. */
  for(i=0;i<image_size;i++){
    u_vector[i] += x_vector[i] - z_vector[i];
  }
}

/*
 * newImage()
 *
 * Sets things up for the analysis of a new image.
 *
 * data - The image data.
 */
void newImage(double *data)
{
  int i;
  double imag, real;

  /* Calculate Atb vector. */
  for(i=0;i<image_size;i++){
    a_vector[i] = data[i];
  }
  fftw_execute(fft_forward);
  for(i=0;i<image_size;i++){
    real = psf_fft[i][0]*a_vector_fft[i][0] + psf_fft[i][1]*a_vector_fft[i][1];
    imag = psf_fft[i][0]*a_vector_fft[i][1] - psf_fft[i][1]*a_vector_fft[i][0];
    a_vector_fft[i][0] = real;
    a_vector_fft[i][1] = imag;
  }
  fftw_execute(fft_backward);
  for(i=0;i<image_size;i++){
    Atb[i] = a_vector[i]*normalization;
  }

  /* zero x, z, u vectors. */
  for(i=0;i<image_size;i++){
    x_vector[i] = 0.0;
    z_vector[i] = 0.0;
    u_vector[i] = 0.0;
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
