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
 * Notes:
 *  1. This has been specialized for the problem of image 
 *     deconvolution and makes heavy use of the FFT. This
 *     program cannot be used to solve the generic lasso problem.
 *
 *  2. Currently limited to 2D. The plan is to add 3D support at
 *     some point in the future.
 *
 * Hazen 02/18
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>


/* ADMM Structure. */
typedef struct{
  int image_size;
  double normalization;
  double rho;

  int *mask;
  double *Atb;
  double *a_vector;
  double *u_vector;
  double *x_vector;
  double *z_vector;

  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *a_vector_fft;
  fftw_complex *m1_fft;
  fftw_complex *psf_fft;
} admmData;


/* Functions. */
void cleanup(admmData *);
void getXVector(admmData *, double *);
admmData* initialize(double *, double, int, int);
void iterate(admmData *, double, int);
void newImage(admmData *, double *);

/*
 * cleanup
 */
void cleanup(admmData *admm_data)
{
  free(admm_data->mask);
  free(admm_data->Atb);
  free(admm_data->a_vector);
  free(admm_data->x_vector);
  free(admm_data->z_vector);
  free(admm_data->u_vector);

  fftw_destroy_plan(admm_data->fft_backward);
  fftw_destroy_plan(admm_data->fft_forward);

  fftw_free(admm_data->a_vector_fft);
  fftw_free(admm_data->m1_fft);
  fftw_free(admm_data->psf_fft);

  free(admm_data);
}

/*
 * getXVector()
 *
 * Returns the current x vector in user supplied storage.
 *
 * xv_copy - Pre-allocated storage to copy the x vector into.
 */
void getXVector(admmData *admm_data, double *xv_copy)
{
  int i;

  for(i=0;i<admm_data->image_size;i++){
    xv_copy[i] = admm_data->x_vector[i];
  }
}

/*
 * initialize()
 *
 * Set things up for analysis.
 *
 * psf - The psf (x_size by y_size)
 * rho_param - The rho parameter.
 * x_size - Size in x of data and psf.
 * y_size - Size in y of data and psf.
 */
admmData* initialize(double *psf, double rho_param, int x_size, int y_size)
{
  int i;
  admmData *admm_data;

  admm_data = (admmData *)malloc(sizeof(admmData));

  admm_data->image_size = x_size * y_size;
  admm_data->normalization = 1.0/((double)(admm_data->image_size));
  admm_data->rho = rho_param;

  /* Allocate storage. */
  admm_data->mask = (int *)malloc(sizeof(int) * admm_data->image_size);
  admm_data->Atb = (double *)malloc(sizeof(double) * admm_data->image_size);
  admm_data->x_vector = (double *)malloc(sizeof(double) * admm_data->image_size);
  admm_data->z_vector = (double *)malloc(sizeof(double) * admm_data->image_size);
  admm_data->u_vector = (double *)malloc(sizeof(double) * admm_data->image_size);

  admm_data->a_vector = (double *)fftw_malloc(sizeof(double) * admm_data->image_size);
  admm_data->a_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * admm_data->image_size);
  admm_data->m1_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * admm_data->image_size);
  admm_data->psf_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * admm_data->image_size);

  /* Create FFT plans. */
  admm_data->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, admm_data->a_vector, admm_data->a_vector_fft, FFTW_MEASURE);
  admm_data->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, admm_data->a_vector_fft, admm_data->a_vector, FFTW_MEASURE);

  /* 
   * Compute FFT of the PSF & save in psf_fft. 
   */
  for(i=0;i<admm_data->image_size;i++){
    admm_data->a_vector[i] = psf[i];
  }
  fftw_execute(admm_data->fft_forward);
  for(i=0;i<admm_data->image_size;i++){
    admm_data->psf_fft[i][0] = admm_data->a_vector_fft[i][0];
    admm_data->psf_fft[i][1] = admm_data->a_vector_fft[i][1];
  }

  /* 
   * Compute m1_fft vector & valid mask. This relies on the psf computation have
   * been performed above so that a_vector_fft contains the fft of the psf.
   */
  for(i=0;i<admm_data->image_size;i++){
    admm_data->a_vector_fft[i][0] = admm_data->a_vector_fft[i][0] * admm_data->a_vector_fft[i][0] + admm_data->a_vector_fft[i][1] * admm_data->a_vector_fft[i][1];
    admm_data->a_vector_fft[i][1] = 0.0;
  }
  fftw_execute(admm_data->fft_backward);
  for(i=0;i<admm_data->image_size;i++){
    admm_data->a_vector[i] = admm_data->a_vector[i] * admm_data->normalization;
  }
  admm_data->a_vector[0] += admm_data->rho;
  fftw_execute(admm_data->fft_forward);
  for(i=0;i<admm_data->image_size;i++){
    admm_data->m1_fft[i][0] = admm_data->a_vector_fft[i][0];
    admm_data->m1_fft[i][1] = admm_data->a_vector_fft[i][1];
    if ((fabs(admm_data->m1_fft[i][0]) > 0.0)||(fabs(admm_data->m1_fft[i][1]) > 0.0)){
      admm_data->mask[i] = 1;
    }
    else{
      admm_data->mask[i] = 0;
    }
  }

  return admm_data;
}

/*
 * iterate()
 *
 * Performs one cycle of improvement.
 *
 * lambda - The lambda value to use.
 * pos_only - If 1, all values of z < lambda/rho are set to 0.0.
 */
void iterate(admmData *admm_data, double lambda, int pos_only)
{
  int i;
  double imag, lr, mag, real, t1, t2;

  /* Calculate fft of (Atb + rho * (z - u)). */
  for(i=0;i<admm_data->image_size;i++){
    admm_data->a_vector[i] = admm_data->Atb[i] + admm_data->rho * (admm_data->z_vector[i] - admm_data->u_vector[i]);
  }
  fftw_execute(admm_data->fft_forward);
  
  /* Divide by m1_fft (where valid). */
  for(i=0;i<admm_data->image_size;i++){
    if(admm_data->mask[i]){
      mag = 1.0/(admm_data->m1_fft[i][0] * admm_data->m1_fft[i][0] + admm_data->m1_fft[i][1] * admm_data->m1_fft[i][1]);
      real = admm_data->a_vector_fft[i][0] * admm_data->m1_fft[i][0] + admm_data->a_vector_fft[i][1] * admm_data->m1_fft[i][1];
      imag = admm_data->a_vector_fft[i][1] * admm_data->m1_fft[i][0] - admm_data->a_vector_fft[i][0] * admm_data->m1_fft[i][1];
      admm_data->a_vector_fft[i][0] = mag * real;
      admm_data->a_vector_fft[i][1] = mag * imag;
    }
  }

  /* Calculate x. */
  fftw_execute(admm_data->fft_backward);
  for(i=0;i<admm_data->image_size;i++){
    admm_data->x_vector[i] = admm_data->a_vector[i] * admm_data->normalization;
  }

  /* Calculate z. */
  lr = lambda/admm_data->rho;
  if (pos_only){
    for(i=0;i<admm_data->image_size;i++){
      t1 = admm_data->x_vector[i] + admm_data->u_vector[i] - lr;
      admm_data->z_vector[i] = (t1 > 0.0) ? t1 : 0.0;
    }
  }
  else{
    for(i=0;i<admm_data->image_size;i++){
      t1 = admm_data->x_vector[i] + admm_data->u_vector[i];
      t2 = fabs(t1) - lr;
      if (t2 <= 0.0){
	admm_data->z_vector[i] = 0.0;
      }
      else {
	admm_data->z_vector[i] = (t1 > 0.0) ? t2 : -t2;
      }
    }
  }

  /* Update u. */
  for(i=0;i<admm_data->image_size;i++){
    admm_data->u_vector[i] += admm_data->x_vector[i] - admm_data->z_vector[i];
  }
}

/*
 * newImage()
 *
 * Sets things up for the analysis of a new image.
 *
 * data - The image data.
 */
void newImage(admmData *admm_data, double *data)
{
  int i;
  double imag, real;

  /* Calculate Atb vector. */
  for(i=0;i<admm_data->image_size;i++){
    admm_data->a_vector[i] = data[i];
  }
  fftw_execute(admm_data->fft_forward);
  for(i=0;i<admm_data->image_size;i++){
    real = admm_data->psf_fft[i][0] * admm_data->a_vector_fft[i][0] + admm_data->psf_fft[i][1] * admm_data->a_vector_fft[i][1];
    imag = admm_data->psf_fft[i][0] * admm_data->a_vector_fft[i][1] - admm_data->psf_fft[i][1] * admm_data->a_vector_fft[i][0];
    admm_data->a_vector_fft[i][0] = real;
    admm_data->a_vector_fft[i][1] = imag;
  }
  fftw_execute(admm_data->fft_backward);
  for(i=0;i<admm_data->image_size;i++){
    admm_data->Atb[i] = admm_data->a_vector[i] * admm_data->normalization;
  }

  /* zero x, z, u vectors. */
  for(i=0;i<admm_data->image_size;i++){
    admm_data->x_vector[i] = 0.0;
    admm_data->z_vector[i] = 0.0;
    admm_data->u_vector[i] = 0.0;
  }
}


/*
 * The MIT License
 *
 * Copyright (c) 2018 Babcock Lab, Harvard University
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
