/*
 * Use the ADMM approach to solve the lasso problem:
 *
 * minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1
 *
 * As described in:
 * Boyd et al., "Distributed Optimization and Statistical Learning 
 * via the Alternating Direction Method of Multipliers", Foundations
 * and Trends in Machine Learning, 2010.
 *
 * Notes:
 *  1. This has been specialized for the problem of image 
 *     deconvolution and makes heavy use of the FFT. This
 *     program cannot be used to solve the generic lasso problem.
 *
 * Hazen 11/19
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>

#include "../sa_library/ft_math.h"


/* ADMM Structure. */
typedef struct{
  int fft_size;
  int image_size;
  int number_psfs;
  int stale_Ax;
  
  double normalization;
  double rho;

  double *Ax;
  double *fft_vector;
  double *image;
  
  double **Atb;
  double **u_vector;
  double **x_vector;
  double **z_vector;

  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *fft_vector_fft;
    
  fftw_complex **A;
  fftw_complex **G_inv;
  fftw_complex **work1;

} admmData;


/* Functions Declarations. */
void calculateAx(admmData *);
void cleanup(admmData *);
void getAx(admmData *, double *);
void getXVector(admmData *, double *);
admmData* initialize2D(double, int, int);
admmData* initialize3D(double, int, int, int);
void initializeA(admmData *, fftw_complex *, int);
void initializeGInv(admmData *, fftw_complex *, int);
double l1Error(admmData *);
double l2Error(admmData *);
void iterate(admmData *, double, int);
void newImage(admmData *, double *);
void run(admmData *, double, int);


/* Functions */

/*
 * calculateAx()
 *
 * admm_data - A pointer to a admmData structure.
 */
void calculateAx(admmData *admm_data)
{
  int i;
  
  if (admm_data->stale_Ax){

    /* Compute Ax. */
    ftmComplexZero(admm_data->work1[0], admm_data->fft_size);
  
    for(i=0;i<admm_data->number_psfs;i++){

      // Compute FFT of x vector for each z plane.
      ftmDoubleCopy(admm_data->x_vector[i], admm_data->fft_vector, admm_data->image_size);
      fftw_execute(admm_data->fft_forward);

      // Multiply FFT of x vector by FFT of the PSF for this z plane.
      ftmComplexMultiplyAccum(admm_data->work1[0], admm_data->fft_vector_fft, admm_data->A[i], admm_data->fft_size, 0);
    }

    ftmComplexCopy(admm_data->work1[0], admm_data->fft_vector_fft, admm_data->fft_size);
    fftw_execute(admm_data->fft_backward);

    /* Store Ax. */
    ftmDoubleCopyNormalize(admm_data->fft_vector, admm_data->Ax, admm_data->normalization, admm_data->image_size);
  }
  
  admm_data->stale_Ax = 0;
}

/*
 * cleanup
 *
 * admm_data - A pointer to a admmData structure.
 */
void cleanup(admmData *admm_data)
{
  int i;

  free(admm_data->Ax);
  free(admm_data->image);
  
  for(i=0;i<admm_data->number_psfs;i++){
    free(admm_data->Atb[i]);
    free(admm_data->x_vector[i]);
    free(admm_data->z_vector[i]);
    free(admm_data->u_vector[i]);
  }
  free(admm_data->Atb);
  free(admm_data->x_vector);
  free(admm_data->z_vector);
  free(admm_data->u_vector);

  fftw_destroy_plan(admm_data->fft_backward);
  fftw_destroy_plan(admm_data->fft_forward);

  fftw_free(admm_data->fft_vector);
  fftw_free(admm_data->fft_vector_fft);
  
  for(i=0;i<admm_data->number_psfs;i++){
    free(admm_data->A[i]);
    free(admm_data->work1[i]);
  }
  for(i=0;i<(admm_data->number_psfs*admm_data->number_psfs);i++){
    free(admm_data->G_inv[i]);
  }

  free(admm_data->A);
  free(admm_data->G_inv);
  free(admm_data->work1);

  free(admm_data);
}

/*
 * getAx()
 *
 * Copies the current Ax vector into user supplied storage.
 *
 * admm_data - A pointer to a admmData structure.
 * data - Storage for the Ax vector.
 */
void getAx(admmData *admm_data, double *data)
{
  calculateAx(admm_data);
  ftmDoubleCopy(admm_data->Ax, data, admm_data->image_size);
}

/*
 * getXVector()
 *
 * Copies the current x vector into user supplied storage. Also
 * converts back from [(x1,y1), (x2,y2), ..] to (x,y,z).
 *
 * admm_data - A pointer to a admmData structure.
 * data - Storage for the x vector.
 */
void getXVector(admmData *admm_data, double *data)
{
  int i,j;
  double *t1;

  for(i=0;i<admm_data->number_psfs;i++){
    t1 = admm_data->x_vector[i];
    for(j=0;j<admm_data->image_size;j++){
      data[j*admm_data->number_psfs+i] = t1[j];
    }
  }
}

/*
 * initialize2D()
 *
 * Set things up for 2D analysis.
 *
 * rho - The rho parameter.
 * x_size - Size in x of data and psf.
 * y_size - Size in y of data and psf.
 */
admmData* initialize2D(double rho, int x_size, int y_size)
{
  return initialize3D(rho, x_size, y_size, 1);
}

/*
 * initialize3D()
 *
 * Set things up for 3D analysis. Note that setting the values in the
 * A and G_inv arrays is done in a separate step.
 *
 * rho - The rho parameter.
 * x_size - Size in x of data and psf.
 * y_size - Size in y of data and psf.
 * z_size - The number of z planes.
 */
admmData* initialize3D(double rho, int x_size, int y_size, int z_size)
{
  int i;
  admmData *admm_data;

  admm_data = (admmData *)malloc(sizeof(admmData));

  admm_data->fft_size = x_size * (y_size/2+1);
  admm_data->image_size = x_size * y_size;
  admm_data->number_psfs = z_size;
  
  admm_data->normalization = 1.0/((double)(admm_data->image_size));
  admm_data->rho = rho;

  /* Allocate storage. */
  admm_data->Ax = (double *)malloc(sizeof(double) * admm_data->image_size);
  admm_data->image = (double *)malloc(sizeof(double) * admm_data->image_size);

  admm_data->Atb = (double **)malloc(sizeof(double *) * z_size);
  admm_data->x_vector = (double **)malloc(sizeof(double *) * z_size);
  admm_data->z_vector = (double **)malloc(sizeof(double *) * z_size);
  admm_data->u_vector = (double **)malloc(sizeof(double *) * z_size);
  
  for(i=0;i<z_size;i++){
    admm_data->Atb[i] = (double *)malloc(sizeof(double) * admm_data->image_size);
    admm_data->x_vector[i] = (double *)malloc(sizeof(double) * admm_data->image_size);
    admm_data->z_vector[i] = (double *)malloc(sizeof(double) * admm_data->image_size);
    admm_data->u_vector[i] = (double *)malloc(sizeof(double) * admm_data->image_size);
  }

  admm_data->fft_vector = (double *)fftw_malloc(sizeof(double) * admm_data->image_size);
  admm_data->fft_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * admm_data->fft_size);

  admm_data->A = (fftw_complex **)malloc(sizeof(fftw_complex *) * z_size);
  admm_data->G_inv = (fftw_complex **)malloc(sizeof(fftw_complex *) * z_size * z_size);
  admm_data->work1 = (fftw_complex **)malloc(sizeof(fftw_complex *) * z_size);

  for(i=0;i<z_size;i++){
    admm_data->A[i] = (fftw_complex *)malloc(sizeof(fftw_complex) * admm_data->fft_size);
    admm_data->work1[i] = (fftw_complex *)malloc(sizeof(fftw_complex) * admm_data->fft_size);
  }

  for(i=0;i<(z_size*z_size);i++){
    admm_data->G_inv[i] = (fftw_complex *)malloc(sizeof(fftw_complex) * admm_data->fft_size);
  }
  
  /* Create FFT plans. */
  admm_data->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, admm_data->fft_vector, admm_data->fft_vector_fft, FFTW_MEASURE);
  admm_data->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, admm_data->fft_vector_fft, admm_data->fft_vector, FFTW_MEASURE);

  return admm_data;
}

/*
 * initializeA()
 *
 * Set initial values in the A array.
 *
 * admm_data - A pointer to an admmData structure.
 * data - Complex values to use for A.
 * index - Which A matrix to store data in.
 */
void initializeA(admmData *admm_data, fftw_complex *data, int index)
{
  ftmComplexCopy(data, admm_data->A[index], admm_data->fft_size);
}

/*
 * initializeGInv()
 *
 * Set initial values in the G_inv array.
 *
 * admm_data - A pointer to an admmData structure.
 * data - Complex values to use for G.
 * index - Which G matrix to store data in.
 */
void initializeGInv(admmData *admm_data, fftw_complex *data, int index)
{
  ftmComplexCopy(data, admm_data->G_inv[index], admm_data->fft_size);
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
  int i,j,n_psfs;
  double lr, t1, t2;
  double *Atb, *uv, *xv, *zv;

  /* Flag that we need to recaculate Ax (if getAx is called). */
  admm_data->stale_Ax = 1;

  n_psfs = admm_data->number_psfs;
  
  /* Calculate fft of (Atb + rho * (z - u)), store in work1. */
  for(i=0;i<n_psfs;i++){
    Atb = admm_data->Atb[i];
    uv = admm_data->u_vector[i];
    zv = admm_data->z_vector[i];
    for(j=0;j<admm_data->image_size;j++){
      admm_data->fft_vector[j] = Atb[j] + admm_data->rho * (zv[j] - uv[j]);
    }
    fftw_execute(admm_data->fft_forward);
    ftmComplexCopy(admm_data->fft_vector_fft, admm_data->work1[i], admm_data->fft_size);
  }

  /* Calculate x(k+1) by multiplying by G_inv. */
  for(i=0;i<n_psfs;i++){
    ftmComplexZero(admm_data->fft_vector_fft, admm_data->fft_size);
    for(j=0;j<n_psfs;j++){
      ftmComplexMultiplyAccum(admm_data->fft_vector_fft, admm_data->G_inv[i*n_psfs+j], admm_data->work1[j], admm_data->fft_size, 0);
    }
    fftw_execute(admm_data->fft_backward);
    ftmDoubleCopyNormalize(admm_data->fft_vector, admm_data->x_vector[i], admm_data->normalization, admm_data->image_size);
  }

  /* Calculate z(k+1). */
  lr = lambda/admm_data->rho;
  for(i=0;i<n_psfs;i++){
    uv = admm_data->u_vector[i];
    xv = admm_data->x_vector[i];
    zv = admm_data->z_vector[i];
    if (pos_only){
      for(j=0;j<admm_data->image_size;j++){
	t1 = xv[j] + uv[j] - lr;
	zv[j] = (t1 > 0.0) ? t1 : 0.0;
      }
    }
    else{
      for(j=0;j<admm_data->image_size;j++){
	t1 = xv[j] + uv[j];
	t2 = fabs(t1) - lr;
	if (t2 <= 0.0){
	  zv[j] = 0.0;
	}
	else {
	  zv[j] = (t1 > 0.0) ? t2 : -t2;
	}
      }
    }
  }
  
  /* Calculate u(k+1). */
  for(i=0;i<n_psfs;i++){
    uv = admm_data->u_vector[i];
    xv = admm_data->x_vector[i];
    zv = admm_data->z_vector[i];
    for(j=0;j<admm_data->image_size;j++){
      uv[j] += xv[j] - zv[j];
    }
  }
}

/*
 * l1Error
 * 
 * Calculate l1 error term (sum of the x vector).
 *
 * admm_data - A pointer to a fistaData structure.
 *
 * Return sum of the x vector.
 */
double l1Error(admmData *admm_data)
{
  int i,j;
  double sum;
  double *xv;
  
  sum = 0.0;
  for(i=0;i<admm_data->number_psfs;i++){
    xv = admm_data->x_vector[i];
    for(j=0;j<admm_data->image_size;j++){
      sum += fabs(xv[j]);
    }
  }

  return sum;
}

/*
 * l2Error
 *
 * Calculate l2 error (Ax - b).
 *
 * admm_data - A pointer to a admmData structure.
 *
 * Return the error in the fit.
 */
double l2Error(admmData *admm_data)
{
  int i;
  double l2_error,t1;

  /* Compute Ax. */
  calculateAx(admm_data);

  /* Compute (Ax - b)^2. */
  l2_error = 0.0;
  for(i=0;i<admm_data->image_size;i++){
    t1 = admm_data->Ax[i] - admm_data->image[i];
    l2_error += t1*t1;
  }

  return sqrt(l2_error);
}

/*
 * newImage()
 *
 * Sets things up for the analysis of a new image.
 *
 * image - The image data (with estimated background subtracted).
 */
void newImage(admmData *admm_data, double *image)
{
  int i;

  /* Flag that we need to recaculate Ax (if getAx is called). */
  admm_data->stale_Ax = 1;
  
  /* Save a copy of the image. */
  ftmDoubleCopy(image, admm_data->image, admm_data->image_size);
    
  /* Calculate Atb vector. */
  ftmDoubleCopy(image, admm_data->fft_vector, admm_data->image_size);
  fftw_execute(admm_data->fft_forward);
  ftmComplexCopy(admm_data->fft_vector_fft, admm_data->work1[0], admm_data->fft_size);
  
  for(i=0;i<admm_data->number_psfs;i++){
    ftmComplexMultiply(admm_data->fft_vector_fft, admm_data->work1[0], admm_data->A[i], admm_data->fft_size, 1);
    fftw_execute(admm_data->fft_backward);
    ftmDoubleCopyNormalize(admm_data->fft_vector, admm_data->Atb[i], admm_data->normalization, admm_data->image_size);
  }

  /* zero u, x, z vectors. */
  for(i=0;i<admm_data->number_psfs;i++){
    ftmDoubleZero(admm_data->u_vector[i], admm_data->image_size);
    ftmDoubleZero(admm_data->x_vector[i], admm_data->image_size);
    ftmDoubleZero(admm_data->z_vector[i], admm_data->image_size);
  }
}

/*
 * run()
 *
 * admm_data - A pointer to a admmData structure.
 * lambda - Lambda value.
 * cycles - Number of iterations to perform.
 *
 * Performs a fixed number of cycles at a fixed lambda.
 */
void run(admmData *admm_data, double lambda, int cycles)
{
  int i;

  for(i=0;i<cycles;i++){
    iterate(admm_data, lambda, 1);
  }
}


/*
 * The MIT License
 *
 * Copyright (c) 2019 Babcock Lab, Harvard University
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
