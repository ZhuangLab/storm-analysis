/*
 * C code for Poisson sparse deconvolution following 3DenseSTORM.
 *
 * As described in:
 * Ovesny et al., "High density 3D localization microscopy using
 * sparse support recovery", Optics Express, 2014.
 *
 * Hazen 11/19
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>

#include "../sa_library/ft_math.h"


/* 3denseSTORM Structure. */
typedef struct{
  int fft_size;
  int image_size;
  int number_psfs;
  int stale_Ax;

  double beta;
  double eta;
  double micro;
  double normalization;

  double *Ax;
  double *b_vec;
  double *d_vec;
  double *fft_vector;
  double *image;
  double *w_vec;
  double *y_vec;
  double *y_vec_h;

  double **e_vec;
  double **work1;  
  double **x_vec_h;
  double **x_vec_t;

  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *fft_vector_fft;
    
  fftw_complex **A;
  fftw_complex **G_inv;
  fftw_complex **work2;

} dsData;


/* Functions Declarations. */
void calculateAx(dsData *);
void cleanup(dsData *);
void getAx(dsData *, double *);
void getXVector(dsData *, double *, int);
dsData* initialize2D(double, double, double, int, int);
dsData* initialize3D(double, double, double, int, int, int);
void initializeA(dsData *, fftw_complex *, int);
void initializeGInv(dsData *, fftw_complex *, int);
void iterate(dsData *);
double l1Error(dsData *);
double l2Error(dsData *);
void newImage(dsData *, double *, double *);
void run(dsData *, int);


/* Functions */

/*
 * calculateAx()
 *
 * ds_data - A pointer to a dsData structure.
 */
void calculateAx(dsData *ds_data)
{
  int i;
  
  if (ds_data->stale_Ax){

    /* Compute Ax. */
    ftmComplexZero(ds_data->work2[0], ds_data->fft_size);
  
    for(i=0;i<ds_data->number_psfs;i++){

      // Compute FFT of x vector for each z plane.
      ftmDoubleCopy(ds_data->x_vec_h[i], ds_data->fft_vector, ds_data->image_size);
      fftw_execute(ds_data->fft_forward);

      // Multiply FFT of x vector by FFT of the PSF for this z plane.
      ftmComplexMultiplyAccum(ds_data->work2[0], ds_data->fft_vector_fft, ds_data->A[i], ds_data->fft_size, 0);
    }

    ftmComplexCopy(ds_data->work2[0], ds_data->fft_vector_fft, ds_data->fft_size);
    fftw_execute(ds_data->fft_backward);

    /* Store Ax. */
    ftmDoubleCopyNormalize(ds_data->fft_vector, ds_data->Ax, ds_data->normalization, ds_data->image_size);
  }
  
  ds_data->stale_Ax = 0;
}

/*
 * cleanup
 *
 * ds_data - A pointer to a dsData structure.
 */
void cleanup(dsData *ds_data)
{
  int i;
  
  free(ds_data->Ax);
  free(ds_data->b_vec);
  free(ds_data->d_vec);
  free(ds_data->image);
  free(ds_data->w_vec);
  free(ds_data->y_vec);
  free(ds_data->y_vec_h);

  for(i=0;i<ds_data->number_psfs;i++){
    free(ds_data->e_vec[i]);
    free(ds_data->work1[i]);
    free(ds_data->x_vec_h[i]);
    free(ds_data->x_vec_t[i]);
  }
  free(ds_data->e_vec);
  free(ds_data->work1);
  free(ds_data->x_vec_h);
  free(ds_data->x_vec_t);

  fftw_destroy_plan(ds_data->fft_backward);
  fftw_destroy_plan(ds_data->fft_forward);

  fftw_free(ds_data->fft_vector);  
  fftw_free(ds_data->fft_vector_fft);

  for(i=0;i<ds_data->number_psfs;i++){
    free(ds_data->A[i]);
    free(ds_data->work2[i]);
  }
  for(i=0;i<(ds_data->number_psfs*ds_data->number_psfs);i++){
    free(ds_data->G_inv[i]);
  }

  free(ds_data->A);
  free(ds_data->G_inv);
  free(ds_data->work2);

  free(ds_data);
}

/*
 * getAx()
 *
 * Copies the current Ax vector into user supplied storage.
 *
 * ds_data - A pointer to a dsData structure.
 * data - Storage for the Ax vector.
 */
void getAx(dsData *ds_data, double *data)
{
  calculateAx(ds_data);
  ftmDoubleCopy(ds_data->Ax, data, ds_data->image_size);
}

/*
 * getXVector()
 *
 * Copies the current x vector into user supplied storage. Also
 * converts back from [(x1,y1), (x2,y2), ..] to (x,y,z).
 *
 * ds_data - A pointer to a dsData structure.
 * data - Storage for the x vector.
 * compressed - Return x_vec_t (1) or x_vec_h (0).
 */
void getXVector(dsData *ds_data, double *data, int compressed)
{
  int i,j;
  double *t1;

  for(i=0;i<ds_data->number_psfs;i++){
    if (compressed){
      t1 = ds_data->x_vec_t[i];
    }
    else{
      t1 = ds_data->x_vec_h[i];
    }
    for(j=0;j<ds_data->image_size;j++){
      data[j*ds_data->number_psfs+i] = t1[j];
    }
  }
}

/*
 * initialize2D()
 *
 * Initialization for 2D analysis.
 *
 * beta - The beta parameter.
 * eta - The eta parameter.
 * micro - The micro parameter.
 * x_size - Size in x of data and psf.
 * y_size - Size in y of data and psf.
 */
dsData* initialize2D(double beta, double eta, double micro, int x_size, int y_size)
{
  return initialize3D(beta, eta, micro, x_size, y_size, 1);
}

/*
 * initialize3D()
 *
 * Initialization for 3D analysis. Note that setting the values 
 * in the A and G_inv arrays is done in a separate step.
 *
 * beta - The beta parameter.
 * eta - The eta parameter.
 * micro - The micro parameter.
 * x_size - Size in x of data and psf.
 * y_size - Size in y of data and psf.
 * z_size - The number of z planes.
 */
dsData* initialize3D(double beta, double eta, double micro, int x_size, int y_size, int z_size)
{
  int i;
  dsData *ds_data;

  ds_data = (dsData *)malloc(sizeof(dsData));

  ds_data->fft_size = x_size * (y_size/2+1);
  ds_data->image_size = x_size * y_size;
  ds_data->number_psfs = z_size;
  
  ds_data->normalization = 1.0/((double)(ds_data->image_size));
  ds_data->beta = beta;
  ds_data->eta = eta;
  ds_data->micro = micro;

  /* Allocate storage. */
  ds_data->Ax = (double *)malloc(sizeof(double) * ds_data->image_size);
  ds_data->b_vec = (double *)malloc(sizeof(double) * ds_data->image_size);
  ds_data->d_vec = (double *)malloc(sizeof(double) * ds_data->image_size);
  ds_data->image = (double *)malloc(sizeof(double) * ds_data->image_size);
  ds_data->w_vec = (double *)malloc(sizeof(double) * ds_data->image_size);
  ds_data->y_vec = (double *)malloc(sizeof(double) * ds_data->image_size);
  ds_data->y_vec_h = (double *)malloc(sizeof(double) * ds_data->image_size);

  ds_data->e_vec = (double **)malloc(sizeof(double *) * z_size);
  ds_data->work1 = (double **)malloc(sizeof(double *) * z_size);
  ds_data->x_vec_h = (double **)malloc(sizeof(double *) * z_size);
  ds_data->x_vec_t = (double **)malloc(sizeof(double *) * z_size);

  for(i=0;i<z_size;i++){
    ds_data->e_vec[i] = (double *)malloc(sizeof(double) * ds_data->image_size);
    ds_data->work1[i] = (double *)malloc(sizeof(double) * ds_data->image_size);
    ds_data->x_vec_h[i] = (double *)malloc(sizeof(double) * ds_data->image_size);
    ds_data->x_vec_t[i] = (double *)malloc(sizeof(double) * ds_data->image_size);
  }
  
  ds_data->fft_vector = (double *)fftw_malloc(sizeof(double) * ds_data->image_size);
  ds_data->fft_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ds_data->fft_size);

  ds_data->A = (fftw_complex **)malloc(sizeof(fftw_complex *) * z_size);
  ds_data->G_inv = (fftw_complex **)malloc(sizeof(fftw_complex *) * z_size * z_size);
  ds_data->work2 = (fftw_complex **)malloc(sizeof(fftw_complex *) * z_size);

  for(i=0;i<z_size;i++){
    ds_data->A[i] = (fftw_complex *)malloc(sizeof(fftw_complex) * ds_data->fft_size);
    ds_data->work2[i] = (fftw_complex *)malloc(sizeof(fftw_complex) * ds_data->fft_size);
  }

  for(i=0;i<(z_size*z_size);i++){
    ds_data->G_inv[i] = (fftw_complex *)malloc(sizeof(fftw_complex) * ds_data->fft_size);
  }
  
  /* Create FFT plans. */
  ds_data->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, ds_data->fft_vector, ds_data->fft_vector_fft, FFTW_MEASURE);
  ds_data->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, ds_data->fft_vector_fft, ds_data->fft_vector, FFTW_MEASURE);

  return ds_data;
}

/*
 * initializeA()
 *
 * Set initial values in the A array.
 *
 * ds_data - A pointer to a dsData structure.
 * data - Complex values to use for A.
 * index - Which A matrix to store data in.
 */
void initializeA(dsData *ds_data, fftw_complex *data, int index)
{
  ftmComplexCopy(data, ds_data->A[index], ds_data->fft_size);
}

/*
 * initializeGInv()
 *
 * Set initial values in the G_inv array.
 *
 * ds_data - A pointer to an admmData structure.
 * data - Complex values to use for G.
 * index - Which G matrix to store data in.
 */
void initializeGInv(dsData *ds_data, fftw_complex *data, int index)
{
  ftmComplexCopy(data, ds_data->G_inv[index], ds_data->fft_size);
}

/*
 * iterate()
 *
 * Performs one cycle of improvement.
 *
 * ds_data - A pointer to a dsData structure.
 */

void iterate(dsData *ds_data)
{
  int i,j,n_psfs;
  double t1, t2, t3;
  double *ev,*w1,*xvh,*xvt;

  n_psfs = ds_data->number_psfs;

  /*
   * Equation 6. 
   */
  /* Calculate At_ydb. */
  for(i=0;i<ds_data->image_size;i++){
    ds_data->fft_vector[i] = ds_data->y_vec_h[i] + ds_data->d_vec[i] - ds_data->b_vec[i];
  }
  fftw_execute(ds_data->fft_forward);
  ftmComplexCopy(ds_data->fft_vector_fft, ds_data->work2[0], ds_data->fft_size);
      
  for(i=0;i<n_psfs;i++){
    ftmComplexMultiply(ds_data->fft_vector_fft, ds_data->work2[0], ds_data->A[i], ds_data->fft_size, 1);
    fftw_execute(ds_data->fft_backward);
    ftmDoubleCopyNormalize(ds_data->fft_vector, ds_data->work1[i], ds_data->normalization, ds_data->image_size);
  }

  /* Calculate At_ydb + micro * (x-tilde + e). */
  for(i=0;i<n_psfs;i++){
    ev = ds_data->e_vec[i];
    w1 = ds_data->work1[i];
    xvt = ds_data->x_vec_t[i];
    for(j=0;j<ds_data->image_size;j++){
      ds_data->fft_vector[j] = w1[j] + ds_data->micro * (xvt[j] + ev[j]);
    }
    fftw_execute(ds_data->fft_forward);
    ftmComplexCopy(ds_data->fft_vector_fft, ds_data->work2[i], ds_data->fft_size);
  }

  /* Calculate updated x-hat. */
  for(i=0;i<n_psfs;i++){
    ftmComplexZero(ds_data->fft_vector_fft, ds_data->fft_size);
    for(j=0;j<n_psfs;j++){
      ftmComplexMultiplyAccum(ds_data->fft_vector_fft, ds_data->G_inv[i*n_psfs+j], ds_data->work2[j], ds_data->fft_size, 0);
    }
    fftw_execute(ds_data->fft_backward);
    ftmDoubleCopyNormalize(ds_data->fft_vector, ds_data->x_vec_h[i], ds_data->normalization, ds_data->image_size);
  }

  /*
   * Equation 7. 
   */
  ds_data->stale_Ax = 1;
  calculateAx(ds_data);
  t1 = 1.0/(2.0*ds_data->eta);
  for(i=0;i<ds_data->image_size;i++){
    t2 = 1.0 - ds_data->eta*(ds_data->Ax[i] + ds_data->b_vec[i] - ds_data->d_vec[i]);
    t3 = sqrt(t2*t2 + 4.0*ds_data->eta*ds_data->y_vec[i]);
    ds_data->y_vec_h[i] = t1*(-t2 + t3);
  }

  /*
   * Equation 8. 
   */
  for(i=0;i<n_psfs;i++){
    xvh = ds_data->x_vec_h[i];
    xvt = ds_data->x_vec_t[i];
    ev = ds_data->e_vec[i];
    for(j=0;j<ds_data->image_size;j++){
      t1 = xvh[j] - ev[j] - ds_data->w_vec[j];
      if (t1 <= 0.0){
	xvt[j] = 0.0;
      }
      else {
	xvt[j] = t1;
      }
    }
  }

  /*
   * Equation 9. 
   */
  for(i=0;i<ds_data->image_size;i++){
    ds_data->d_vec[i] = ds_data->d_vec[i] - (ds_data->Ax[i] + ds_data->b_vec[i] - ds_data->y_vec_h[i]);
  }

  /*
   * Equation 10. 
   */
  for(i=0;i<n_psfs;i++){
    xvh = ds_data->x_vec_h[i];
    xvt = ds_data->x_vec_t[i];
    ev = ds_data->e_vec[i];
    for(j=0;j<ds_data->image_size;j++){
      ev[j] = ev[j] - (xvh[j] - xvt[j]);
    }
  }

  /* Flag that we need to recaculate Ax (if getAx is called). */
  /* ds_data->stale_Ax = 1; */
}

/*
 * l1Error
 * 
 * Calculate l1 error term (sum of the x vector).
 *
 * ds_data - A pointer to a dsData structure.
 *
 * Return sum of the x vector.
 */
double l1Error(dsData *ds_data)
{
  int i,j;
  double sum;
  double *xv;
  
  sum = 0.0;
  for(i=0;i<ds_data->number_psfs;i++){
    xv = ds_data->x_vec_h[i];
    for(j=0;j<ds_data->image_size;j++){
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
 * ds_data - A pointer to a dsData structure.
 *
 * Return the error in the fit.
 */
double l2Error(dsData *ds_data)
{
  int i;
  double l2_error,t1;

  /* Compute Ax. */
  calculateAx(ds_data);

  /* Compute (Ax - b)^2. */
  l2_error = 0.0;
  for(i=0;i<ds_data->image_size;i++){
    t1 = ds_data->Ax[i] - ds_data->image[i];
    l2_error += t1*t1;
  }

  return sqrt(l2_error);
}

/*
 * newImage()
 *
 * Sets things up for the analysis of a new image.
 *
 * ds_data - A pointer to a dsData structure.
 * image - The image data.
 * background - Image background estimate.
 */
void newImage(dsData *ds_data, double *image, double *background)
{
  int i;

  /* Flag that we need to recaculate Ax (if getAx is called). */
  ds_data->stale_Ax = 1;
  
  /* Store image - background. */
  for(i=0;i<ds_data->image_size;i++){
    ds_data->image[i] = image[i] - background[i];
  }
  
  ftmDoubleCopy(image, ds_data->y_vec, ds_data->image_size);
  ftmDoubleCopy(background, ds_data->b_vec, ds_data->image_size);

  for(i=0;i<ds_data->image_size;i++){
    ds_data->w_vec[i] = ds_data->beta * sqrt(background[i]) / ds_data->micro;
  }
    
  /* d, y-hat vectors. */
  ftmDoubleCopy(image, ds_data->y_vec_h, ds_data->image_size);
  ftmDoubleZero(ds_data->d_vec, ds_data->image_size);

  /* zero e, x-hat, x-tidle vectors. */
  for(i=0;i<ds_data->number_psfs;i++){
    ftmDoubleZero(ds_data->e_vec[i], ds_data->image_size);
    ftmDoubleZero(ds_data->x_vec_h[i], ds_data->image_size);
    ftmDoubleZero(ds_data->x_vec_t[i], ds_data->image_size);
  }
}

/*
 * run()
 *
 * ds_data - A pointer to a dsData structure.
 * cycles - Number of iterations to perform.
 *
 * Performs a fixed number of cycles.
 */
void run(dsData *ds_data, int cycles)
{
  int i;

  for(i=0;i<cycles;i++){
    iterate(ds_data);
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
