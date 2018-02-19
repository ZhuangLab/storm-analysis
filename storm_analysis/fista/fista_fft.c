/*
 * C library for FISTA (using FFT to compute Ax).
 *
 * Notes:
 *   1. Image size should probably be a power of 2.
 *
 * Hazen 2/16
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>


/* FISTA Structure */
typedef struct{
  int fft_size;
  int image_size;
  int number_psfs;
  double normalization;
  double time_step;
  double tk;
  double *fft_vector;
  double *image;
  double *x_vector;
  double *x_vector_old;
  double *y_vector;

  fftw_plan fft_forward;
  fftw_plan fft_backward;

  fftw_complex *Ax_fft;
  fftw_complex *fft_vector_fft;
  fftw_complex *image_fft;
  fftw_complex *psf_fft;
} fistaData;


/* Function Declarations */
void cleanup(fistaData *);
void getXVector(fistaData *, double *);
fistaData* initialize2D(double *, double, int, int);
fistaData* initialize3D(double *, double, int, int, int);
void iterate(fistaData *, double);
double l1Error(fistaData *);
double l2Error(fistaData *);
void newImage(fistaData *, double *);
void run(fistaData *, double, int);


/* Functions */

/*
 * cleanup()
 *
 * fista_data - A pointer to a fistaData structure.
 */
void cleanup(fistaData *fista_data)
{
  free(fista_data->image);
  free(fista_data->x_vector);
  free(fista_data->x_vector_old);
  free(fista_data->y_vector);

  fftw_destroy_plan(fista_data->fft_forward);
  fftw_destroy_plan(fista_data->fft_backward);

  fftw_free(fista_data->Ax_fft);
  fftw_free(fista_data->fft_vector);
  fftw_free(fista_data->fft_vector_fft);
  fftw_free(fista_data->image_fft);
  fftw_free(fista_data->psf_fft);

  free(fista_data);
}

/*
 * getXVector()
 *
 * Copies the current x vector into user supplied storage. Also
 * converts back from (z,x,y) to (x,y,z).
 *
 * fista_data - A pointer to a fistaData structure.
 * data - Storage for the x vector.
 */
void getXVector(fistaData *fista_data, double *data)
{
  int i,j,n;

  for(i=0;i<fista_data->number_psfs;i++){
    n = i*fista_data->image_size;
    for(j=0;j<fista_data->image_size;j++){
      data[j*fista_data->number_psfs+i] = fista_data->x_vector[n+j];
    }
  }
}

/*
 * initialize2D()
 *
 * Set things up for 2D analysis.
 *
 * psf - The psf (also x_size x y_size)
 * t_step - The time step to use.
 * x_size - Size in x of data and psf.
 * y_size - Size in y of data and psf.
 */
fistaData* initialize2D(double *psf, double t_step, int x_size, int y_size)
{
  return initialize3D(psf, t_step, x_size, y_size, 1);
}

/*
 * initialize3D()
 *
 * Set things up for 3D analysis.
 *
 * psf - The psf (x_size, y_size, z_size)
 * t_step - the time step to use.
 * x_size - Size in x of data and psf.
 * y_size - Size in y of data and psf.
 * z_size - The z dimension of the psf.
 */
fistaData* initialize3D(double *psf, double t_step, int x_size, int y_size, int z_size)
{
  int i,j,n;
  double temp;
  fistaData *fista_data;

  fista_data = (fistaData *)malloc(sizeof(fistaData));

  /* Initialize some variables. */
  fista_data->fft_size = x_size * (y_size/2 + 1);
  fista_data->image_size = x_size * y_size;
  fista_data->normalization = 1.0/((double)(fista_data->image_size));
  fista_data->number_psfs = z_size;
  
  /* Allocate storage. */
  fista_data->image = (double *)malloc(sizeof(double) * fista_data->image_size);
  fista_data->x_vector = (double *)malloc(sizeof(double) * z_size * fista_data->image_size);
  fista_data->x_vector_old = (double *)malloc(sizeof(double) * z_size * fista_data->image_size);
  fista_data->y_vector = (double *)malloc(sizeof(double) * z_size * fista_data->image_size);

  fista_data->fft_vector = (double *)fftw_malloc(sizeof(double) * fista_data->image_size);
  fista_data->Ax_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fista_data->fft_size);
  fista_data->fft_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fista_data->fft_size);
  fista_data->image_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fista_data->fft_size);
  fista_data->psf_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * z_size * fista_data->fft_size);

  /* Create FFT plans. */
  fista_data->fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, fista_data->fft_vector, fista_data->fft_vector_fft, FFTW_MEASURE);
  fista_data->fft_backward = fftw_plan_dft_c2r_2d(x_size, y_size, fista_data->fft_vector_fft, fista_data->fft_vector, FFTW_MEASURE);

  /* 
     Compute FFTs of the psfs and save in psf_fft. 
     Note: The input psfs are indexed (x,y,z) but internally we
           use (z,x,y) as this is more convenient.
  */
  for(i=0;i<z_size;i++){
    n = i*fista_data->fft_size;
    temp = 0.0;
    for(j=0;j<fista_data->image_size;j++){
      fista_data->fft_vector[j] = psf[j*z_size+i];
      temp += psf[j*z_size+i];
    }
    fftw_execute(fista_data->fft_forward);
    for(j=0;j<fista_data->fft_size;j++){
      fista_data->psf_fft[n+j][0] = fista_data->fft_vector_fft[j][0];
      fista_data->psf_fft[n+j][1] = fista_data->fft_vector_fft[j][1];
    }
  }
  
  /* Copy psf into x_vector for debugging. */
  if (0){
    printf("\n");
    for(i=0;i<z_size;i++){
      n = i*fista_data->image_size;
      for(j=0;j<fista_data->image_size;j++){
	fista_data->x_vector[n+j] = psf[j*z_size+i];
      }
    }
  }

  /* Calculate optimal time step. */
  fista_data->time_step = 0.0;
  for(i=0;i<(z_size*fista_data->fft_size);i++){
    temp = fista_data->psf_fft[i][0] * fista_data->psf_fft[i][0] + fista_data->psf_fft[i][1] * fista_data->psf_fft[i][1];
    if(temp>fista_data->time_step){
      fista_data->time_step = temp;
    }
  }
  
  fista_data->time_step = 1.0/(2.0*fista_data->time_step);
  printf("Optimal time step is: %f\n", fista_data->time_step);
  
  fista_data->time_step = t_step;

  return fista_data;
}

/*
 * iterate()
 *
 * Performs a single cycle of optimization.
 *
 * fista_data - A pointer to a fistaData structure.
 * lambda - l1 error term weigth.
 */
void iterate(fistaData *fista_data, double lambda)
{
  int i,j,n,o;
  double lt,new_tk,t1,t2;
  
  /* Copy current x vector into old x vector. */
  for(i=0;i<(fista_data->number_psfs * fista_data->image_size);i++){
    fista_data->x_vector_old[i] = fista_data->x_vector[i];
  }

  /*
   * This is the plk(y) step in the FISTA algorithm.
   */
  
  /* Compute Ax_fft (n,n). x is generic here, and does not
     implicitly refer to x_vector. */
  for(i=0;i<fista_data->fft_size;i++){
    fista_data->Ax_fft[i][0] = 0.0;
    fista_data->Ax_fft[i][1] = 0.0;
  }

  for(i=0;i<fista_data->number_psfs;i++){

    // Compute FFT of y vector for each z plane.
    n = i*fista_data->image_size;
    for(j=0;j<fista_data->image_size;j++){
      fista_data->fft_vector[j] = fista_data->y_vector[n+j];
    }
    fftw_execute(fista_data->fft_forward);

    // Multiply FFT of y vector by FFT of the PSF for this z plane.
    o = i*fista_data->fft_size;
    for(j=0;j<fista_data->fft_size;j++){
      fista_data->Ax_fft[j][0] += fista_data->fft_vector_fft[j][0] * fista_data->psf_fft[o+j][0] - fista_data->fft_vector_fft[j][1] * fista_data->psf_fft[o+j][1];
      fista_data->Ax_fft[j][1] += fista_data->fft_vector_fft[j][0] * fista_data->psf_fft[o+j][1] + fista_data->fft_vector_fft[j][1] * fista_data->psf_fft[o+j][0];
    }
  }
  
  /* Compute Ax_fft - b_fft (image_fft) (n,n). */
  for(i=0;i<fista_data->fft_size;i++){
    fista_data->Ax_fft[i][0] -= fista_data->image_fft[i][0];
    fista_data->Ax_fft[i][1] -= fista_data->image_fft[i][1];
  }

  /* Compute x = y - At(Ax-b) (z,n,n). */
  for(i=0;i<fista_data->number_psfs;i++){

    // Compute inverse FFT of At(Ax-b) for each image plane.
    o = i*fista_data->fft_size;    
    for(j=0;j<fista_data->fft_size;j++){
      fista_data->fft_vector_fft[j][0] = fista_data->psf_fft[o+j][0] * fista_data->Ax_fft[j][0] + fista_data->psf_fft[o+j][1] * fista_data->Ax_fft[j][1];
      fista_data->fft_vector_fft[j][1] = fista_data->psf_fft[o+j][0] * fista_data->Ax_fft[j][1] - fista_data->psf_fft[o+j][1] * fista_data->Ax_fft[j][0];
    }
    fftw_execute(fista_data->fft_backward);

    // Update x vector.
    n = i*fista_data->image_size;    
    t1 = 2.0*fista_data->time_step*fista_data->normalization;
    for(j=0;j<fista_data->image_size;j++){
      fista_data->x_vector[n+j] = fista_data->y_vector[n+j] - t1*fista_data->fft_vector[j];
    }
  }
  
  /* Shrink x vector. */
  lt = fista_data->time_step*lambda;
  for(i=0;i<(fista_data->number_psfs * fista_data->image_size);i++){
    t1 = fista_data->x_vector[i];
    t2 = fabs(fista_data->x_vector[i]) - lt;
    if (t2 <= 0.0){
      fista_data->x_vector[i] = 0.0;
    }
    else{
      fista_data->x_vector[i] = (t1 > 0.0) ? t2 : -t2;
    }
  }

  /* 
   * Update tk term step.
   */
  new_tk = 0.5*(1.0 + sqrt(1.0 + 4.0 * fista_data->tk * fista_data->tk));

  /* 
   * Compute new y vector step.
   */
  t1 = (fista_data->tk - 1.0)/new_tk;
  for(i=0;i<(fista_data->number_psfs * fista_data->image_size);i++){
    fista_data->y_vector[i] = fista_data->x_vector[i] + t1*(fista_data->x_vector[i] - fista_data->x_vector_old[i]);
  }

  /* 
   * Update tk term step.
   */
  fista_data->tk = new_tk;
}

/*
 * l1Error
 *
 * fista_data - A pointer to a fistaData structure.
 * Return sum of the x vector.
 */
double l1Error(fistaData *fista_data)
{
  int i;
  double sum;

  sum = 0.0;
  for(i=0;i<(fista_data->number_psfs * fista_data->image_size);i++){
    sum += fabs(fista_data->x_vector[i]);
  }

  return sum;
}

/*
 * l2Error
 *
 * fista_data - A pointer to a fistaData structure.
 * Return the error in the fit.
 */
double l2Error(fistaData *fista_data)
{
  int i,j,n,o;
  double l2_error,t1;
    
  /* Compute Ax_fft (n,n). */
  for(i=0;i<fista_data->fft_size;i++){
    fista_data->Ax_fft[i][0] = 0.0;
    fista_data->Ax_fft[i][1] = 0.0;
  }
  
  for(i=0;i<fista_data->number_psfs;i++){

    // Compute FFT of x vector for each z plane.
    n = i*fista_data->image_size;
    for(j=0;j<fista_data->image_size;j++){
      fista_data->fft_vector[j] = fista_data->x_vector[n+j];
    }
    fftw_execute(fista_data->fft_forward);

    // Multiply FFT of x vector by FFT of the PSF for this z plane.
    o = i*fista_data->fft_size;
    for(j=0;j<fista_data->fft_size;j++){
      fista_data->Ax_fft[j][0] += fista_data->fft_vector_fft[j][0] * fista_data->psf_fft[o+j][0] - fista_data->fft_vector_fft[j][1] * fista_data->psf_fft[o+j][1];
      fista_data->Ax_fft[j][1] += fista_data->fft_vector_fft[j][0] * fista_data->psf_fft[o+j][1] + fista_data->fft_vector_fft[j][1] * fista_data->psf_fft[o+j][0];
    }
  }

  /* Compute Ax. */
  for(i=0;i<fista_data->fft_size;i++){
    fista_data->fft_vector_fft[i][0] = fista_data->Ax_fft[i][0];
    fista_data->fft_vector_fft[i][1] = fista_data->Ax_fft[i][1];
  }
  fftw_execute(fista_data->fft_backward);

  /* Compute (Ax - b)^2. */
  l2_error = 0.0;
  for(i=0;i<fista_data->image_size;i++){
    t1 = fista_data->fft_vector[i] * fista_data->normalization - fista_data->image[i];
    l2_error += t1*t1;
  }

  return sqrt(l2_error);
}

/*
 * newImage()
 *
 * Initialize with a new image.
 *
 * fista_data - A pointer to a fistaData structure.
 * data - The image data.
 */
void newImage(fistaData *fista_data, double *data)
{
  int i;
  
  fista_data->tk = 1.0;
  
  /* Save a copy of the image. */
  for(i=0;i<fista_data->image_size;i++){
    fista_data->image[i] = data[i];
  }
  
  /* Compute FFT of image. */
  for(i=0;i<fista_data->image_size;i++){
    fista_data->fft_vector[i] = data[i];
  }
  fftw_execute(fista_data->fft_forward);
  for(i=0;i<fista_data->fft_size;i++){
    fista_data->image_fft[i][0] = fista_data->fft_vector_fft[i][0];
    fista_data->image_fft[i][1] = fista_data->fft_vector_fft[i][1];
  }

  /* Initialize x_vector, y_vectors */
  for(i=0;i<(fista_data->number_psfs * fista_data->image_size);i++){
    fista_data->x_vector[i] = 0.0;
    fista_data->x_vector_old[i] = 0.0;
    fista_data->y_vector[i] = 0.0;
  }
}

/*
 * run()
 *
 * fista_data - A pointer to a fistaData structure.
 * Performs a fixed number of cycles at a fixed lambda.
 */
void run(fistaData *fista_data, double lambda, int cycles)
{
  int i;

  for(i=0;i<cycles;i++){
    iterate(fista_data, lambda);
  }
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
