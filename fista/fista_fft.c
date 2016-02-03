/*
 * C library for FISTA (using FFT to compute Ax).
 *
 * Notes:
 *   1. Only works on square images.
 *   2. Image size should probably be a power of 2.
 *
 * Hazen 2/16
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall fista_fft.c
 *  gcc -shared -Wl,-soname,fista_fft.so.1 -o fista_fft.so.1.0.1 fista_fft.o -lc -lfftw3
 *  ln -s fista_fft.so.1.0.1 fista_fft.so
 *
 * Windows:
 *  gcc -c fista_fft.c
 *  gcc -shared -o fista_fft.dll -lfftw3
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>

/* Function Declarations */
void cleanup(void);
void getXVector(double *);
void initialize2D(double *, double, int);
void initialize3D(double *, double, int, int);
void iterate(double);
double l1Error(void);
double l2Error(void);
void newImage(double *);
void run(double, int);

/* Global Variables */
static int fft_size;
static int image_size;
static int number_psfs;
static double normalization;
static double time_step;
static double tk;
static double *fft_vector;
static double *image;
static double *x_vector;
static double *x_vector_old;
static double *y_vector;

static fftw_plan fft_forward;
static fftw_plan fft_backward;

static fftw_complex *Ax_fft;
static fftw_complex *fft_vector_fft;
static fftw_complex *image_fft;
static fftw_complex *psf_fft;

/* Functions */

/*
 * cleanup()
 */
void cleanup(void)
{
  free(fft_vector);
  free(image);
  free(x_vector);
  free(x_vector_old);
  free(y_vector);

  fftw_destroy_plan(fft_forward);
  fftw_destroy_plan(fft_backward);

  fftw_free(Ax_fft);
  fftw_free(fft_vector_fft);
  fftw_free(image_fft);
  fftw_free(psf_fft);
}

/*
 * getXVector()
 *
 * Copies the current x vector into user supplied storage. Also
 * converts back from (z,x,y) to (x,y,z).
 *
 * data - storage for the x vector.
 */
void getXVector(double *data)
{
  int i,j,n;

  for(i=0;i<number_psfs;i++){
    n = i*image_size;
    for(j=0;j<image_size;j++){
      data[j*number_psfs+i] = x_vector[n+j];
    }
  }
}

/*
 * initialize2D()
 *
 * Set things up for 2D analysis.
 *
 * psf - the psf (also size x size)
 * t_step - the time step to use.
 * xy_size - the x and y dimensions of data and psf
 */
void initialize2D(double *psf, double t_step, int xy_size)
{
  initialize3D(psf, t_step, xy_size, 1);
}

/*
 * initialize3D()
 *
 * Set things up for 3D analysis.
 *
 * psf - the psf (xy_size, xy_size, z_size)
 * t_step - the time step to use.
 * xy_size - the x and y dimensions of data and psf
 * z_size - the z dimension of the psf.
 */
void initialize3D(double *psf, double t_step, int xy_size, int z_size)
{
  int i,j,n;
  double temp;
  
  /* Initialize some variables. */
  fft_size = xy_size * (xy_size/2 + 1);
  image_size = xy_size * xy_size;
  normalization = 1.0/((double)(image_size));
  number_psfs = z_size;
  
  /* Allocate storage. */
  image = (double *)malloc(sizeof(double)*image_size);
  x_vector = (double *)malloc(sizeof(double)*z_size*image_size);
  x_vector_old = (double *)malloc(sizeof(double)*z_size*image_size);
  y_vector = (double *)malloc(sizeof(double)*z_size*image_size);

  fft_vector = (double *)fftw_malloc(sizeof(double)*image_size);
  Ax_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*fft_size);
  fft_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*fft_size);
  image_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*fft_size);
  psf_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*z_size*fft_size);

  /* Create FFT plans. */
  fft_forward = fftw_plan_dft_r2c_2d(xy_size, xy_size, fft_vector, fft_vector_fft, FFTW_MEASURE);
  fft_backward = fftw_plan_dft_c2r_2d(xy_size, xy_size, fft_vector_fft, fft_vector, FFTW_MEASURE);

  /* 
     Compute FFTs of the psfs and save in psf_fft. 
     Note: The input psfs are indexed (x,y,z) but internally we
           use (z,x,y) as this is more convenient.
  */
  for(i=0;i<z_size;i++){
    n = i*fft_size;
    temp = 0.0;
    for(j=0;j<image_size;j++){
      fft_vector[j] = psf[j*z_size+i];
      temp += psf[j*z_size+i];
    }
    fftw_execute(fft_forward);
    for(j=0;j<fft_size;j++){
      psf_fft[n+j][0] = fft_vector_fft[j][0];
      psf_fft[n+j][1] = fft_vector_fft[j][1];
    }
  }

  /* Copy psf into x_vector for debugging. */
  if (0){
    printf("\n");
    for(i=0;i<z_size;i++){
      n = i*image_size;
      for(j=0;j<image_size;j++){
	x_vector[n+j] = psf[j*z_size+i];
      }
    }
  }

  /* Calculate optimal time step. */
  time_step = 0.0;
  for(i=0;i<(z_size*fft_size);i++){
    temp = psf_fft[i][0]*psf_fft[i][0] + psf_fft[i][1]*psf_fft[i][1];
    if(temp>time_step){
      time_step = temp;
    }
  }
  time_step = 1.0/(2.0*time_step);
  printf("Optimal time step is: %f\n", time_step);
  
  time_step = t_step;
}

/*
 * iterate()
 *
 * Performs a single cycle of optimization.
 *
 * lambda - l1 error term weigth.
 */
void iterate(double lambda)
{
  int i,j,n,o;
  double lt,new_tk,t1,t2;
  
  /* Copy current x vector into old x vector. */
  for(i=0;i<(number_psfs*image_size);i++){
    x_vector_old[i] = x_vector[i];
  }

  /*
   * This is the plk(y) step in the FISTA algorithm.
   */
  
  /* Compute Ax_fft (n,n). x is generic here, and does not
     implicitly refer to x_vector. */
  for(i=0;i<fft_size;i++){
    Ax_fft[i][0] = 0.0;
    Ax_fft[i][1] = 0.0;
  }

  for(i=0;i<number_psfs;i++){

    // Compute FFT of y vector for each z plane.
    n = i*image_size;
    for(j=0;j<image_size;j++){
      fft_vector[j] = y_vector[n+j];
    }
    fftw_execute(fft_forward);

    // Multiply FFT of y vector by FFT of the PSF for this z plane.
    o = i*fft_size;
    for(j=0;j<fft_size;j++){
      Ax_fft[j][0] += fft_vector_fft[j][0]*psf_fft[o+j][0] - fft_vector_fft[j][1]*psf_fft[o+j][1];
      Ax_fft[j][1] += fft_vector_fft[j][0]*psf_fft[o+j][1] + fft_vector_fft[j][1]*psf_fft[o+j][0];
    }
  }
  
  /* Compute Ax_fft - b_fft (image_fft) (n,n). */
  for(i=0;i<fft_size;i++){
    Ax_fft[i][0] -= image_fft[i][0];
    Ax_fft[i][1] -= image_fft[i][1];
  }

  /* Compute x = y - At(Ax-b) (z,n,n). */
  for(i=0;i<number_psfs;i++){

    // Compute inverse FFT of At(Ax-b) for each image plane.
    o = i*fft_size;    
    for(j=0;j<fft_size;j++){
      fft_vector_fft[j][0] = psf_fft[o+j][0]*Ax_fft[j][0] + psf_fft[o+j][1]*Ax_fft[j][1];
      fft_vector_fft[j][1] = psf_fft[o+j][0]*Ax_fft[j][1] - psf_fft[o+j][1]*Ax_fft[j][0];
    }
    fftw_execute(fft_backward);

    // Update x vector.
    n = i*image_size;    
    t1 = 2.0*time_step*normalization;
    for(j=0;j<image_size;j++){
      x_vector[n+j] = y_vector[n+j] - t1*fft_vector[j];
    }
  }
  
  /* Shrink x vector. */
  lt = time_step*lambda;
  for(i=0;i<(number_psfs*image_size);i++){
    t1 = x_vector[i];
    t2 = fabs(x_vector[i]) - lt;
    if (t2 <= 0.0){
      x_vector[i] = 0.0;
    }
    else{
      x_vector[i] = (t1 > 0.0) ? t2 : -t2;
    }
  }

  /* 
   * Update tk term step.
   */
  new_tk = 0.5*(1.0 + sqrt(1.0 + 4.0*tk*tk));

  /* 
   * Compute new y vector step.
   */
  t1 = (tk - 1.0)/new_tk;
  for(i=0;i<(number_psfs*image_size);i++){
    y_vector[i] = x_vector[i] + t1*(x_vector[i] - x_vector_old[i]);
  }

  /* 
   * Update tk term step.
   */
  tk = new_tk;
}

/*
 * l1Error
 *
 * Return sum of the x vector.
 */
double l1Error(void)
{
  int i;
  double sum;

  sum = 0.0;
  for(i=0;i<(number_psfs*image_size);i++){
    sum += fabs(x_vector[i]);
  }

  return sum;
}

/*
 * l2Error
 *
 * Return the error in the fit.
 */
double l2Error(void)
{
  int i,j,n,o;
  double l2_error,t1;
    
  /* Compute Ax_fft (n,n). */
  for(i=0;i<fft_size;i++){
    Ax_fft[i][0] = 0.0;
    Ax_fft[i][1] = 0.0;
  }
  
  for(i=0;i<number_psfs;i++){

    // Compute FFT of x vector for each z plane.
    n = i*image_size;
    for(j=0;j<image_size;j++){
      fft_vector[j] = x_vector[n+j];
    }
    fftw_execute(fft_forward);

    // Multiply FFT of x vector by FFT of the PSF for this z plane.
    o = i*fft_size;
    for(j=0;j<fft_size;j++){
      Ax_fft[j][0] += fft_vector_fft[j][0]*psf_fft[o+j][0] - fft_vector_fft[j][1]*psf_fft[o+j][1];
      Ax_fft[j][1] += fft_vector_fft[j][0]*psf_fft[o+j][1] + fft_vector_fft[j][1]*psf_fft[o+j][0];
    }
  }

  /* Compute Ax. */
  for(i=0;i<fft_size;i++){
    fft_vector_fft[i][0] = Ax_fft[i][0];
    fft_vector_fft[i][1] = Ax_fft[i][1];
  }
  fftw_execute(fft_backward);

  /* Compute (Ax - b)^2. */
  l2_error = 0.0;
  for(i=0;i<image_size;i++){
    t1 = fft_vector[i]*normalization - image[i];
    l2_error += t1*t1;
  }

  return sqrt(l2_error);
}

/*
 * newImage()
 *
 * Initialize with a new image.
 *
 * data - The image data.
 */
void newImage(double *data)
{
  int i;
  
  tk = 1.0;
  
  /* Save a copy of the image. */
  for(i=0;i<image_size;i++){
    image[i] = data[i];
  }
  
  /* Compute FFT of image. */
  for(i=0;i<image_size;i++){
    fft_vector[i] = data[i];
  }
  fftw_execute(fft_forward);
  for(i=0;i<fft_size;i++){
    image_fft[i][0] = fft_vector_fft[i][0];
    image_fft[i][1] = fft_vector_fft[i][1];
  }

  /* Initialize x_vector, y_vectors */
  for(i=0;i<(number_psfs*image_size);i++){
    x_vector[i] = 0.0;
    x_vector_old[i] = 0.0;
    y_vector[i] = 0.0;
  }
}

/*
 * run()
 *
 * Performs a fixed number of cycles at a fixed lambda.
 */
void run(double lambda, int cycles)
{
  int i;

  for(i=0;i<cycles;i++){
    iterate(lambda);
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
