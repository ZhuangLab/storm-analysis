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

  int fft_x_size;
  
  int x_size;
  int y_size;
  int z_size;
  
  double normalization;

  double *forward_vector;
  double *x_shift_c;
  double *x_shift_r;
  double *y_shift_c;
  double *y_shift_r;
  double *z_shift_c;
  double *z_shift_r;

  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *backward_vector;
  fftw_complex *psf_fft;
};
typedef struct psf_fft_struct psf_fft;

/* Function Declarations */
void calcShiftVector(double *, double *, double, int);
void cleanup(psf_fft *);
void getPSF(psf_fft *, double *, double, double, double);
psf_fft *initialize(double *, int, int, int);

/* Functions */

/*
 * calcShiftVector()
 *
 * Calculate the FFT shift vector.
 *
 * sr - The real part of the shift vector.
 * sc - The complex part of the shift vector.
 * dx - Shift delta (in pixels).
 * size - The size of the vector.
 */
void calcShiftVector(double *sr, double *sc, double dx, int size)
{
  int i, max_i;
  double t1, t2, tc, tr;

  sr[0] = 1.0;
  sc[0] = 0.0;
  t1 = M_PI * 2.0 * dx/((double)size);
  for(i=1;i<(size/2+1);i++){
    t2 = t1 * (double)i;
    tr = cos(t2);
    tc = sin(t2);
    sr[size-i] = tr;
    sr[i] = tr;
    sc[size-i] = -tc;
    sc[i] = tc;
  }
}

/*
 * cleanup()
 *
 * pfft - A pointer to a psf_fft structure.
 */
void cleanup(psf_fft *pfft)
{
  free(pfft->forward_vector);
  free(pfft->x_shift_c);
  free(pfft->x_shift_r);
  free(pfft->y_shift_c);
  free(pfft->y_shift_r);
  free(pfft->z_shift_c);
  free(pfft->z_shift_r);
 
  fftw_destroy_plan(pfft->fft_backward);
  fftw_destroy_plan(pfft->fft_forward);

  fftw_free(pfft->backward_vector);
  fftw_free(pfft->psf_fft);
}

/*
 * getPsf()
 *
 * Return the PSF translated by dx, dy, dz.
 *
 * pfft - A pointer to a psf_fft structure.
 * psf - Pre-allocated storage for the result.
 * dz - Displacement in z (in pixels).
 * dy - Displacement in y (in pixels).
 * dx - Displacement in x (in pixels).
 */
void getPSF(psf_fft *pfft, double *psf, double dz, double dy, double dx)
{
  int i, j, k, t1, t2, t3;
  int mid_z, size_xy;
  double c1, c2, r1, r2, *sxc, *sxr, *syc, *syr, *szc, *szr;

  /* Calculate FFT translation vectors. */
  sxc = pfft->x_shift_c;
  sxr = pfft->x_shift_r;
  calcShiftVector(sxr, sxc, dx, pfft->x_size);
  
  syc = pfft->y_shift_c;
  syr = pfft->y_shift_r;
  calcShiftVector(syr, syc, dy, pfft->y_size);
  
  szc = pfft->z_shift_c;
  szr = pfft->z_shift_r;
  calcShiftVector(szr, szc, dz, pfft->z_size);

  /* Translate */
  for(i=0;i<pfft->z_size;i++){
    t1 = i * (pfft->y_size * pfft->fft_x_size);
    for(j=0;j<pfft->y_size;j++){
      t2 = j * pfft->fft_x_size;
      for(k=0;k<pfft->fft_x_size;k++){
	t3 = t1 + t2 + k;
	
	r1 = pfft->psf_fft[t3][0];
	c1 = pfft->psf_fft[t3][1];

	/* z shift. */
	r2 = r1 * szr[i] - c1 * szc[i];
	c2 = r1 * szc[i] + c1 * szr[i];

	/* y shift. */
	r1 = r2 * syr[j] - c2 * syc[j];
	c1 = r2 * syc[j] + c2 * syr[j];

	/* x shift. */
	r2 = r1 * sxr[k] - c1 * sxc[k];
	c2 = r1 * sxc[k] + c1 * sxr[k];		

	pfft->backward_vector[t3][0] = r2;
	pfft->backward_vector[t3][1] = c2;
      }
    }
  }
  
  /* Do reverse transform. */
  fftw_execute(pfft->fft_backward);

  /* The 2D PSF is the middle plane of the 3D PSF. */
  size_xy = pfft->x_size*pfft->y_size;
  mid_z = size_xy * (pfft->z_size/2);
  
  for(i=0;i<size_xy;i++){
    psf[i] = pfft->forward_vector[mid_z+i] * pfft->normalization;
  }
  
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
 *
 * return - A psf_fft structure.
 */
psf_fft *initialize(double *psf, int z_size, int y_size, int x_size)
{
  int i;
  psf_fft *pfft;

  pfft = (psf_fft *)malloc(sizeof(psf_fft));
  
  /* Initialize some variables. */
  pfft->fft_size = z_size * y_size * (x_size/2 + 1);
  pfft->psf_size = z_size * y_size * x_size;

  pfft->fft_x_size = (x_size/2 + 1);
  pfft->x_size = x_size;
  pfft->y_size = y_size;
  pfft->z_size = z_size;
  
  pfft->normalization = 1.0/((double)(pfft->psf_size));

  /* Allocate storage. */
  pfft->x_shift_c = (double *)malloc(sizeof(double)*x_size);
  pfft->x_shift_r = (double *)malloc(sizeof(double)*x_size);
  pfft->y_shift_c = (double *)malloc(sizeof(double)*y_size);
  pfft->y_shift_r = (double *)malloc(sizeof(double)*y_size);
  pfft->z_shift_c = (double *)malloc(sizeof(double)*z_size);
  pfft->z_shift_r = (double *)malloc(sizeof(double)*z_size);

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
