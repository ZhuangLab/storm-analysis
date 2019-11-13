/*
 * Fourier transform and related math.
 *
 * Hazen 11/19
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>

#include "ft_math.h"


/*
 * ftmComplexCopy()
 *
 * Copy fftw_complex vectors (if they are not already the same).
 *
 * s1 - Pointer to source fftw_complex vector.
 * d1 - Pointer to destination fftw_complex vector.
 * size - Size of d1 and s1.
 */
void ftmComplexCopy(fftw_complex *s1, fftw_complex *d1, int size)
{
  memcpy(d1, s1, (sizeof(fftw_complex) * size));
}


/*
 * ftmComplexCopyNormalize()
 *
 * Copy fftw_complex vectors (if they are not already the same).
 *
 * s1 - Pointer to source fftw_complex vector.
 * d1 - Pointer to destination fftw_complex vector.
 * norm - Normalization constant (s1 will be multiplied by this).
 * size - Size of d1 and s1.
 */
void ftmComplexCopyNormalize(fftw_complex *s1, fftw_complex *d1, double norm, int size)
{
  int i;

  for(i=0;i<size;i++){
    d1[i][0] = s1[i][0]*norm;
    d1[i][1] = s1[i][1]*norm;
  }
}


/*
 * ftmComplexMultiply()
 *
 * Multiply two fftw_complex vectors together.
 *
 * d1 - Pointer to fftw_complex vector to save results in.
 * s1 - Pointer to first fftw_complex vector.
 * s2 - Pointer to second fftw_complex vector.
 * size - Size of d1, s1, s2.
 * conj - Multiply using complex conjugate of s2;
 */
void ftmComplexMultiply(fftw_complex *d1, fftw_complex *s1, fftw_complex *s2, int size, int conj)
{
  int i;
  double c,r;

  /* Use c,r intermediate variables so that this will also work inplace. */
  if (conj){
    for(i=0;i<size;i++){
      r = s1[i][0]*s2[i][0] + s1[i][1]*s2[i][1];
      c = s1[i][1]*s2[i][0] - s1[i][0]*s2[i][1];
      d1[i][0] = r;
      d1[i][1] = c;
    }
  }
  else{
    for(i=0;i<size;i++){
      r = s1[i][0]*s2[i][0] - s1[i][1]*s2[i][1];
      c = s1[i][1]*s2[i][0] + s1[i][0]*s2[i][1];
      d1[i][0] = r;
      d1[i][1] = c;
    }
  }
}


/*
 * ftmComplexMultiplyAccum()
 *
 * Multiply two fftw_complex vectors together and add to d1.
 *
 * d1 - Pointer to fftw_complex vector to save results in.
 * s1 - Pointer to first fftw_complex vector.
 * s2 - Pointer to second fftw_complex vector.
 * size - Size of d1, s1, s2.
 * conj - Multiply using complex conjugate of s2;
 */
void ftmComplexMultiplyAccum(fftw_complex *d1, fftw_complex *s1, fftw_complex *s2, int size, int conj)
{
  int i;

  if (conj){
    for(i=0;i<size;i++){
      d1[i][0] += s1[i][0]*s2[i][0] + s1[i][1]*s2[i][1];
      d1[i][1] += s1[i][1]*s2[i][0] - s1[i][0]*s2[i][1];
    }
  }
  else{
    for(i=0;i<size;i++){
      d1[i][0] += s1[i][0]*s2[i][0] - s1[i][1]*s2[i][1];
      d1[i][1] += s1[i][1]*s2[i][0] + s1[i][0]*s2[i][1];
    }
  }
}


/*
 * ftmComplexZero()
 *
 * Set all the values of a complex vector to zero.
 *
 * v1 - Pointer to fftw_complex vector.
 * size - Size of v1.
 */
void ftmComplexZero(fftw_complex *v1, int size)
{
  int i;

  for(i=0;i<size;i++){
    v1[i][0] = 0.0;
    v1[i][1] = 0.0;
  }
}


/*
 * ftmDoubleCopy()
 *
 * Copy double vectors (if they are not already the same).
 *
 * s1 - Pointer to source double vector.
 * d1 - Pointer to destination double vector.
 * size - Size of d1 and s1.
 */
void ftmDoubleCopy(double *s1, double *d1, int size)
{
  memcpy(d1, s1, (sizeof(double) * size));
}


/*
 * ftmDoubleCopyNormalize()
 *
 * Copy and normalize a double vector.
 *
 * s1 - Pointer to source double vector.
 * d1 - Pointer to destination double vector.
 * norm - Normalization constant (s1 will be multiplied by this).
 * size - Size of d1 and s1.
 */
void ftmDoubleCopyNormalize(double *s1, double *d1, double norm, int size)
{
  int i;

  for(i=0;i<size;i++){
    d1[i] = s1[i]*norm;
  }
}


/*
 * ftmDoubleNormalize()
 *
 * Normalize a double vector.
 *
 * v1 - Pointer to double vector to normalize.
 * norm - Normalization constant (v1 will be multiplied by this).
 * size - Size of v1.
 */
/*
void ftmDoubleNormalize(double *v1, double norm, int size)
{
  int i;

  for(i=0;i<size;i++){
    v1[i] = v1[i]*norm;
  }
}
*/


/*
 * ftmDoubleZero()
 *
 * Set all the values of a double vector to zero.
 *
 * v1 - Pointer to a double vector.
 * size - Size of v1.
 */
void ftmDoubleZero(double *v1, int size)
{
  int i;

  for(i=0;i<size;i++){
    v1[i] = 0.0;
  }
}


/*
 * ftmForwardCheck2D()
 *
 * This function is used for figuring out how FFTW orders frequencies (in 2D).
 *
 * data_fft - A fftw_complex pointer.
 * data - A double pointer to a 2D array.
 * x_size - Data x size.
 * y_size - Data y size.
 */
void ftmForwardCheck2D(fftw_complex *data_fft, double *data, int x_size, int y_size)
{
  fftw_plan fft_forward;

  fft_forward = fftw_plan_dft_r2c_2d(x_size, y_size, data, data_fft, FFTW_ESTIMATE);
  fftw_execute(fft_forward);  
  fftw_destroy_plan(fft_forward);
}
