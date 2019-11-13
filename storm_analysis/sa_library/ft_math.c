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
 * ftmBackward()
 *
 * Perform inverse fourier transform.
 *
 * fft_backward - A FFTW plan.
 * fftw_vec - Pointer to memory used by fft_backward plan.
 * vec - Pointer to vector to do inverse FFT of.
 * size - Size of vec.
 */
void ftmBackward(fftw_plan fft_backward, fftw_complex *fftw_vec, fftw_complex *vec, int size)
{
  ftmComplexCopy(vec, fftw_vec, size);
  fftw_execute(fft_backward);
}


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
  if (d1 != s1){
    memcpy(d1, s1, (sizeof(fftw_complex) * size));
  }
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
  if (d1 != s1){
    memcpy(d1, s1, (sizeof(double) * size));
  }
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
 * ftmForward()
 *
 * Perform forward fourier transform.
 *
 * fft_forward - A FFTW plan.
 * fftw_vec - Pointer to memory used by fft_forward plan.
 * vec - Pointer to vector to do FFT of.
 * size - Size of vec.
 */
void ftmForward(fftw_plan fft_forward, double *fftw_vec, double *vec, int size)
{
  ftmDoubleCopy(vec, fftw_vec, size);
  fftw_execute(fft_forward);
}
