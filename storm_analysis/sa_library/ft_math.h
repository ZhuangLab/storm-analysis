/* 
 * Header file for fourier transform and related math.
 *
 * Hazen 11/19
 */

#ifndef FT_MATH_H
#define FT_MATH_H

void ftmBackward(fftw_plan, fftw_complex *, fftw_complex *, int);
void ftmComplexCopy(fftw_complex *, fftw_complex *, int);
void ftmComplexCopyNormalize(fftw_complex *, fftw_complex *, double, int);
void ftmComplexMultiply(fftw_complex *, fftw_complex *, fftw_complex *, int, int);
void ftmComplexMultiplyAccum(fftw_complex *, fftw_complex *, fftw_complex *, int, int);
void ftmComplexZero(fftw_complex *, int);
void ftmDoubleCopy(double *, double *, int);
void ftmDoubleCopyNormalize(double *, double *, double, int);
/* void ftmDoubleNormalize(double *, double, int); */
void ftmForward(fftw_plan, double *, double *, int);

#endif
