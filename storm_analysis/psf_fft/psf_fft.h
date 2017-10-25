/*
 * Header file for psf_fft.c.
 *
 * Hazen 10/17
 */

#ifndef PSF_FFT_H
#define PSF_FFT_H

#include <fftw3.h>

typedef struct psfFFT
{  
  int x_size;              /* PSF size in x. */
  int y_size;              /* PSF size in y. */
  int z_size;              /* PSF size in z. */
  int psf_size;            /* x_size * y_size * z_size. */

  int fft_x_size;          /* FFT size in x (fast dimension), x_size/2 + 1 */
  int fft_size;            /* (x_size/2 + 1) * y_size * z_size. */
  
  double *kx_c;            /* These are all working storage. */
  double *kx_r;
  double *ky_c;
  double *ky_r;
  double *kz_c;
  double *kz_r;

  double *fftw_real;      /* Real space vector for FFTW. */
  fftw_complex *fftw_fft; /* FFT space vector for FFTW. */
  
  fftw_complex *ws;       /* Working storage. */
  fftw_complex *psf;      /* Fourier transform of the PSF. */

  fftw_plan fft_backward;
} psfFFT;

void pFTCalcShiftVector(double *, double *, double, int);
void pFTCalcShiftVectorDerivative(double *, int);
void pFTCleanup(psfFFT *);
void pFTGetPSF(psfFFT *, double *);
void pFTGetPSFdx(psfFFT *, double *);
void pFTGetPSFdy(psfFFT *, double *);
void pFTGetPSFdz(psfFFT *, double *);
int pFTGetXSize(psfFFT *);
int pFTGetYSize(psfFFT *);
int pFTGetZSize(psfFFT *);
psfFFT *pFTInitialize(double *, int, int, int);
void pFTTranslate(psfFFT *, double, double, double);

#endif
