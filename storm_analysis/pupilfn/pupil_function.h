/*
 * Structures and functions in pupil_function.c
 *
 * Hazen 10/17
 */

#ifndef PUPIL_FUNCTION_H
#define PUPIL_FUNCTION_H

#include <fftw3.h>

typedef struct pupilData
{
  int size;      /* The size of pupil function in X/Y in pixels. */

  double *kx;    /* X translation multiplier. */
  double *ky;    /* Y translation multiplier. */
  double *kz;    /* Z translation multiplier. */

  double *kx_c;  /* These are all working storage. */
  double *kx_r;
  double *ky_c;
  double *ky_r;
  double *kz_c;
  double *kz_r;

  /* 
   * These are all for a vectorial PSF. At present they are
   * not used for PF fitting, they are just used by the
   * simulator/pupil_math.py for creating a vectorial PSF.
   */
  double *px_ex;
  double *px_ey;
  double *py_ex;
  double *py_ey;
  double *pz_ex;
  double *pz_ey;
  
  fftw_complex *pf;  /* Pupil function. */
  fftw_complex *ws;  /* Working storage. */

  fftw_complex *fftw_pf;
  fftw_complex *fftw_psf;

  fftw_plan fft_backward;
} pupilData;

void pfnCleanup(pupilData *);
void pfnGetPSF(pupilData *, double *, double *);
void pfnGetPSFIntensity(pupilData *, double *);
void pfnGetPSFdx(pupilData *, double *, double *);
void pfnGetPSFdy(pupilData *, double *, double *);
void pfnGetPSFdz(pupilData *, double *, double *);
void pfnGETPNEN(pupilData *, double *, double *, double *);
void pfnGetPXEX(pupilData *, double *, double *);
void pfnGetPXEY(pupilData *, double *, double *);
void pfnGetPYEX(pupilData *, double *, double *);
void pfnGetPYEY(pupilData *, double *, double *);
void pfnGetPZEX(pupilData *, double *, double *);
void pfnGetPZEY(pupilData *, double *, double *);
int pfnGetSize(pupilData *);
pupilData *pfnInitialize(double *, double *, double *, int);
void pfnSetPF(pupilData *, double *, double *);
void pfnSetPNEN(pupilData *, double *, double *, double *, double *, double *, double *);
void pfnTranslate(pupilData *, double, double, double);
void pfnTranslateZ(pupilData *, double);

#endif

