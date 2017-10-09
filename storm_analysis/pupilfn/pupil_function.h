/*
 * Structures and functions in pupil_function.c
 *
 * Hazen 10/17
 */

typedef struct pupilData
{
  int size;      /* The size of pupil function in X/Y in pixels. */

  fftw_complex *pf;  /* Pupil function. */
  fftw_complex *kx;  /* X translation multiplier. */
  fftw_complex *ky;  /* Y translation multiplier. */
  fftw_complex *kz;  /* Z translation multiplier. */

  fftw_complex *fftw_pf;
  fftw_complex *fftw_psf;

  fftw_plan fft_backward;
} pupilData;

void pfCleanup(pupilData *);
void pfGetPSF(pupilData *, double *);
pupilData *pfInitialize(double *, double *, double *, double *, double *, double *, int);
void pfSetPf(pupilData *, double *, double *);
