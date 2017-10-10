/*
 * Structures and functions in pupil_function.c
 *
 * Hazen 10/17
 */

typedef struct pupilData
{
  int size;      /* The size of pupil function in X/Y in pixels. */

  double *kx;    /* X translation multiplier. */
  double *ky;    /* Y translation multiplier. */
  double *kz;    /* Z translation multiplier. */

  fftw_complex *pf;  /* Pupil function. */
  fftw_complex *ws;  /* Working storage. */

  fftw_complex *fftw_pf;
  fftw_complex *fftw_psf;

  fftw_plan fft_backward;
} pupilData;

void pfCleanup(pupilData *);
void pfGetPSF(pupilData *, double *);
void pfGetPSFdx(pupilData *, double *);
pupilData *pfInitialize(double *, double *, double *, int);
void pfSetPF(pupilData *, double *, double *);
