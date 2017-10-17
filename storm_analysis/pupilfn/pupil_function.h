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

  double *kx_c;  /* These are all working storage. */
  double *kx_r;
  double *ky_c;
  double *ky_r;
  double *kz_c;
  double *kz_r; 
  
  fftw_complex *pf;  /* Pupil function. */
  fftw_complex *ws;  /* Working storage. */

  fftw_complex *fftw_pf;
  fftw_complex *fftw_psf;

  fftw_plan fft_backward;
} pupilData;

void pfnCleanup(pupilData *);
void pfnGetPSF(pupilData *, double *, double *);
void pfnGetPSFdx(pupilData *, double *, double *);
void pfnGetPSFdy(pupilData *, double *, double *);
void pfnGetPSFdz(pupilData *, double *, double *);
int pfnGetSize(pupilData *);
pupilData *pfnInitialize(double *, double *, double *, int);
void pfnSetPF(pupilData *, double *, double *);
void pfnTranslate(pupilData *, double, double, double);
