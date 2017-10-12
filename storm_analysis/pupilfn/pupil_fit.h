/*
 * Pupil function fitting API.
 *
 * Hazen 10/17
 */

/* Structures */
typedef struct pupilPeak
{
  double dx;      /* Peak delta x (difference between actual and integer position). */
  double dy;      /* Peak delta y (difference between actual and integer position). */
  double dz;      /* Peak delta z (difference between actual and integer position). */
  
  double *psf_c;  /* The complex part of peak shape. */
  double *psf_r;  /* The real part of peak shape. */
} pupilPeak;

  
typedef struct pupilFit
{
  int pupil_size; /* The size in X/Y of the pupil function. */
  
  double *dx_c;   /* Temporary storage for x derivative (complex part). */
  double *dx_r;   /* Temporary storage for x derivative (real part). */
  double *dy_c;   /* Temporary storage for y derivative (complex part). */
  double *dy_r;   /* Temporary storage for y derivative (real part). */
  double *dz_c;   /* Temporary storage for z derivative (complex part). */
  double *dz_r;   /* Temporary storage for z derivative (real part). */
  
  pupilData *pupil_data;    /* Pupil function data structure. */
} pupilFit;

void pfitAddPeak(fitData *);
void pfitCalcJH3D(fitData *, double *, double *);
void pfitCleanup(fitData *);
void pfitCopyPeak(peakData *, peakData *);
fitData* pfitInitialize(pupilData *, double *, double *, double, int, int);
void pfitNewPeaks(fitData *, double *, int);
void pfitSubtractPeak(fitData *);
void pfitUpdate3D(fitData *, double *);
