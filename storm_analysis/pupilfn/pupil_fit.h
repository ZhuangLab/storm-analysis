/*
 * Pupil function fitting API.
 *
 * Hazen 10/17
 */

#ifndef PUPIL_FIT_H
#define PUPIL_FIT_H

#include "pupil_function.h"
#include "../sa_library/multi_fit.h"

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

  double max_z;   /* Maximum allowed z value (in microns). */
  double min_z;   /* Minimum allowed z value (in microns). */
  
  double *dx_c;   /* Temporary storage for x derivative (complex part). */
  double *dx_r;   /* Temporary storage for x derivative (real part). */
  double *dy_c;   /* Temporary storage for y derivative (complex part). */
  double *dy_r;   /* Temporary storage for y derivative (real part). */
  double *dz_c;   /* Temporary storage for z derivative (complex part). */
  double *dz_r;   /* Temporary storage for z derivative (real part). */
  
  pupilData *pupil_data;    /* Pupil function data structure. */
} pupilFit;

void pfitAllocPeaks(peakData *, int);
void pfitCalcJH3D(fitData *, double *, double *);
void pfitCalcPeakShape(fitData *);
void pfitCleanup(fitData *);
void pfitCopyPeak(fitData *, peakData *, peakData *);
void pfitFreePeaks(peakData *, int);
fitData* pfitInitialize(pupilData *, double *, double *, double, int, int);
void pfitNewPeaks(fitData *, double *, char *, int);
void pfitSetZRange(fitData *, double, double);
void pfitUpdate3D(fitData *, double *);
void pfitZRangeCheck(fitData *);

#endif
