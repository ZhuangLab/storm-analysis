/*
 * Spline fitting API.
 *
 * Hazen 06/17
 */

#ifndef CUBIC_FIT_H
#define CUBIC_FIT_H

#include "cubic_spline.h"
#include "../sa_library/multi_fit.h"

/* Structures */
typedef struct
{
  int zi;                     /* Location of the spline in z. */

  int x_start;                /* Spline offset in x (0 or 1). */
  int y_start;                /* Spline offset in y (0 or 1). */
  
  double x_delta;             /* Peak x delta (0.0 - 1.0). */
  double y_delta;             /* Peak y delta (0.0 - 1.0). */
  double z_delta;             /* Peak z delta (0.0 - 1.0). */
  
} splinePeak;

  
typedef struct
{
  int fit_type;               /* 2D or 3D spline fit, (S2D or S3D). */

  int spline_size_x;          /* The size of the spline in x (in pixels). */
  int spline_size_y;          /* The size of the spline in y (in pixels). */
  int spline_size_z;          /* The size of the spline in z. */
  
  splineData *spline_data;    /* Spline data structure. */
} splineFit;


/* Functions */
void cfAllocPeaks(peakData *, int);
void cfCalcJH2D(fitData *, double *, double *);
void cfCalcJH3D(fitData *, double *, double *);
void cfCalcJH3DALS(fitData *, double *, double *);
void cfCalcPeakShape(fitData *);
void cfCleanup(fitData *);
void cfCopyPeak(peakData *, peakData *);
void cfFreePeaks(peakData *, int);
fitData* cfInitialize(splineData *, double *, double *, double *, double, int, int);
void cfInitialize2D(fitData *);
void cfInitialize3D(fitData *, int);
void cfNewPeaks(fitData *, double *, char *, int);
void cfUpdate2D(fitData *, double *);
void cfUpdate3D(fitData *, double *);
void cfZRangeCheck(fitData *);

#endif
