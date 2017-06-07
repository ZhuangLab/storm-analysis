/*
 * Spline fitting API.
 *
 * Hazen 06/17
 */

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
  
  double *peak_values;        /* The peak shape. */
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
void cfAddPeak(fitData *, peakData *);
void cfCleanup(fitData *);
void cfFitDataUpdate(fitData *, peakData *, double *);
fitData* cfInitialize(splineData *, double *, double *, double, int, int);
void cfIterateSpline(fitData *);
void cfNewPeaks(fitData *, double *, int);
void cfSubtractPeak(fitData *, peakData *);
void cfUpdateSpline2D(fitData *, peakData *);
void cfUpdateSpline3D(fitData *, peakData *);
