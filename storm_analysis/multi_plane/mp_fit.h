/*
 * API for the multi-plane fitting library.
 *
 * Hazen 01/18
 */

#ifndef MP_FIT_H
#define MP_FIT_H

#include "../sa_library/multi_fit.h"

typedef struct mpFit
{
  int im_size_x;                /* Image size in x (the fast axis). */
  int im_size_y;                /* Image size in y (the slow axis). */
  
  int independent_heights;      /* Flag for independent peak heights. */
  
  int n_channels;               /* The number of different channels / image planes. */
  int n_weights;                /* The number of (z) weight values. */

  double tolerance;             /* Fit tolerance. */
    
  double w_z_offset;            /* Offset value to convert peak z to a weight index. */
  double w_z_scale;             /* Scale value to convert peak z to a weight index. */
    
  double zmin;                  /* Minimum allowed z value, units are fitter dependent. */
  double zmax;                  /* Maximum allowed z value, units are fitter dependent. */

  double *xt_0toN;              /* Transform x coordinate from channel 0 to channel N. */
  double *yt_0toN;              /* Transform y coordinate from channel 0 to channel N. */
  double *xt_Nto0;              /* Transform x coordinate from channel N to channel 0. */
  double *yt_Nto0;              /* Transform y coordinate from channel N to channel 0. */

  double *w_bg;                 /* Per channel z dependent weighting for the background parameter. */
  double *w_h;                  /* Per channel z dependent weighting for the height parameter. */
  double *w_x;                  /* Per channel z dependent weighting for the x parameter. */
  double *w_y;                  /* Per channel z dependent weighting for the y parameter. */
  double *w_z;                  /* Per channel z dependent weighting for the z parameter. */
  double *heights;              /* Per channel heights for parameter weighting. */

  double **jacobian;            /* Storage for the jacobian calculations. */
  double **w_jacobian;          /* Storage for copies of the jacobians. */
  double **hessian;             /* Storage for the hessian calculations. */
  double **w_hessian;           /* Storage for copies of the jacobians. */
  
  fitData **fit_data;           /* Array of pointers to fitData structures. */

  void (*fn_cleanup)(struct fitData *);                         /* Function for cleaning up fitting for a particular channel. */
  void (*fn_newpeaks)(struct fitData *, double *, char *, int); /* Function for adding new peaks for a particular channel. */
  void (*fn_peak_xi_yi)(struct peakData *);                     /* Function for updating peak xi, yi values. */
  void (*fn_update)(struct mpFit *);                            /* Function for updating the parameters of the working peaks. */
  void (*fn_zrange)(struct fitData *);                          /* Function for enforcing the z range. */
  
} mpFit;


void mpCheckError(mpFit *, int);
void mpCleanup(mpFit *);
void mpCopyFromWorking(mpFit *, int, int);
mpFit *mpInitialize(double, int, int, int, int);
void mpIterateLM(mpFit *);
void mpResetWorkingPeaks(mpFit *, int);
void mpSetTransforms(mpFit *, double *, double *, double *, double *);
void mpSetWeights(mpFit *, double *, double *, double *, double *, double *, int);
void mpSetWeightsIndexing(mpFit *, double, double);
void mpUpdate(mpFit *);
void mpWeightIndex(mpFit *, double *, int *);
double mpWeightInterpolate(double *, double, int, int, int);

#endif
