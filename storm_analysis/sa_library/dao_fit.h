/*
 * 3D-DAOSTORM fitting API.
 *
 * Hazen 01/18
 */

#ifndef DAO_FIT_H
#define DAO_FIT_H

#include "multi_fit.h"

/* Structures */
typedef struct daoPeak
{
  double wx_term;
  double wy_term;
  
  double *xt;
  double *ext;
  double *yt;
  double *eyt;
} daoPeak;

typedef struct daoFit
{
  int roi_size;                 /* The size of the fitting ROI. */
  int zfit;                     /* This is a flag for the 'Z' fitting model. */
  
  double min_z;                 /* minimum z value. */
  double max_z;                 /* maximum z value. */  

  double width_min;             /* Minimum allowed width for Gaussian fitting. */
  double width_max;             /* Maximum allowed width for Gaussian fitting. */
  
  double wx_z_params[5];        /* x width versus z parameters. */
  double wy_z_params[5];        /* y width versus z parameters. */
} daoFit;


/* Functions */
void daoAllocPeaks(peakData *, int);
void daoCalcJH2DFixed(fitData *, double *, double *);
void daoCalcJH2DFixedALS(fitData *, double *, double *);
void daoCalcJH2DFixedLS(fitData *, double *, double *);
void daoCalcJH2DFixedDWLS(fitData *, double *, double *);
void daoCalcJH2DFixedFWLS(fitData *, double *, double *);
void daoCalcJH2D(fitData *, double *, double *);
void daoCalcJH2DALS(fitData *, double *, double *);
void daoCalcJH2DLS(fitData *, double *, double *);
void daoCalcJH2DDWLS(fitData *, double *, double *);
void daoCalcJH2DFWLS(fitData *, double *, double *);
void daoCalcJH3D(fitData *, double *, double *);
void daoCalcJHZ(fitData *, double *, double *);
void daoCalcPeakShape(fitData *);
void daoCalcWidthsFromZ(fitData *, peakData *);
int daoCheck2DFixed(fitData *);
int daoCheck2D(fitData *);
int daoCheck3D(fitData *);
int daoCheckZ(fitData *);
void daoCleanup(fitData *);
void daoCopyPeak(fitData *, peakData *, peakData *);
void daoFreePeaks(peakData *, int);
fitData* daoInitialize(double *, double *, double, int, int, int);
void daoInitialize2DFixed(fitData *);
void daoInitialize2DFixedALS(fitData *);
void daoInitialize2DFixedLS(fitData *);
void daoInitialize2DFixedDWLS(fitData *);
void daoInitialize2DFixedFWLS(fitData *);
void daoInitialize2D(fitData *, double, double);
void daoInitialize2DALS(fitData *, double, double);
void daoInitialize2DLS(fitData *, double, double);
void daoInitialize2DDWLS(fitData *, double, double);
void daoInitialize2DFWLS(fitData *, double, double);
void daoInitialize3D(fitData *, double, double);
void daoInitializeZ(fitData *, double *, double *, double, double);
void daoNewPeaks(fitData *, double *, char *, int);
void daoUpdate2DFixed(fitData *, double *);
void daoUpdate2D(fitData *, double *);
void daoUpdate3D(fitData *, double *);
void daoUpdateZ(fitData *, double *);

#endif
