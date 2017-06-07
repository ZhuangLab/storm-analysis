/*
 * Fit multiple, possibly overlapping, cubic splines
 * to image data from multiple planes.
 *
 * Most of the work is done using the spliner/cubic_fit.c.
 *
 * The expectation is that there will be n_channels copies of
 * each input peak, organized by channel, so for example
 * if there are 3 peaks and 2 channels the peak array would 
 * be [peak1_c1, peak2_c1, peak3_c1, peak1_c2, peak2_c2, peak2_c3].
 * This analysis will then keep the groups of peaks in sync,
 * i.e. peak1_c1 and peak1_c2 will have the same peak status
 * (RUNNING, CONVERGED, ERROR), height and z value. And their
 * x, y coordinates will be the same after affine transformation.
 * BACKGROUND is the only parameter that is independent.
 *
 * Hazen 06/17
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../spliner/cubic_fit.h"


typedef struct
{
  int im_size_x;                /* Image size in x (the fast axis). */
  int im_size_y;                /* Image size in y (the slow axis). */

  int n_channels;               /* The number of different channels / image planes. */
  int nfit;                     /* The number of peaks to fit per channel. The total 
				   number of peaks is n_channels * nfit. */

  double tolerance;             /* Fit tolerance. */
  double clamp_start[NFITTING]; /* Starting value for the peak clamp values. */

  double *xt_0toN;              /* Transform x coordinate from channel 0 to channel N. */
  double *yt_0toN;              /* Transform y coordinate from channel 0 to channel N. */
  double *xt_Nto0;              /* Transform x coordinate from channel N to channel 0. */
  double *yt_Nto0;              /* Transform y coordinate from channel N to channel 0. */
  
  fitData **fit_data;           /* Array of pointers to fitData structures. */
} mpFit;


void mpCleanup(mpFit *);
void mpGetFitImage(mpFit *, double *, int);
void mpGetResults(mpFit *, double *);
int mpGetUnconverged(mpFit *);
mpFit *mpInitialize(double *, double, int, int, int);
void mpInitializeChannel(mpFit *, splineData *, double *, int);
void mpIterate(mpFit *);
void mpNewImage(mpFit *, double *, int);
void mpNewPeaks(mpFit *, double *, int);
void mpSetTransforms(mpFit *, double *, double *, double *, double *);
void mpUpdateSpline3D(fitData *, peakData *, double *, int *);


/* LAPACK Functions */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);


/*
 * mpCleanup()
 *
 * Clean up at the end of the analysis.
 */
void mpCleanup(mpFit *mp_fit)
{
  int i;

  /* Call cubic spline cleanup. */
  for(i=0;i<mp_fit->n_channels;i++){
    cfCleanup(mp_fit->fit_data[i]);
  }

  /* Free affine transform arrays. */
  free(mp_fit->xt_0toN);
  free(mp_fit->yt_0toN);
  free(mp_fit->xt_Nto0);
  free(mp_fit->yt_Nto0);

  free(mp_fit->fit_data);
  free(mp_fit);
}

/*
 * mpGetFitImage()
 *
 * Return an image created from the current best fit peaks.
 */
void mpGetFitImage(mpFit *mp_fit, double *fit_image, int channel)
{
  int i;
  fitData *fdata;

  fdata = mp_fit->fit_data[channel];
  
  for(i=0;i<(fdata->image_size_x * fdata->image_size_y);i++){
    fit_image[i] = fdata->f_data[i];
  }
}

/*
 * mpGetResults()
 *
 * Return the current fitting results.
 */
void mpGetResults(mpFit *mp_fit, double *peak_params)
{
  int i,j;

  for(i=0;i<mp_fit->n_channels;i++){
    j = i*mp_fit->nfit;
    mFitGetResults(mp_fit->fit_data[i], &peak_params[j]);
  }  
}

/*
 * mpGetUnconverged()
 *
 * Return how many peak fits have not converged. 
 */
int mpGetUnconverged(mpFit *mp_fit)
{
  /*
   * We only need to check the peaks for channel 0 as the
   * peaks for other channels will have the same state.
   */
  return mFitGetUnconverged(mp_fit->fit_data[0]);
}

/*
 * mpInitialize()
 *
 * Create and return the mpFit structure to use for fitting.
 */
mpFit *mpInitialize(double *clamp, double tolerance, int n_channels, int im_size_x, int im_size_y)
{
  mpFit *mp_fit;

  mp_fit = (mpFit *)malloc(sizeof(mpFit));

  mp_fit->im_size_x = im_size_x;
  mp_fit->im_size_y = im_size_y;
  mp_fit->n_channels = n_channels;
  mp_fit->tolerance = tolerance;

  mp_fit->xt_0toN = (double *)malloc(3*n_channels*sizeof(double));
  mp_fit->yt_0toN = (double *)malloc(3*n_channels*sizeof(double));
  mp_fit->xt_Nto0 = (double *)malloc(3*n_channels*sizeof(double));
  mp_fit->yt_Nto0 = (double *)malloc(3*n_channels*sizeof(double));

  mp_fit->fit_data = (fitData **)malloc(n_channels*sizeof(fitData*));

  return mp_fit;
}

/*
 * mpInitializeChannel()
 *
 * Initialize a single channel / plane.
 */
void mpInitializeChannel(mpFit *mp_fit, splineData *spline_data, double *variance, int channel)
{
  mp_fit->fit_data[channel] = cfInitialize(spline_data,
					   variance,
					   mp_fit->clamp_start,
					   mp_fit->tolerance,
					   mp_fit->im_size_x,
					   mp_fit->im_size_y);					   
}

/*
 * mpIterate()
 *
 * Perform a single cycle of fitting.
 */
void mpIterate(mpFit *mp_fit)
{
}

/*
 * mpNewImage()
 *
 * Copy in a new image to fit. Note that this should be once
 * per channel with the image for the channel.
 */
void mpNewImage(mpFit *mp_fit, double *new_image, int channel)
{
  mFitNewImage(mp_fit->fit_data[channel], new_image);
}

/*
 * mpNewPeaks()
 *
 * New peaks to fit.
 *
 * n_peaks is the number of peaks per channel.
 */
void mpNewPeaks(mpFit *mp_fit, double *peak_params, int n_peaks)
{
  int i,j;
  
  mp_fit->nfit = n_peaks;
  
  for(i=0;i<mp_fit->n_channels;i++){
    j = i*mp_fit->nfit;
    cfNewPeaks(mp_fit->fit_data[i], &peak_params[j], n_peaks);
  }
}

/*
 * mpSetTransforms()
 *
 * Set affine transform arrays that describe how to change
 * the coordinates between channels.
 */
void mpSetTransforms(mpFit *mp_fit, double *xt_0toN, double *yt_0toN, double *xt_Nto0, double *yt_Nto0)
{
  int i,m;

  m = mp_fit->n_channels*3;
  
  for(i=0;i<m;i++){
    mp_fit->xt_0toN[i] = xt_0toN[i];
    mp_fit->yt_0toN[i] = yt_0toN[i];
    mp_fit->xt_Nto0[i] = xt_Nto0[i];
    mp_fit->yt_Nto0[i] = yt_Nto0[i];
  }
}

/*
 * mpUpdateSpline3D()
 *
 * Determine delta for updating a single peak in a single plane.
 */
void mpUpdateSpline3D(fitData *fit_data, peakData *peak, double *delta, int *good)
{
}
