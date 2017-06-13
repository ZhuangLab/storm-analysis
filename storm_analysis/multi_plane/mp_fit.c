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

  double *w_bg;                 /* Per channel z dependent weighting for the background parameter. */
  double *w_h;                  /* Per channel z dependent weighting for the height parameter. */
  double *w_x;                  /* Per channel z dependent weighting for the x parameter. */
  double *w_y;                  /* Per channel z dependent weighting for the y parameter. */
  double *w_z;                  /* Per channel z dependent weighting for the z parameter. */
  
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
void mpSetWeights(mpFit *, double *, double *, double *, double *, double *);
void mpUpdateSpline3D(fitData *, peakData *, double *, int *);
int mpWeightedDelta(mpFit *, peakData *, double *, double *, int *);

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

  /* Free weight arrays. */
  free(mp_fit->w_bg);
  free(mp_fit->w_h);
  free(mp_fit->w_x);
  free(mp_fit->w_y);
  free(mp_fit->w_z);

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
    j = i*mp_fit->nfit*NPEAKPAR;
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
 * Perform a single cycle of fitting for each localization.
 */
void mpIterate(mpFit *mp_fit)
{
  int i,j,np;
  int *good;
  double *ch0_delta;
  double *deltas;

  good = (int *)malloc(sizeof(int)*mp_fit->n_channels);
  ch0_delta = (double *)malloc(sizeof(double)*NPEAKPAR);
  deltas = (double *)malloc(sizeof(double)*mp_fit->n_channels*NPEAKPAR);

  /* Iterate over localizations. */
  for(i=0;i<mp_fit->nfit;i++){
    
    /* Calculate updates in each channel. */
    for(j=0;j<mp_fit->n_channels;j++){
      mpUpdateSpline3D(mp_fit->fit_data[j],
		       &mp_fit->fit_data[j]->fit[i],
		       &delta[NPEAKPAR*j],
		       &good[j])
    }

    /* Subtract peaks from each channel. */
    for(j=0;j<mp_fit->n_channels;j++){
      cfSubtractPeak(mp_fit->fit_data[j],
		     &mp_fit->fit_data[j]->fit[i])
    }

    /* Calculate how to update channel 0 peak. */
    mpWeightedDelta(mp_fit,
		    &mp_fit->fit_data[0]->fit[i],
		    deltas,
		    ch0_delta,
		    good)

    /* Update peaks. */

    /* Check that the peaks are still in the image, etc. */

    /* Add peaks. */
    for(j=0;j<mp_fit->n_channels;j++){
      cfAddPeak(mp_fit->fit_data[j],
		&mp_fit->fit_data[j]->fit[i])
    }
  }

  /* Update fitting error. */
  for(i=0;i<mp_fit->nfit;i++){
    for(j=0;j<mp_fit->n_channels;j++){
      mFitCalcErr(mp_fit->fit_data[j],
		  &mp_fit->fit_data[j]->fit[i]);
    }
  }
  
  free(good);
  free(ch0_delta);
  free(deltas);
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
    j = i*mp_fit->nfit*NPEAKPAR;
    cfNewPeaks(mp_fit->fit_data[i], &peak_params[j], n_peaks);
  }
}

/*
 * mpSetTransforms()
 *
 * Set affine transform arrays that describe how to change
 * the coordinates between channels. 
 *
 * These are expected to be by channel, then by coefficient.
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
 * mpSetWeights()
 *
 * Set values to use when averaging the per-channel updates. For
 * now the background parameter is independent for each channel,
 * though this may change, so we set it anyway.
 *
 * These are expected to be indexed by z, then channel, so the z
 * value is the slow axis and the channel is the fast axis.
 *
 * The overall size is the number of channels times the spline
 * size in z.
 *
 * This cannot be called before mpInitializeChannel().
 */
void mpSetWeights(mpFit *mp_fit, double *w_bg, double *w_h, double *w_x, double *w_y, double *w_z)
{
  int i,n,z_size;

  /* Figure out spline z size. */
  z_size = ((splineFit *)mp_fit->fit_data[0]->fit_model)->spline_size_z;
  
  /* Allocate storage. */
  n = mp_fit->n_channels*z_size;
  mp_fit->w_bg = (double *)malloc(sizeof(double)*n);
  mp_fit->w_h = (double *)malloc(sizeof(double)*n);
  mp_fit->w_x = (double *)malloc(sizeof(double)*n);
  mp_fit->w_y = (double *)malloc(sizeof(double)*n);
  mp_fit->w_z = (double *)malloc(sizeof(double)*n);
  
  /* Copy values. */
  for(i=0;i<n;i++){
    mp_fit->w_bg[i] = w_bg[i];
    mp_fit->w_h[i] = w_h[i];
    mp_fit->w_x[i] = w_x[i];
    mp_fit->w_y[i] = w_y[i];
    mp_fit->w_z[i] = w_z[i];
  }
}

/*
 * mpUpdateSpline3D()
 *
 * Determine delta for updating a single peak in a single plane. This is pretty
 * much an exact copy for spliner/cubic_fit.c except that we return the delta
 * and whether or not LAPACKs dposv_() function failed.
 */
void mpUpdateSpline3D(fitData *fit_data, peakData *peak, double *delta, int *good)
{
  /* These are for Lapack */
  int o = 5, nrhs = 1, lda = 5, ldb = 5, info;

  int i,j,k,l,m,n,zi;
  int x_start, y_start;
  double height,fi,t1,t2,xi;
  double jt[5];
  double jacobian[5];
  double hessian[25];  
  splinePeak *spline_peak;
  splineFit *spline_fit;
  
  /* Initializations. */
  spline_peak = (splinePeak *)peak->peak_model;
  spline_fit = (splineFit *)fit_data->fit_model;
      
  x_start = spline_peak->x_start;
  y_start = spline_peak->y_start;
  zi = spline_peak->zi;

  *good = 0;
  
  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }
  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }

  /* Calculate values x, y, z, xx, xy, yy, etc. terms for a 3D spline. */
  computeDelta3D(spline_fit->spline_data, spline_peak->z_delta, spline_peak->y_delta, spline_peak->x_delta);
  
  /*
   * Calculate jacobian and hessian.
   */
  height = peak->params[HEIGHT];
  i = peak->yi * fit_data->image_size_x + peak->xi;
  for(j=0;j<peak->size_y;j++){
    for(k=0;k<peak->size_x;k++){
      l = i + j * fit_data->image_size_x + k;
      fi = fit_data->f_data[l] + fit_data->bg_data[l] / ((double)fit_data->bg_counts[l]);
      xi = fit_data->x_data[l];

      /*
       * The derivative in x and y is multiplied by 0.5 as 
       * this is 1.0/(spline up-sampling, i.e. 2x).
       */
      jt[0] = spline_peak->peak_values[j*peak->size_x + k];
      jt[1] = -0.5*height*dxfAt3D(spline_fit->spline_data,zi,2*j+y_start,2*k+x_start);
      jt[2] = -0.5*height*dyfAt3D(spline_fit->spline_data,zi,2*j+y_start,2*k+x_start);
      jt[3] = height*dzfAt3D(spline_fit->spline_data,zi,2*j+y_start,2*k+x_start);
      jt[4] = 1.0;

      /* Calculate jacobian. */
      t1 = 2.0*(1.0 - xi/fi);
      for(m=0;m<5;m++){
	jacobian[m] += t1*jt[m];
      }
	  
      /* Calculate hessian. */
      t2 = 2.0*xi/(fi*fi);
      for(m=0;m<5;m++){
	for(n=m;n<5;n++){
	  hessian[m*5+n] += t2*jt[m]*jt[n];
	}
      }
    }
  }
      
  /* Use Lapack to solve AX=B to calculate update vector. */
  dposv_( "Lower", &o, &nrhs, hessian, &lda, jacobian, &ldb, &info );

  if(info!=0){
    peak->status = ERROR;
    fit_data->n_dposv++;
    if(TESTING){
      printf("fitting error! %d %d %d\n", peak->index, info, ERROR);
    }
  }
  else{
    
    /* Update params. */
    delta[HEIGHT]     = jacobian[0];
    delta[XCENTER]    = jacobian[1];
    delta[YCENTER]    = jacobian[2];
    delta[ZCENTER]    = jacobian[3];
    delta[BACKGROUND] = jacobian[4];

    *good = 1;
  }
}

/*
 * mpWeightedDelta()
 *
 * Given the delta for each channel, figure out the (hopefully)
 * optimal delta for the localization.
 */
int mpWeightedDelta(mpFit *mp_fit, peakData *peak, double *deltas, double *ch0_delta, int *good)
{
  int i,nc,zi;
  double delta, p_ave, p_total;

  nc = mp_fit->n_channels;
  zi = ((splinePeak *)peak->peak_model)->zi;

  /* Height parameter is a simple weighted average. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(good[i]){
      p_ave += deltas[NPEAKPAR*i+HEIGHT] * mp_fit->w_h[zi*nc+i];
      p_total += mp_fit->w_h[zi*nc+i];
    }
  }

  /*
   * We are assuming that there are no zero / negative weight 
   * values and only doing this check once.
   */
  if(p_total>0.0){
    ch0_delta[HEIGHT] = p_ave/p_total;
  }
  else {
    return 0;
  }

  /* X parameters depends on the mapping. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(good[i]){
      delta = mp_fit->xt_Nto0[i*3+1] * deltas[NPEAKPAR*i+XCENTER];
      delta += mp_fit->xt_Nto0[i*3+2] * deltas[NPEAKPAR*i+YCENTER];
      p_ave += delta * mp_fit->w_x[zi*nc+i];
      p_total += mp_fit->w_x[zi*nc+i];
    }
  }
  ch0_delta[XCENTER] = p_ave/p_total;

  /* Y parameters depends on the mapping. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(good[i]){
      delta = mp_fit->yt_Nto0[i*3+1] * deltas[NPEAKPAR*i+XCENTER];
      delta += mp_fit->yt_Nto0[i*3+2] * deltas[NPEAKPAR*i+YCENTER];
      p_ave += delta * mp_fit->w_y[zi*nc+i];
      p_total += mp_fit->w_y[zi*nc+i];
    }
  }
  ch0_delta[YCENTER] = p_ave/p_total;

  /* Z parameter is a simple weighted average. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(good[i]){
      p_ave += deltas[NPEAKPAR*i+ZCENTER] * mp_fit->w_z[zi*nc+i];
      p_total += mp_fit->w_z[zi*nc+i];
    }
  }
  ch0_delta[ZCENTER] = p_ave/p_total;

  /* Background terms float independently. */
  ch0_delta[BACKGROUND] = deltas[BACKGROUND];

  return 1;
}
