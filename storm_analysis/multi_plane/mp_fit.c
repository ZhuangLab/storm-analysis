/*
 * Fit multiple, possibly overlapping, cubic splines
 * to image data from multiple planes.
 *
 * Most of the work is done using spliner/cubic_fit.c.
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
void mpFitDataUpdate(mpFit *, double *, int *, int);
mpFit *mpInitialize(double *, double, int, int, int);
void mpInitializeChannel(mpFit *, splineData *, double *, int);
void mpIterate(mpFit *);
void mpNewImage(mpFit *, double *, int);
void mpNewPeaks(mpFit *, double *, int);
void mpSetTransforms(mpFit *, double *, double *, double *, double *);
void mpSetWeights(mpFit *, double *, double *, double *, double *, double *);
void mpUpdateParameter(peakData *, double *, int, int);
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
 * mpFitDataUpdate()
 *
 * Update peak parameters for peaks in other channels.
 */
void mpFitDataUpdate(mpFit *mp_fit, double *delta, int *good, int pn)
{
  int i,mx,my,xi,yi;
  double t,xoff,yoff;
  double *ch0_params, *params;
  peakData *peak;

  xoff = mp_fit->fit_data[0]->xoff;
  yoff = mp_fit->fit_data[0]->yoff;
  ch0_params = mp_fit->fit_data[0]->fit[pn].params;

  for(i=1;i<mp_fit->n_channels;i++){
    params = mp_fit->fit_data[i]->fit[pn].params;
    
    /* Height and z parameters are the same. */
    params[HEIGHT] = ch0_params[HEIGHT];
    params[ZCENTER] = ch0_params[ZCENTER];

    /* 
     * X and Y need to be mapped first.
     * 
     * Note: The meaning of x and y is transposed here compared to in the
     *       mapping.
     *
     * Note: Spliner uses the upper left corner as 0,0 so we need to adjust
     *       to the center, transform, then adjust back. This is particularly
     *       important if one channel is inverted relative to another.
     */
    t = mp_fit->yt_0toN[i*3];
    t += mp_fit->yt_0toN[i*3+1] * (ch0_params[YCENTER]+yoff);
    t += mp_fit->yt_0toN[i*3+2] * (ch0_params[XCENTER]+xoff);
    params[XCENTER] = t-xoff;

    t = mp_fit->xt_0toN[i*3];
    t += mp_fit->xt_0toN[i*3+1] * (ch0_params[YCENTER]+yoff);
    t += mp_fit->xt_0toN[i*3+2] * (ch0_params[XCENTER]+xoff);
    params[YCENTER] = t-yoff;
    
    /* 
     * Background is updated in the normal way, but only if we
     * have a valid value for delta, otherwise we don't change
     * it.
     *
     * This is probably not the best approach? Should re-try 
     * the fit with only the height and background parameters?
     */
    peak = &mp_fit->fit_data[i]->fit[pn];
    if(good[i]){
      mpUpdateParameter(peak, delta, i, BACKGROUND);
    }

    /* Update peak (integer) location with hysteresis. */
    if(fabs(peak->params[XCENTER] - (double)peak->xi - 0.5) > HYSTERESIS){
      peak->xi = (int)peak->params[XCENTER];
    }
    if(fabs(peak->params[YCENTER] - (double)peak->yi - 0.5) > HYSTERESIS){
      peak->yi = (int)peak->params[YCENTER];
    }

    /*
     * Check that the peak hasn't moved to close to the 
     * edge of the image. Flag the peak as bad if it has.
     */
    xi = peak->xi;
    yi = peak->yi;
    mx = mp_fit->fit_data[i]->image_size_x - peak->size_x;
    my = mp_fit->fit_data[i]->image_size_y - peak->size_y;
    if((xi < 0)||(xi >= mx)||(yi < 0)||(yi >= my)){
      peak->status = ERROR;
      mp_fit->fit_data[i]->n_margin++;
      if(TESTING){
	printf("object outside margins, %d, %d, %d, %d\n", i, peak->index, xi, yi);
      }
    }

    /* 
     * Update peak (integer) z location.
     */
    ((splinePeak *)peak->peak_model)->zi = (int)(peak->params[ZCENTER]);
  }
  
}

/*
 * mpInitialize()
 *
 * Create and return the mpFit structure to use for fitting.
 */
mpFit *mpInitialize(double *clamp, double tolerance, int n_channels, int im_size_x, int im_size_y)
{
  int i;
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

  for(i=0;i<NFITTING;i++){
    mp_fit->clamp_start[i] = clamp[i];
  }

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
  int converged,has_error,i,j;
  int *good;
  double *ch0_delta;
  double *deltas;

  good = (int *)malloc(sizeof(int)*mp_fit->n_channels);
  ch0_delta = (double *)malloc(sizeof(double)*NPEAKPAR);
  deltas = (double *)malloc(sizeof(double)*mp_fit->n_channels*NPEAKPAR);

  for(i=0;i<NPEAKPAR;i++){
    ch0_delta[i] = 0.0;
  }
  for(i=0;i<(mp_fit->n_channels*NPEAKPAR);i++){
    deltas[i] = 0.0;
  }
  
  /* Iterate over localizations. */
  for(i=0;i<mp_fit->nfit;i++){

    /* Skip if this peak is CONVERGED or ERROR. */
    if(mp_fit->fit_data[0]->fit[i].status != RUNNING){
      continue;
    }
    
    /* Calculate updates in each channel. */
    for(j=0;j<mp_fit->n_channels;j++){
      mpUpdateSpline3D(mp_fit->fit_data[j],
		       &mp_fit->fit_data[j]->fit[i],
		       &deltas[NPEAKPAR*j],
		       &good[j]);
    }

    /*
     * Make sure no garbage has snuck into the deltas array.
     */
    if (TESTING){
      for(j=0;j<mp_fit->n_channels;j++){
	if((deltas[NPEAKPAR*j+XWIDTH] != 0.0) || (deltas[NPEAKPAR*j+YWIDTH] != 0.0)){
	  printf("Non-zero x/y width delta detected!\n");
	}
      }
    }

    /* Subtract peaks from each channel. */
    for(j=0;j<mp_fit->n_channels;j++){
      cfSubtractPeak(mp_fit->fit_data[j], &mp_fit->fit_data[j]->fit[i]);
    }

    /* Calculate how to update channel 0 peak. */
    if(mpWeightedDelta(mp_fit, &mp_fit->fit_data[0]->fit[i], deltas, ch0_delta, good)){
      
      /* 
       * Update channel 0 peak. 
       */
      cfFitDataUpdate(mp_fit->fit_data[0], &mp_fit->fit_data[0]->fit[i], ch0_delta);

      if (TESTING){
	if((ch0_delta[XWIDTH] != 0.0) || (ch0_delta[YWIDTH] != 0.0)){
	  printf("Non-zero channel 0 x/y width delta detected!\n");
	}	
      }
      
      /*
	cfFitDataUpdate(mp_fit->fit_data[0], &mp_fit->fit_data[0]->fit[i], deltas);
	cfFitDataUpdate(mp_fit->fit_data[1], &mp_fit->fit_data[1]->fit[i], &deltas[NPEAKPAR]);
      */
      /*
       * Mark all bad if channel 0 peak is bad.
       */
      if(mp_fit->fit_data[0]->fit[i].status == ERROR){
	for(j=1;j<mp_fit->n_channels;j++){
	  mp_fit->fit_data[j]->fit[i].status = ERROR;
	}
      }      
      
      /*
       * Update peaks in the other channels, check that the peaks 
       * are still in the image, etc. 
       */
      mpFitDataUpdate(mp_fit, deltas, good, i);

      /*
       * Mark all bad if any are bad in the other channels.
       */
      has_error = 0;
      for(j=1;j<mp_fit->n_channels;j++){
	if(mp_fit->fit_data[j]->fit[i].status == ERROR){
	  has_error = 1;
	}
      }

      if(has_error){
	for(j=0;j<mp_fit->n_channels;j++){
	  mp_fit->fit_data[j]->fit[i].status = ERROR;
	}
      }
      
      /* Add peaks. */
      if(mp_fit->fit_data[0]->fit[i].status != ERROR){
	for(j=0;j<mp_fit->n_channels;j++){
	  cfAddPeak(mp_fit->fit_data[j], &mp_fit->fit_data[j]->fit[i]);
	}
      }
    }
    else{
      for(j=0;j<mp_fit->n_channels;j++){
	mp_fit->fit_data[j]->fit[i].status = ERROR;
      }
    }
  }

  /* 
   * Update fitting error. 
   *
   * This function also flags bad peaks, so we need to check
   * again for bad peaks after we run it.
   */
  for(i=0;i<mp_fit->nfit;i++){
    for(j=0;j<mp_fit->n_channels;j++){
      mFitCalcErr(mp_fit->fit_data[j],
		  &mp_fit->fit_data[j]->fit[i]);
    }
  }

  /* 
   * If any peak is in a group is in an error state mark all
   * members as being in an error state.
   *
   * If any peak has converged in a group has converged mark them all 
   * as converged. But background may still be incorrect?
   *
   * In the event of CONVERGED and ERROR, ERROR gets priority.
   */
  for(i=0;i<mp_fit->nfit;i++){

    converged = 0;
    has_error = 0;
    for(j=0;j<mp_fit->n_channels;j++){
      if(mp_fit->fit_data[j]->fit[i].status == CONVERGED){
	converged = 1;
      }
      if(mp_fit->fit_data[j]->fit[i].status == ERROR){
	has_error = 1;
      }
    }
    
    if(converged){
      for(j=0;j<mp_fit->n_channels;j++){
	mp_fit->fit_data[j]->fit[i].status = CONVERGED;
      }
    }

    if(has_error){
      for(j=0;j<mp_fit->n_channels;j++){
	mp_fit->fit_data[j]->fit[i].status = ERROR;
      }
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
 * mpUpdateParameter()
 *
 * Update a single parameter of a peak.
 */
void mpUpdateParameter(peakData *peak, double *delta, int channel, int param_index)
{
  double d;
  
  d = delta[channel*NPEAKPAR+param_index];
  if (d != 0.0){

    /* update sign & clamp if the solution appears to be oscillating. */
    if (peak->sign[param_index] != 0){
      if ((peak->sign[param_index] == 1) && (d < 0.0)){
	peak->clamp[param_index] *= 0.5;
      }
      else if ((peak->sign[param_index] == -1) && (d > 0.0)){
	peak->clamp[param_index] *= 0.5;
      }
    }
    if (d > 0.0){
      peak->sign[param_index] = 1;
    }
    else {
      peak->sign[param_index] = -1;
    }
    
    peak->params[param_index] -= d/(1.0 + fabs(d)/peak->clamp[param_index]);
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

  /*
   * X parameters depends on the mapping.
   *
   * Note: The meaning of x and y is transposed here compared to in the
   *       mapping. This is also true for the y parameter below.
   */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(good[i]){
      delta = mp_fit->yt_Nto0[i*3+1] * deltas[NPEAKPAR*i+YCENTER];
      delta += mp_fit->yt_Nto0[i*3+2] * deltas[NPEAKPAR*i+XCENTER];
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
      delta = mp_fit->xt_Nto0[i*3+1] * deltas[NPEAKPAR*i+YCENTER];
      delta += mp_fit->xt_Nto0[i*3+2] * deltas[NPEAKPAR*i+XCENTER];
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
