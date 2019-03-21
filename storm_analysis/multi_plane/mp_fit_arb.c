/*
 * Multi-plane fitting with 'arbitrary' PSFs.
 *
 * Most of the work is done using one of:
 *  1. psf_fft/fft_fit.c
 *  2. pupilfn/pupil_fit.c
 *  3. spliner/cubic_fit.c
 *
 * Hazen 10/17
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mp_fit.h"

#include "../psf_fft/fft_fit.h"
#include "../pupilfn/pupil_fit.h"
#include "../spliner/cubic_fit.h"


void mpArbInitializePSFFFTChannel(mpFit *, psfFFT *, double *, double *, int);
void mpArbInitializePupilFnChannel(mpFit *, pupilData *, double *, double *, double, double, int);
void mpArbInitializeSplineChannel(mpFit *, splineData *, double *, double *, int);
void mpArbNewPeaks(mpFit *, double *, char *, int);
void mpArbUpdateZ(mpFit *);
void mpArbUpdateFixed(mpFit *);
void mpArbUpdateIndependent(mpFit *);


/*
 * mpArbInitializePSFFFTChannel()
 *
 * Initialize a single channel / plane for 3D PSFFFT fitting.
 */
void mpArbInitializePSFFFTChannel(mpFit *mp_fit, psfFFT *psf_fft_data, double *rqe, double *variance, int channel)
{
  int jac_size;

  /* Specify how to add new peaks and how to cleanup. */
  if(channel == 0){
    mp_fit->fn_cleanup = &ftFitCleanup;
    mp_fit->fn_newpeaks = &ftFitNewPeaks;
    mp_fit->fn_peak_xi_yi = &mFitUpdate;
    mp_fit->fn_zrange = &ftFitZRangeCheck;

    /* Specify how to handle heights in each channel. */
    if(mp_fit->independent_heights){
      mp_fit->fn_update = &mpArbUpdateIndependent;
    }
    else{
      mp_fit->fn_update = &mpArbUpdateFixed;
    }
  }
  
  /*
   * Initialize pupil function fitting for this channel / plane.
   */
  mp_fit->fit_data[channel] = ftFitInitialize(psf_fft_data,
					      rqe,
					      variance,
					      mp_fit->tolerance,
					      mp_fit->im_size_x,
					      mp_fit->im_size_y);
  mp_fit->fit_data[channel]->minimum_height = 1.0;
  
  /*
   * Allocate storage for jacobian and hessian calculations.
   */
  jac_size = mp_fit->fit_data[channel]->jac_size;
  mp_fit->jacobian[channel] = (double *)malloc(jac_size*sizeof(double));
  mp_fit->w_jacobian[channel] = (double *)malloc(jac_size*sizeof(double));
  mp_fit->hessian[channel] = (double *)malloc(jac_size*jac_size*sizeof(double));
  mp_fit->w_hessian[channel] = (double *)malloc(jac_size*jac_size*sizeof(double));

}


/*
 * mpArbInitializePupilFnChannel()
 *
 * Initialize a single channel / plane for 3D pupil function fitting.
 */
void mpArbInitializePupilFnChannel(mpFit *mp_fit, pupilData *pupil_data, double *rqe, double *variance, double zmin, double zmax, int channel)
{
  int jac_size;

  /* Specify how to add new peaks and how to cleanup. */
  if(channel == 0){
    mp_fit->fn_cleanup = &pfitCleanup;
    mp_fit->fn_newpeaks = &pfitNewPeaks;
    mp_fit->fn_peak_xi_yi = &mFitUpdate;
    mp_fit->fn_zrange = &pfitZRangeCheck;

    /* Specify how to handle heights in each channel. */
    if(mp_fit->independent_heights){
      mp_fit->fn_update = &mpArbUpdateIndependent;
    }
    else{
      mp_fit->fn_update = &mpArbUpdateFixed;
    }
  }
  
  /*
   * Initialize pupil function fitting for this channel / plane.
   */
  mp_fit->fit_data[channel] = pfitInitialize(pupil_data,
					     rqe,
					     variance,
					     mp_fit->tolerance,
					     mp_fit->im_size_x,
					     mp_fit->im_size_y);
  mp_fit->fit_data[channel]->minimum_height = 1.0;
  
  pfitSetZRange(mp_fit->fit_data[channel], zmin, zmax);
  
  /*
   * Allocate storage for jacobian and hessian calculations.
   */
  jac_size = mp_fit->fit_data[channel]->jac_size;
  mp_fit->jacobian[channel] = (double *)malloc(jac_size*sizeof(double));
  mp_fit->w_jacobian[channel] = (double *)malloc(jac_size*sizeof(double));
  mp_fit->hessian[channel] = (double *)malloc(jac_size*jac_size*sizeof(double));
  mp_fit->w_hessian[channel] = (double *)malloc(jac_size*jac_size*sizeof(double));
}


/*
 * mpArbInitializeSplineChannel()
 *
 * Initialize a single channel / plane for 3D spline fitting.
 */
void mpArbInitializeSplineChannel(mpFit *mp_fit, splineData *spline_data, double *rqe, double *variance, int channel)
{
  int jac_size;
  
  /* Specify how to add new peaks and how to cleanup. */
  if(channel == 0){
    mp_fit->fn_cleanup = &cfCleanup;
    mp_fit->fn_newpeaks = &cfNewPeaks;
    mp_fit->fn_peak_xi_yi = &mFitUpdate;
    mp_fit->fn_zrange = &cfZRangeCheck;

    /* Specify how to handle heights in each channel. */
    if(mp_fit->independent_heights){
      mp_fit->fn_update = &mpArbUpdateIndependent;
    }
    else{
      mp_fit->fn_update = &mpArbUpdateFixed;
    }
  }
  
  /*
   * Initialize spliner fitting for this channel / plane.
   */
  mp_fit->fit_data[channel] = cfInitialize(spline_data,
					   rqe,
					   variance,
					   mp_fit->tolerance,
					   mp_fit->im_size_x,
					   mp_fit->im_size_y);
  mp_fit->fit_data[channel]->minimum_height = 1.0;
  
  cfInitialize3D(mp_fit->fit_data[channel]);
  
  /*
   * Allocate storage for jacobian and hessian calculations.
   */
  jac_size = mp_fit->fit_data[channel]->jac_size;
  mp_fit->jacobian[channel] = (double *)malloc(jac_size*sizeof(double));
  mp_fit->w_jacobian[channel] = (double *)malloc(jac_size*sizeof(double));
  mp_fit->hessian[channel] = (double *)malloc(jac_size*jac_size*sizeof(double));
  mp_fit->w_hessian[channel] = (double *)malloc(jac_size*jac_size*sizeof(double));
}


/*
 * mpArbNewPeaks()
 *
 * New peaks to fit.
 *
 * n_peaks is the number of peaks per channel.
 */
void mpArbNewPeaks(mpFit *mp_fit, double *peak_params, char *p_type, int n_peaks)
{
  int i,j,k;
  int start,stop;
  double height,tx,ty;
  double *mapped_peak_params;
  fitData *fit_data;
  
  if(VERBOSE){
    printf("mpNP %d\n", n_peaks);
  }

  start = mp_fit->fit_data[0]->nfit;
  stop = start + n_peaks;
  
  if(!strcmp(p_type, "finder") || !strcmp(p_type, "testing")){

    /* We'll use this to pass the mapped peak positions. */
    mapped_peak_params = (double *)malloc(sizeof(double)*n_peaks*3);

    /* Map peak positions & pass to each of the fitters. */
    for(i=0;i<mp_fit->n_channels;i++){
      if(i>0){
	for(j=0;j<n_peaks;j++){
	  k = 3*j;
	  tx = peak_params[k];
	  ty = peak_params[k+1];
	  mapped_peak_params[k] = mp_fit->xt_0toN[i*3] + tx*mp_fit->xt_0toN[i*3+1] + ty*mp_fit->xt_0toN[i*3+2];
	  mapped_peak_params[k+1] = mp_fit->yt_0toN[i*3] + tx*mp_fit->yt_0toN[i*3+1] + ty*mp_fit->yt_0toN[i*3+2];
	  mapped_peak_params[k+2] = peak_params[k+2];
	}
	mp_fit->fn_newpeaks(mp_fit->fit_data[i], mapped_peak_params, p_type, n_peaks);
      }
      else{
	mp_fit->fn_newpeaks(mp_fit->fit_data[i], peak_params, p_type, n_peaks);
      }
    }

    /* Correct heights and errors when peaks are not independent. */
    if(!mp_fit->independent_heights){
      for(i=start;i<stop;i++){
	
	/* 
	 * Copy current peaks into working peaks and calculate
	 * average height.
	 */
	height = 0.0;
	for(j=0;j<mp_fit->n_channels;j++){
	  fit_data = mp_fit->fit_data[j];
	  fit_data->fn_copy_peak(fit_data, &fit_data->fit[i], fit_data->working_peak);
	  height += fit_data->working_peak->params[HEIGHT];
	}
	height = height/((double)mp_fit->n_channels);

	/* Subtract current peaks from the image. */
	for(j=0;j<mp_fit->n_channels;j++){
	  fit_data = mp_fit->fit_data[j];
	  if(fit_data->working_peak->status != ERROR){
	    mFitSubtractPeak(fit_data);
	  }
	}

	/* 
	 * Set all peaks to have the same height, add back into fit 
	 * image, calculate their error & copy back from the working peak.
	 *
	 * Note: We don't have to re-calculate the peak shape because
	 *       it has not changed, just the height.
	 */
	for(j=0;j<mp_fit->n_channels;j++){
	  fit_data = mp_fit->fit_data[j];
	  fit_data->working_peak->params[HEIGHT] = height;
	  if(fit_data->working_peak->status != ERROR){
	    mFitAddPeak(fit_data);
	    mFitCalcErr(fit_data);	    
	  }
	  fit_data->fn_copy_peak(fit_data, fit_data->working_peak, &fit_data->fit[i]);
	}
      }
    }

    free(mapped_peak_params);
  }
  else{
    mapped_peak_params = (double *)malloc(sizeof(double)*n_peaks*5);

    /* Map peak positions & pass to each of the fitters. */
    for(i=0;i<mp_fit->n_channels;i++){
      if(i>0){
	for(j=0;j<n_peaks;j++){
	  k = 5*j;
	  tx = peak_params[k];
	  ty = peak_params[k+1];
	  mapped_peak_params[k] = mp_fit->xt_0toN[i*3] + tx*mp_fit->xt_0toN[i*3+1] + ty*mp_fit->xt_0toN[i*3+2];
	  mapped_peak_params[k+1] = mp_fit->yt_0toN[i*3] + tx*mp_fit->yt_0toN[i*3+1] + ty*mp_fit->yt_0toN[i*3+2];
	  mapped_peak_params[k+2] = peak_params[k+2];
	  mapped_peak_params[k+3] = peak_params[k+3];
	  mapped_peak_params[k+4] = peak_params[k+4];
	}
	mp_fit->fn_newpeaks(mp_fit->fit_data[i], mapped_peak_params, p_type, n_peaks);
      }
      else{
	mp_fit->fn_newpeaks(mp_fit->fit_data[i], peak_params, p_type, n_peaks);
      }
    }
    
    free(mapped_peak_params);
  }

  /* 
   * Check for error peaks & synchronize status. This can happen for example
   * because the peak in one channel is outside the image.
   */
  for(i=start;i<stop;i++){
    mpCheckError(mp_fit, i);
  }
}


/*
 * mpArbUpdateZ()
 *
 * This updates ZCENTER parameter.
 *
 * Note: This assumes that the fitting library is using the 
 *       following convention:
 *
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = ZCENTER;
 *  delta[4] = BACKGROUND;
 */
void mpArbUpdateZ(mpFit *mp_fit)
{
  int i,nc,zi;
  double delta,dz,p_ave,p_total,w;
  double *heights;
  peakData *peak;

  heights = mp_fit->heights;
  nc = mp_fit->n_channels;

  /* Calculate index into z-dependent weight values and do some range checking. */
  mpWeightIndex(mp_fit, &dz, &zi);
  
  /* Z parameter is a simple weighted average. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    w = mpWeightInterpolate(mp_fit->w_z, dz, zi, nc, i);
    p_ave += mp_fit->w_jacobian[i][3] * w * heights[i];
    p_total += w * heights[i];
  }
  delta = p_ave/p_total;

  for(i=0;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    mFitUpdateParam(peak, delta, ZCENTER);
    
    /* Force z value to stay in range. */
    mp_fit->fn_zrange(mp_fit->fit_data[i]);
  }
}


/*
 * mpArbUpdateFixed()
 *
 * Calculate weighted delta and update each channel for fitting
 * with the peak height fixed to be the same in all the channels.
 *
 * Note: This allows negative heights, which will get removed by fn_check().
 *
 * Note: This assumes that the fitting library is using the 
 *       following convention:
 *
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = ZCENTER;
 *  delta[4] = BACKGROUND;
 */
void mpArbUpdateFixed(mpFit *mp_fit)
{
  int i,nc,zi;
  double delta, dz, p_ave, p_total, w;
  fitData *fit_data_ch0;

  fit_data_ch0 = mp_fit->fit_data[0];
  nc = mp_fit->n_channels;

  mpWeightIndex(mp_fit, &dz, &zi);

  /* Height, this is a simple weighted average. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(VERBOSE){
      printf("mpAUF h %d %.3e\n", i, mp_fit->w_jacobian[i][0]);
    }
    w = mpWeightInterpolate(mp_fit->w_h, dz, zi, nc, i);
    p_ave += mp_fit->w_jacobian[i][0] * w;
    p_total += w;
  }
  delta = p_ave/p_total;
  
  mFitUpdateParam(fit_data_ch0->working_peak, delta, HEIGHT);
  for(i=1;i<nc;i++){
    mp_fit->fit_data[i]->working_peak->params[HEIGHT] = fit_data_ch0->working_peak->params[HEIGHT];
  }

  mpUpdate(mp_fit);
  mpArbUpdateZ(mp_fit);
}


/*
 * mpArbUpdateIndependent()
 *
 * Calculate weighted delta and update each channel for fitting
 * with independently adjustable peak heights.
 *
 * Note: This does not allow negative peak heights, so peaks will
 *       never be culled for being to short.
 *
 * Note: This assumes that the PSF fitting library is using the 
 *       following convention:
 *
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = ZCENTER;
 *  delta[4] = BACKGROUND;
 */
void mpArbUpdateIndependent(mpFit *mp_fit)
{
  int i,nc;
  peakData *peak;

  nc = mp_fit->n_channels;
  for(i=0;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    mFitUpdateParam(peak, mp_fit->w_jacobian[i][0], HEIGHT);

    /* Prevent small/negative peak heights. */
    if(peak->params[HEIGHT] < 0.01){
      peak->params[HEIGHT] = 0.01;
    }
    
    mp_fit->heights[i] = peak->params[HEIGHT];
  }

  mpUpdate(mp_fit);
  mpArbUpdateZ(mp_fit);
}
