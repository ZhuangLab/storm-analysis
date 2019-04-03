/*
 * Multi-plane fitting with Gaussian PSFs. For now the assumption is
 * this would be used for standard biplane / multiplane imaging so
 * the PSFs will be (more or less) symmetric.
 *
 * Most of the work is done sa_library/dao_fit.c
 *
 * Hazen 01/18
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../multi_plane/mp_fit.h"

#include "../sa_library/dao_fit.h"


void mpDaoInitialize2DChannel(mpFit *, double *, double *, double, double, int, int);
void mpDaoNewPeaks(mpFit *, double *, char *, int);
void mpDaoUpdate(mpFit *);


/*
 * mpDaoInitialize2DChannel()
 *
 * Initialize a single channel / plane for 2D Gaussian fitting.
 */
void mpDaoInitialize2DChannel(mpFit *mp_fit, double *rqe, double *variance, double width_min, double width_max, int roi_size, int channel)
{
  int jac_size;

  /* Specify how to add new peaks and how to cleanup. */
  if(channel == 0){
    mp_fit->fn_cleanup = &daoCleanup;
    mp_fit->fn_newpeaks = &daoNewPeaks;
    mp_fit->fn_peak_xi_yi = &mFitUpdate;
    mp_fit->fn_update = &mpDaoUpdate;
    mp_fit->fn_zrange = NULL;             /* This should never get called. */
  }
  
  /*
   * Initialize 2D Gaussian fitting for this channel / plane.
   */
  mp_fit->fit_data[channel] = daoInitialize(rqe,
					    variance,
					    mp_fit->tolerance,
					    mp_fit->im_size_x,
					    mp_fit->im_size_y,
					    roi_size);
  daoInitialize2D(mp_fit->fit_data[channel],
		  width_min,
		  width_max);

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
 * mpDaoNewPeaks()
 *
 * New peaks to fit.
 *
 * n_peaks is the number of peaks per channel.
 */
void mpDaoNewPeaks(mpFit *mp_fit, double *peak_params, char *p_type, int n_peaks)
{
  int i,j,k;
  int start,stop;
  double tx,ty;
  double *mapped_peak_params;
  
  if(VERBOSE){
    printf("mpNP %d\n", n_peaks);
  }

  start = mp_fit->fit_data[0]->nfit;
  stop = start + n_peaks;
  
  if(!strcmp(p_type, "finder") || !strcmp(p_type, "testing")){

    /* We'll use this to pass the mapped peak positions. */
    mapped_peak_params = (double *)malloc(sizeof(double)*n_peaks*4);

    /* Map peak positions & pass to each of the fitters. */
    for(i=0;i<mp_fit->n_channels;i++){
      if(i>0){
	for(j=0;j<n_peaks;j++){
	  k = 4*j;
	  tx = peak_params[k];
	  ty = peak_params[k+1];
	  mapped_peak_params[k] = mp_fit->xt_0toN[i*3] + tx*mp_fit->xt_0toN[i*3+1] + ty*mp_fit->xt_0toN[i*3+2];
	  mapped_peak_params[k+1] = mp_fit->yt_0toN[i*3] + tx*mp_fit->yt_0toN[i*3+1] + ty*mp_fit->yt_0toN[i*3+2];
	  mapped_peak_params[k+2] = peak_params[k+2];
	  mapped_peak_params[k+3] = peak_params[k+3];
	}
	mp_fit->fn_newpeaks(mp_fit->fit_data[i], mapped_peak_params, p_type, n_peaks);
      }
      else{
	mp_fit->fn_newpeaks(mp_fit->fit_data[i], peak_params, p_type, n_peaks);
      }
    }

    free(mapped_peak_params);
  }
  else{
    mapped_peak_params = (double *)malloc(sizeof(double)*n_peaks*7);

    /* Map peak positions & pass to each of the fitters. */
    for(i=0;i<mp_fit->n_channels;i++){
      if(i>0){
	for(j=0;j<n_peaks;j++){
	  k = 7*j;
	  tx = peak_params[k];
	  ty = peak_params[k+1];
	  mapped_peak_params[k] = mp_fit->xt_0toN[i*3] + tx*mp_fit->xt_0toN[i*3+1] + ty*mp_fit->xt_0toN[i*3+2];
	  mapped_peak_params[k+1] = mp_fit->yt_0toN[i*3] + tx*mp_fit->yt_0toN[i*3+1] + ty*mp_fit->yt_0toN[i*3+2];
	  mapped_peak_params[k+2] = peak_params[k+2];
	  mapped_peak_params[k+3] = peak_params[k+3];
	  mapped_peak_params[k+4] = peak_params[k+4];
	  mapped_peak_params[k+5] = peak_params[k+5];
	  mapped_peak_params[k+6] = peak_params[k+6];
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
 * mpDaoUpdate()
 *
 * Calculate weighted delta for XCENTER, YCENTER and update all the
 * parameters in each channel.
 *
 * Note: This does not allow negative peak heights, so peaks will
 *       never be culled for being to short.
 *
 * Note: This assumes that the fitting library is using the 
 *       following convention:
 *
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = XWIDTH;
 *  delta[4] = BACKGROUND;
 */
void mpDaoUpdate(mpFit *mp_fit)
{
  int i,nc;
  peakData *peak;

  /* Update heights, these float independently. */
  nc = mp_fit->n_channels;
  for(i=0;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    mFitUpdateParam(peak, mp_fit->w_jacobian[i][0], HEIGHT);

    /* Prevent small/negative peak heights. */
    if(peak->params[HEIGHT] < 0.01){
      peak->params[HEIGHT] = 0.01;
    }

    /* This is used for weighting the per channel XCENTER, YCENTER values. */
    mp_fit->heights[i] = peak->params[HEIGHT];
  }
  
  /* Widths float independently. */
  for(i=0;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    mFitUpdateParam(peak, mp_fit->w_jacobian[i][3], XWIDTH);
    
    /* 
     * Range clamp. 
     *
     * Basically the problem is that in some channels the peak intensity will 
     * be really low so fitting the width is basically impossible anyway.
     */
    daoCheck2D(mp_fit->fit_data[i]);
  }

  mpUpdate(mp_fit);
}
