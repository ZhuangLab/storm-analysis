/*
 * Fit multiple, possibly overlapping, PSFs to image data 
 * from multiple planes.
 *
 * The expectation is that there will be n_channels copies of
 * each input peak, organized by channel, so for example
 * if there are 3 peaks and 2 channels the peak array would 
 * be [peak1_c1, peak2_c1, peak3_c1, peak1_c2, peak2_c2, peak2_c3].
 * This analysis will then keep the groups of peaks in sync,
 * i.e. peak1_c1 and peak1_c2 will have the same peak status
 * (RUNNING, CONVERGED, ERROR), z value and possibly height. And 
 * their x, y coordinates will be the same after affine 
 * transformation.
 *
 * Proper initialization involves multiple steps:
 *  1. mpInitialize()
 *  2. mpInitializeXXChannel() for each channel.
 *  3. mpSetTransforms() to configure affine transforms between 
 *       channels.
 *  4. mpSetWeights() to set z dependent channel parameter
 *       weighting factors.
 *  5. mpSetWeightsIndexing() to set how to go from a peaks
 *       z value to the correct index in the weighting array.
 *
 * Hazen 10/17
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mp_fit.h"


/*
 * mpCheckError()
 *
 * Check if any peaks in a group are in the ERROR state. If one is
 * remove all the other peaks in the group and mark them as also
 * being in the ERROR state.
 */
void mpCheckError(mpFit *mp_fit, int peak_index)
{
  int is_bad,j;
  fitData *fit_data;

  is_bad = 0;
  for(j=0;j<mp_fit->n_channels;j++){
    if(mp_fit->fit_data[j]->fit[peak_index].status == ERROR){
      is_bad = 1;
      break;
    }
  }
    
  if(is_bad){
    for(j=0;j<mp_fit->n_channels;j++){

      /* Check if we need to subtract this peak out of the image. */
      if(mp_fit->fit_data[j]->fit[peak_index].status != ERROR){
	fit_data = mp_fit->fit_data[j];
	fit_data->fn_copy_peak(fit_data, &fit_data->fit[peak_index], fit_data->working_peak);
	mFitSubtractPeak(fit_data);
	fit_data->fn_copy_peak(fit_data, fit_data->working_peak, &fit_data->fit[peak_index]);
      }
      mp_fit->fit_data[j]->fit[peak_index].status = ERROR;
    }
  } 
}


/*
 * mpCleanup()
 *
 * Clean up at the end of the analysis.
 */
void mpCleanup(mpFit *mp_fit)
{
  int i;

  /* Call PSF specific cleanup. */
  for(i=0;i<mp_fit->n_channels;i++){
    mp_fit->fn_cleanup(mp_fit->fit_data[i]);
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
  free(mp_fit->heights);

  /* Free jacobian / hessian storage. */
  for(i=0;i<mp_fit->n_channels;i++){
    free(mp_fit->jacobian[i]);
    free(mp_fit->w_jacobian[i]);
    free(mp_fit->hessian[i]);
    free(mp_fit->w_hessian[i]);
  }
  free(mp_fit->jacobian);
  free(mp_fit->w_jacobian);
  free(mp_fit->hessian);
  free(mp_fit->w_hessian);
  
  free(mp_fit->fit_data);
  free(mp_fit);
}


/*
 * mpCopyFromWorking()
 *
 * Copy the working peak into the indicated peak. This will also
 * set the status of all the (paired) peaks to the same value.
 */
void mpCopyFromWorking(mpFit *mp_fit, int index, int status)
{
  int i;
  fitData *fit_data;

  for(i=0;i<mp_fit->n_channels;i++){
    fit_data = mp_fit->fit_data[i];
    fit_data->working_peak->status = status;
    fit_data->fn_copy_peak(fit_data, fit_data->working_peak, &fit_data->fit[index]);
  }
}


/*
 * mpInitialize()
 *
 * Create and return the mpFit structure to use for fitting.
 */
mpFit *mpInitialize(double tolerance, int n_channels, int independent_heights, int im_size_x, int im_size_y)
{
  mpFit *mp_fit;

  mp_fit = (mpFit *)malloc(sizeof(mpFit));

  mp_fit->im_size_x = im_size_x;
  mp_fit->im_size_y = im_size_y;
  mp_fit->independent_heights = independent_heights;
  mp_fit->n_channels = n_channels;
  mp_fit->tolerance = tolerance;
  mp_fit->w_z_offset = 0.0;
  mp_fit->w_z_scale = 0.0;
  mp_fit->zmin = 0.0;
  mp_fit->zmax = 0.0;

  mp_fit->xt_0toN = (double *)malloc(3*n_channels*sizeof(double));
  mp_fit->yt_0toN = (double *)malloc(3*n_channels*sizeof(double));
  mp_fit->xt_Nto0 = (double *)malloc(3*n_channels*sizeof(double));
  mp_fit->yt_Nto0 = (double *)malloc(3*n_channels*sizeof(double));

  mp_fit->fit_data = (fitData **)malloc(n_channels*sizeof(fitData*));

  mp_fit->jacobian = (double **)malloc(n_channels*sizeof(double *));
  mp_fit->w_jacobian = (double **)malloc(n_channels*sizeof(double *));
  mp_fit->hessian = (double **)malloc(n_channels*sizeof(double *));
  mp_fit->w_hessian = (double **)malloc(n_channels*sizeof(double *));
    
  return mp_fit;
}


/*
 * mpIterateLM()
 *
 * Perform a single cycle of fitting for each localization using the
 * Levenberg-Marquardt algorithm.
 */
void mpIterateLM(mpFit *mp_fit)
{
  int i,j,k,l,m,n,nfit;
  int info,is_bad;
  int n_add;
  double current_error,starting_error;
  fitData *fit_data;

  nfit = mp_fit->fit_data[0]->nfit;
  
  if(VERBOSE){
    printf("mpILM, nfit = %d\n", nfit);
  }

  for(i=0;i<nfit;i++){
      
    /* Skip ahead if this peak is not RUNNING. */
    if(mp_fit->fit_data[0]->fit[i].status != RUNNING){
      continue;
    }

    if(VERBOSE){
      printf("mpILM index = %d\n", i);
    }
    
    /* 
     * This is for debugging, to make sure that we are adding and 
     * subtracting the right number of times.
     */    
    n_add = mp_fit->n_channels;

    /*
     * Copy peak, calculate jacobian and hessian and subtract.
     */
    starting_error = 0.0;
    for(j=0;j<mp_fit->n_channels;j++){
      fit_data = mp_fit->fit_data[j];
      
      /* Copy current peak into working peak. */
      fit_data->fn_copy_peak(fit_data, &fit_data->fit[i], fit_data->working_peak);

      /* Calculate current error. */
      mFitCalcErr(fit_data);
      starting_error += fit_data->working_peak->error;
      
      /* Calculate Jacobian and Hessian. This is expected to use 'working_peak'. */
      fit_data->fn_calc_JH(fit_data, mp_fit->jacobian[j], mp_fit->hessian[j]);
    
      /* Subtract current peak out of image. This is expected to use 'working_peak'. */
      mFitSubtractPeak(fit_data);
      n_add--;
    }
    
    /*
     * Try and improve paired peak parameters.
     */
    j = 0;
    while(1){
      j++;

      if(VERBOSE){
	printf("mpILM cycle %d %d %d\n", i, j, n_add);
      }
      
      is_bad = 0;

      /* 0. Check if we are stuck on this peak, error it out if we are. */
      if(mp_fit->fit_data[0]->working_peak->lambda > LAMBDAMAX){
	for(k=0;k<mp_fit->n_channels;k++){
	  fit_data = mp_fit->fit_data[k];
	  fit_data->n_lost++;
	  fit_data->working_peak->status = ERROR;
	}
	break;
      }

      /* 1. Reset status, as it may have changed on a previous pass through this loop. */
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
	fit_data->working_peak->status = RUNNING;
      }
      
      /* 2. Solve for the update vectors. */
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];

	/* Update peak iterations counter. */
	fit_data->working_peak->iterations++;
	
	/* Update total fitting iterations counter. */
	fit_data->n_iterations++;

	/* Copy Jacobian and Hessian. */
	for(l=0;l<fit_data->jac_size;l++){
	  mp_fit->w_jacobian[k][l] = mp_fit->jacobian[k][l];
	  m = l*fit_data->jac_size;
	  for(n=0;n<fit_data->jac_size;n++){
	    if (l == n){
	      mp_fit->w_hessian[k][m+n] = (1.0 + fit_data->working_peak->lambda) * mp_fit->hessian[k][m+n];
	    }
	    else{
	      mp_fit->w_hessian[k][m+n] = mp_fit->hessian[k][m+n];
	    }
	  }
	}
      
	/*  Solve for update. Note that this also changes jacobian. */
	info = mFitSolve(mp_fit->w_hessian[k], mp_fit->w_jacobian[k], fit_data->jac_size);
	
	/* If the solver failed, set is_bad = 1 and exit this loop. */
	if(info!=0){
	  is_bad = 1;
	  fit_data->n_dposv++;
	  if(VERBOSE){
	    printf("mpILM mFitSolve() failed %d %d\n", i, info);
	  }
	  break;
	}
      }

      /* If the solver failed then start over again with a higher lambda for all paired peaks. */
      if(is_bad){
	for(k=0;k<mp_fit->n_channels;k++){
	  fit_data = mp_fit->fit_data[k];
	  fit_data->working_peak->status = ERROR;
	  fit_data->working_peak->lambda = fit_data->working_peak->lambda * LAMBDAUP;
	}
	continue;
      }

      /* 3. Update working peaks. This will use the deltas in w_jacobian. */
      mp_fit->fn_update(mp_fit);
      
      /* 4. Check that the peaks are still in the image, etc.. */
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
	if(fit_data->fn_check(fit_data)){
	  is_bad = 1;
	  if(VERBOSE){
	    printf("mpILM fn_check() failed %d\n", i);
	  }
	}
      }

      /* If the peak parameter check failed start over again with a higher lambda for all paired peaks. */
      if(is_bad){
	mpResetWorkingPeaks(mp_fit, i);
	continue;
      }

      /* 5. Add working peaks back to the fit image. */
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
	fit_data->fn_calc_peak_shape(fit_data);
	mFitAddPeak(fit_data);
	n_add++;
      }

      /* 6. Calculate updated error. */
      current_error = 0.0;
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];

	/* 
	 * mFitCalcErr() updates fit_data->working_peak->error, and returns 0 for
	 * success, 1 for failure.
	 */
	if(mFitCalcErr(fit_data)){
	  is_bad = 1;
	  if(VERBOSE){
	    printf("mpILM mFitCalcErr() failed\n");
	  }
	}
	else{
	  current_error += fit_data->working_peak->error;
	}
      }

      /* If the peak error calculation failed start over again with a higher lambda for all paired peaks. */
      if(is_bad){
	
	/* Undo peak addition. */
	for(k=0;k<mp_fit->n_channels;k++){
	  fit_data = mp_fit->fit_data[k];
	  mFitSubtractPeak(fit_data);
	  n_add--;
	}
	
	/* Reset working peaks. */
	mpResetWorkingPeaks(mp_fit, i);

	/* Try again. */
	continue;
      }

      /* 7. Check that the error is decreasing. */

      /* If the peak error has increased then start over again with a higher lambda for all paired peaks. */
      if(current_error > starting_error){
	
	if(VERBOSE){
	  printf("mpILM increasing error %.6e %.6e\n", current_error, starting_error);
	}
	
	/* 
	 * Check for error convergence. 
	 *
	 * Usually this will happen because the lambda term has gotten so 
	 * large that the peak will barely move in the update.
	 */
      	if (((current_error - starting_error)/starting_error) < mp_fit->tolerance){
	  for(k=0;k<mp_fit->n_channels;k++){
	    mp_fit->fit_data[k]->working_peak->status = CONVERGED;
	  }
	  break;
	}
	else{
	  
	  /* Undo peak addition, and increment counter. */
	  for(k=0;k<mp_fit->n_channels;k++){
	    fit_data = mp_fit->fit_data[k];
	    fit_data->n_non_decr++;
	    mFitSubtractPeak(fit_data);
	    n_add--;
	  }
	
	  /* Reset working peaks. */
	  mpResetWorkingPeaks(mp_fit, i);

	  /* Try again. */
	  continue;
	}
      }
      else{
	
	if(VERBOSE){
	  printf("mpILM decreasing error %.6e %.6e\n", current_error, starting_error);
	}

	/* Check for error convergence. */
      	if (((starting_error - current_error)/starting_error) < mp_fit->tolerance){
	  for(k=0;k<mp_fit->n_channels;k++){
	    mp_fit->fit_data[k]->working_peak->status = CONVERGED;
	  }
	}

	/* Otherwise reduce lambda. */
	else if(mp_fit->fit_data[0]->working_peak->lambda > LAMBDAMIN){
	  for(k=0;k<mp_fit->n_channels;k++){
	    fit_data = mp_fit->fit_data[k];
	    fit_data->working_peak->lambda = fit_data->working_peak->lambda * LAMBDADOWN;
	  }
	}
	break;
      }
    }
	  
    /* We expect n_add to be n_channels if there were no errors, 0 otherwise. */
    if(TESTING){
      if(mp_fit->fit_data[0]->working_peak->status == ERROR){
	if(n_add != 0){
	  printf("Problem detected in peak addition / subtraction logic, status == ERROR, counts = %d\n", n_add);
	  exit(EXIT_FAILURE);
	}
      }
      else{
	if(n_add != mp_fit->n_channels){
	  printf("Problem detected in peak addition / subtraction logic, status != ERROR, counts = %d\n", n_add);
	  exit(EXIT_FAILURE);
	}
      }
    }
    
    /* Copy updated working peak back into current peak. */
    mpCopyFromWorking(mp_fit, i, mp_fit->fit_data[0]->working_peak->status);
  }

  /* 
   * Recenter peaks. This may throw peaks into the ERROR state, so after
   * we do this we need to synchronize the peak ERROR state across all the 
   * channels.
   */
  for(j=0;j<mp_fit->n_channels;j++){
    mFitRecenterPeaks(mp_fit->fit_data[j]);
  }

  for(i=0;i<nfit;i++){
    mpCheckError(mp_fit, i);
    
    /* Sanity check. */
    if(TESTING){
      for(j=1;j<mp_fit->n_channels;j++){
	if(mp_fit->fit_data[j]->fit[i].status != mp_fit->fit_data[0]->fit[i].status){
	  printf("Peak channel statuses do not agree for peak %d!\n", i);
	  exit(EXIT_FAILURE);
	}
      }
    }
  }
}


/*
 * mpResetWorkingPeaks()
 *
 * Restore working peaks to their original state, but with larger lambda
 * and status ERROR. This is used by mpIterateLM().
 */
void mpResetWorkingPeaks(mpFit *mp_fit, int index)
{
  int i,tmp_added;
  double tmp_lambda;
  fitData *fit_data;
  
  for(i=0;i<mp_fit->n_channels;i++){
    fit_data = mp_fit->fit_data[i];
    
    /* Reset peak (it was changed by fn_update()), increase lambda. */
    tmp_added = fit_data->working_peak->added;
    tmp_lambda = fit_data->working_peak->lambda;
    fit_data->fn_copy_peak(fit_data, &fit_data->fit[index], fit_data->working_peak);
    fit_data->working_peak->added = tmp_added;
    fit_data->working_peak->lambda = tmp_lambda * LAMBDAUP;

    /* Set status to ERROR in case this is the last iteration. */
    fit_data->working_peak->status = ERROR;
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
 */
void mpSetWeights(mpFit *mp_fit, double *w_bg, double *w_h, double *w_x, double *w_y, double *w_z, int z_size)
{
  int i,n;
  
  mp_fit->n_weights = z_size;
  
  /* Allocate storage. */
  n = mp_fit->n_channels*z_size;
  mp_fit->w_bg = (double *)malloc(sizeof(double)*n);
  mp_fit->w_h = (double *)malloc(sizeof(double)*n);
  mp_fit->w_x = (double *)malloc(sizeof(double)*n);
  mp_fit->w_y = (double *)malloc(sizeof(double)*n);
  mp_fit->w_z = (double *)malloc(sizeof(double)*n);
  mp_fit->heights = (double *)malloc(sizeof(double)*mp_fit->n_channels);
  
  /* Copy values. */
  for(i=0;i<n;i++){
    mp_fit->w_bg[i] = w_bg[i];
    mp_fit->w_h[i] = w_h[i];
    mp_fit->w_x[i] = w_x[i];
    mp_fit->w_y[i] = w_y[i];
    mp_fit->w_z[i] = w_z[i];
  }

  /* Set initial height weighting values to 1.0 for fixed (relative) height fitting. */
  for(i=0;i<mp_fit->n_channels;i++){
    mp_fit->heights[i] = 1.0;
  }
}

/*
 * mpSetWeightsIndexing()
 *
 * Set the values to use for conversion of a peak Z position to 
 * an index into the weights arrays.
 */
void mpSetWeightsIndexing(mpFit *mp_fit, double z_offset, double z_scale)
{
  mp_fit->w_z_offset = z_offset;
  mp_fit->w_z_scale = z_scale;
}

/*
 * mpUpdate()
 *
 * This updates the XCENTER, YCENTER and BACKGROUND parameters.
 *
 * Calculate weighted delta and update each channel. The weights 
 * cover the z range of the PSF and include the two end-points.
 * We use linear interpolation between the points in the weights
 * array.
 *
 * mp_fit->heights should be all 1.0 for fixed (relative) heights.
 *
 * Note: This assumes that the fitting library is using the 
 *       following convention:
 *
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = ZCENTER or XWIDTH;
 *  delta[4] = BACKGROUND;
 */
void mpUpdate(mpFit *mp_fit)
{
  int i,nc,zi;
  double delta,dz,p_ave,p_total,t,w,xoff,yoff;
  double *params_ch0,*heights;
  peakData *peak;
  fitData *fit_data_ch0;

  heights = mp_fit->heights;
  fit_data_ch0 = mp_fit->fit_data[0];
  params_ch0 = fit_data_ch0->working_peak->params;
  xoff = fit_data_ch0->xoff;
  yoff = fit_data_ch0->yoff;
  
  nc = mp_fit->n_channels;

  /* Calculate index into z-dependent weight values and do some range checking. */
  mpWeightIndex(mp_fit, &dz, &zi);
  
  /*
   * X parameters depends on the mapping.
   */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(VERBOSE){
      printf("mpU x %d %.3e %.3e", i, heights[i], mp_fit->xt_Nto0[i*3+1]);
      printf("mpU %.3e %.3e %.3e\n", mp_fit->w_jacobian[i][2], mp_fit->xt_Nto0[i*3+2], mp_fit->w_jacobian[i][1]);
    }
    w = mpWeightInterpolate(mp_fit->w_x, dz, zi, nc, i);
    delta = mp_fit->xt_Nto0[i*3+1] * mp_fit->w_jacobian[i][1];
    delta += mp_fit->xt_Nto0[i*3+2] * mp_fit->w_jacobian[i][2];
    p_ave += delta * w * heights[i];
    p_total += w * heights[i];
  }
  delta = p_ave/p_total;
  mFitUpdateParam(fit_data_ch0->working_peak, delta, XCENTER);

  /* Y parameters also depend on the mapping. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(VERBOSE){
      printf("mpU y %d %.3e %.3e", i, heights[i], mp_fit->yt_Nto0[i*3+1]);
      printf("mpU %.3e %.3e %.3e\n", mp_fit->w_jacobian[i][2], mp_fit->yt_Nto0[i*3+2], mp_fit->w_jacobian[i][1]);
    }
    w = mpWeightInterpolate(mp_fit->w_y, dz, zi, nc, i);
    delta = mp_fit->yt_Nto0[i*3+1] * mp_fit->w_jacobian[i][1];
    delta += mp_fit->yt_Nto0[i*3+2] * mp_fit->w_jacobian[i][2];
    p_ave += delta * w * heights[i];
    p_total += w * heights[i];
  }
  delta = p_ave/p_total;
  mFitUpdateParam(fit_data_ch0->working_peak, delta, YCENTER);  

  /* 
   * Use mapping to update peak locations in the remaining channels.
   *
   * Note: Spliner uses the upper left corner as 0,0 so we need to adjust
   *       to the center, transform, then adjust back. This is particularly
   *       important if one channel is inverted relative to another.
   */
  for(i=1;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    
    t = mp_fit->xt_0toN[i*3];
    t += mp_fit->xt_0toN[i*3+1] * (params_ch0[XCENTER]+xoff);
    t += mp_fit->xt_0toN[i*3+2] * (params_ch0[YCENTER]+yoff);
    peak->params[XCENTER] = t-xoff;

    t = mp_fit->yt_0toN[i*3];
    t += mp_fit->yt_0toN[i*3+1] * (params_ch0[XCENTER]+xoff);
    t += mp_fit->yt_0toN[i*3+2] * (params_ch0[YCENTER]+yoff);
    peak->params[YCENTER] = t-yoff;
  }

  /* Update peak (integer) location with hysteresis. */
  for(i=0;i<nc;i++){
    mp_fit->fn_peak_xi_yi(mp_fit->fit_data[i]->working_peak);
  }

  /* Backgrounds float independently. */
  for(i=0;i<nc;i++){
    mFitUpdateParam(mp_fit->fit_data[i]->working_peak, mp_fit->w_jacobian[i][4], BACKGROUND);
  }
}


/*
 * mpWeightIndex()
 *
 * Determine which weight value to use.
 */
void mpWeightIndex(mpFit *mp_fit, double *dz, int *zi)
{
  fitData *fit_data;

  fit_data = mp_fit->fit_data[0];
    
  *zi = (int)(mp_fit->w_z_scale * (fit_data->working_peak->params[ZCENTER] - mp_fit->w_z_offset));
  *dz = (mp_fit->w_z_scale * (fit_data->working_peak->params[ZCENTER] - mp_fit->w_z_offset)) - (double)*zi;
  
  if(*zi<0){
    if(TESTING){
      printf("Negative weight index detected %d (%.2f)\n!", *zi, fit_data->working_peak->params[ZCENTER]);
      exit(EXIT_FAILURE);
    }
    *zi = 0;
  }
  if(*zi>(mp_fit->n_weights-2)){
    if(TESTING){
      printf("Out of range weight index detected %d (%.2f)\n!", *zi, fit_data->working_peak->params[ZCENTER]);
      exit(EXIT_FAILURE);
    }
    *zi = mp_fit->n_weights-2;
  }
}

/*
 * mpWeightInterpolate()
 *
 * Linearly interpolates a weights array.
 */
double mpWeightInterpolate(double *weights, double dz, int z_index, int n_channels, int channel)
{
  double w1,w2,dw;

  /* Get weight as the two end-points. */
  w1 = weights[n_channels*z_index+channel];
  w2 = weights[n_channels*(z_index+1)+channel];

  dw = w2 - w1;

  return w1 + dw*dz;
  /* return w1; */
}
