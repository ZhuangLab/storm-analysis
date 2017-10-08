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
 * (RUNNING, CONVERGED, ERROR), z value and possibly height. And 
 * their x, y coordinates will be the same after affine 
 * transformation.
 *
 * Hazen 06/17
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../spliner/cubic_fit.h"

typedef struct mpFit
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
  double *heights;              /* Per channel heights for parameter weighting. */

  double **jacobian;            /* Storage for the jacobian calculations. */
  double **w_jacobian;          /* Storage for copies of the jacobians. */
  double **hessian;             /* Storage for the hessian calculations. */
  double **w_hessian;           /* Storage for copies of the jacobians. */
  
  fitData **fit_data;           /* Array of pointers to fitData structures. */

  void (*fn_update)(struct mpFit *); /* Function for updating the parameters of the working peaks. */
  
} mpFit;


void mpCleanup(mpFit *);
void mpCopyFromWorking(mpFit *, int, int);
//void mpCopyToWorking(mpFit *, int);
void mpGetFitImage(mpFit *, double *, int);
void mpGetResults(mpFit *, double *);
int mpGetUnconverged(mpFit *);
mpFit *mpInitialize(double *, double, int, int, int, int);
void mpInitializeChannel(mpFit *, splineData *, double *, int);
void mpIterateLM(mpFit *);
void mpIterateOriginal(mpFit *);
void mpNewImage(mpFit *, double *, int);
void mpNewPeaks(mpFit *, double *, int);
void mpResetWorkingPeaks(mpFit *, int);
void mpSetTransforms(mpFit *, double *, double *, double *, double *);
void mpSetWeights(mpFit *, double *, double *, double *, double *, double *);
void mpUpdate(mpFit *);
void mpUpdateFixed(mpFit *);
void mpUpdateIndependent(mpFit *);


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
    fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[index]);
  }
}


/*
 * mpCopyToWorking()
 *
 * Copy the indicated peak to the working peak.
 */
/*
void mpCopyToWorking(mpFit *mp_fit, int index)
{
  int i;
  fitData *fit_data;

  for(i=0;i<mp_fit->n_channels;i++){
    fit_data = mp_fit->fit_data[i];
    fit_data->fn_copy_peak(&fit_data->fit[index], fit_data->working_peak);
  }
}
*/


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
  int i;
  /*
   * We only need to check the peaks for channel 0 as the
   * peaks for other channels will have the same state.
   */
  if(VERBOSE){
    i = mFitGetUnconverged(mp_fit->fit_data[0]);
    printf("mpGU %d\n", i);
  }
  return mFitGetUnconverged(mp_fit->fit_data[0]);
}


/*
 * mpInitialize()
 *
 * Create and return the mpFit structure to use for fitting.
 */
mpFit *mpInitialize(double *clamp, double tolerance, int n_channels, int independent_heights, int im_size_x, int im_size_y)
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

  mp_fit->jacobian = (double **)malloc(n_channels*sizeof(double *));
  mp_fit->w_jacobian = (double **)malloc(n_channels*sizeof(double *));
  mp_fit->hessian = (double **)malloc(n_channels*sizeof(double *));
  mp_fit->w_hessian = (double **)malloc(n_channels*sizeof(double *));

  if(independent_heights){
    mp_fit->fn_update = &mpUpdateIndependent;
  }
  else{
    mp_fit->fn_update = &mpUpdateFixed;
  }
    
  return mp_fit;
}


/*
 * mpInitializeChannel()
 *
 * Initialize a single channel / plane for 3D spline fitting.
 */
void mpInitializeChannel(mpFit *mp_fit, splineData *spline_data, double *variance, int channel)
{
  int jac_size;
  
  /*
   * Initialize spliner fitting for this channel / plane.
   */
  mp_fit->fit_data[channel] = cfInitialize(spline_data,
					   variance,
					   mp_fit->clamp_start,
					   mp_fit->tolerance,
					   mp_fit->im_size_x,
					   mp_fit->im_size_y);
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
 * mpIterateLM()
 *
 * Perform a single cycle of fitting for each localization using the
 * Levenberg-Marquardt algorithm.
 */
void mpIterateLM(mpFit *mp_fit)
{
  int i,j,k,l,m,n;
  int info,is_bad,is_converged;
  int n_add;
  double error, error_old;
  fitData *fit_data;

  if(VERBOSE){
    printf("mpILM, nfit = %d\n", mp_fit->nfit);
  }

  for(i=0;i<mp_fit->nfit;i++){
      
    /* Skip ahead if this peak is not RUNNING. */
    if(mp_fit->fit_data[0]->fit[i].status != RUNNING){
      continue;
    }

    if(VERBOSE){
      printf("mpILM index = %d\n", i);
    }
    
    /* 
     * This is for debugging, to make sure that we adding and subtracting
     * the right number of times.
     */    
    n_add = mp_fit->n_channels;

    /*
     * Copy peak, calculate jacobian and hessian and subtract.
     */
    for(j=0;j<mp_fit->n_channels;j++){
      fit_data = mp_fit->fit_data[j];
      
      /* Copy current peak into working peak. */
      fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);

      /* Calculate Jacobian and Hessian. This is expected to use 'working_peak'. */
      fit_data->fn_calc_JH(fit_data, mp_fit->jacobian[j], mp_fit->hessian[j]);
    
      /* Subtract current peak out of image. This is expected to use 'working_peak'. */
      fit_data->fn_subtract_peak(fit_data);
      n_add--;
    }
    
    /*
     * Try and improve paired peak parameters.
     */
    for(j=0;j<=MAXCYCLES;j++){
      is_bad = 0;

      /* 1. Reset status, as it may have changed on a previous pass through this loop. */
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
	fit_data->working_peak->status = RUNNING;
      }
      
      /* 2. Solve for the update vectors. */
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
    
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
	    printf(" mFitSolve() failed %d %d\n", i, info);
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
	    printf(" fn_check() failed %d\n", i);
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
	fit_data->fn_add_peak(fit_data);
	n_add++;
      }

      /* 6. Calculate updated error. */
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
	if(mFitCalcErr(fit_data)){
	  is_bad = 1;
	  if(VERBOSE){
	    printf(" mFitCalcErr() failed\n");
	  }
	}
      }

      /* If the peak error calculation failed start over again with a higher lambda for all paired peaks. */
      if(is_bad){
	
	/* Undo peak addition. */
	for(k=0;k<mp_fit->n_channels;k++){
	  fit_data = mp_fit->fit_data[k];
	  fit_data->fn_subtract_peak(fit_data);
	  n_add--;
	}
	
	/* Reset working peaks. */
	mpResetWorkingPeaks(mp_fit, i);

	/* Try again. */
	continue;
      }

      /* 7. Check if all the peaks have converged. */
      is_converged = 1;
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
	if(fit_data->working_peak->status != CONVERGED){
	  is_converged = 0;
	}
      }

      /* If all the peaks have converged then go to the next peak. */
      if(is_converged){
	if(VERBOSE){
	  printf("All converged\n");
	}
	break;
      }

      /* 8. Check that the error is decreasing. */
      error = 0.0;
      error_old = 0.0;
      for(k=0;k<mp_fit->n_channels;k++){
	fit_data = mp_fit->fit_data[k];
	error += fit_data->working_peak->error;
	error_old += fit_data->working_peak->error_old;
      }

      /* If the peak error has increased then start over again with a higher lambda for all paired peaks. */
      if(error > error_old){
	if(j<MAXCYCLES){
	  
	  if(VERBOSE){
	    printf("Error did not decrease %.2f %.2f\n", error, error_old);
	  }
	  
	  /* Undo peak addition, and increment counter. */
	  for(k=0;k<mp_fit->n_channels;k++){
	    fit_data = mp_fit->fit_data[k];
	    fit_data->n_non_decr++;
	    fit_data->fn_subtract_peak(fit_data);
	    n_add--;
	  }
	
	  /* Reset working peaks. */
	  mpResetWorkingPeaks(mp_fit, i);

	  /* Try again. */
	  continue;
	}
	else{
	  /* Note that even though the peak has not improved we are leaving it where it is and moving on. */
	  if(TESTING){
	    printf("Reached max cycles with no improvement in peak error for %d\n", i);
	  }
	}
      }
      else{
	/* Reduce lambda and exit the loop. */
	for(k=0;k<mp_fit->n_channels;k++){
	  fit_data = mp_fit->fit_data[k];
	  fit_data->working_peak->lambda = fit_data->working_peak->lambda * LAMBDADOWN;
	}
	break;
      }
    }
	  
    /* We expect n_add to be n_channels if there were no errors, 0 otherwise. */
    if(TESTING){
      if(mp_fit->fit_data[0]->working_peak->status == ERROR){
	if(n_add != 0){
	  printf("Problem detected in peak addition / subtraction logic, status == ERROR, counts = %d\n", n_add);
	  printf("Exiting now\n");
	  exit(1);
	}
      }
      else{
	if(n_add != mp_fit->n_channels){
	  printf("Problem detected in peak addition / subtraction logic, status != ERROR, counts = %d\n", n_add);
	  printf("Exiting now\n");
	  exit(1);
	}
      }
    }
    
    /* Copy updated working peak back into current peak. */
    mpCopyFromWorking(mp_fit, i, mp_fit->fit_data[0]->working_peak->status);
  }
}


/*
 * mpIterateOriginal()
 *
 * Perform a single cycle of fitting for each localization using
 * the original 3D-DAOSTORM like algorithm.
 */
void mpIterateOriginal(mpFit *mp_fit)
{
  int i,j;
  int info,is_bad,is_converged;
  fitData *fit_data;

  if(VERBOSE){
    printf("mpIO %d\n", mp_fit->nfit);
  }

  /*
   * 1. Calculate updated peaks.
   */
  for(i=0;i<mp_fit->nfit;i++){
      
    /* Skip ahead if this peak is not RUNNING. */
    if(mp_fit->fit_data[0]->fit[i].status != RUNNING){
      continue;
    }

    if(VERBOSE){
      printf("mpIO %d\n", i);
    }

    /*
     * Calculate update vector for each channel.
     */
    is_bad = 0;
    for(j=0;j<mp_fit->n_channels;j++){
      fit_data = mp_fit->fit_data[j];
      
      /* Copy current peak into working peak. */
      fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);

      /* Calculate Jacobian and Hessian. This is expected to use 'working_peak'. */
      fit_data->fn_calc_JH(fit_data, mp_fit->w_jacobian[j], mp_fit->w_hessian[j]);
    
      /* Subtract current peak out of image. This is expected to use 'working_peak'. */
      fit_data->fn_subtract_peak(fit_data);
    
      /* Update total fitting iterations counter. */
      fit_data->n_iterations++;

      /*  Solve for update. Note that this also changes jacobian. */
      info = mFitSolve(mp_fit->w_hessian[j], mp_fit->w_jacobian[j], fit_data->jac_size);

      /* If the solver failed, set is_bad = 1 and exit this loop. */
      if(info!=0){
	is_bad = 1;
	fit_data->n_dposv++;
	if(VERBOSE){
	  printf(" mFitSolve() failed %d %d\n", i, info);
	}
	break;
      }
    }

    /* 
     * If the solver failed for any peak, mark them all bad and go to
     * the next peak.
     */
    if(is_bad){
      mpCopyFromWorking(mp_fit, i, ERROR);
      continue;
    }

    /* 
     * Update parameters of working peaks. This will use the deltas 
     * in w_jacobian.
     */
    mp_fit->fn_update(mp_fit);

    /* 
     * Check that peaks are still in the image, etc.. The fn_check function
     * should return 0 if everything is okay.
     */    
    for(j=0;j<mp_fit->n_channels;j++){
      fit_data = mp_fit->fit_data[j];
      if(fit_data->fn_check(fit_data)){
	is_bad = 1;
	if(VERBOSE){
	  printf(" fn_check() failed %d\n", i);
	}
      }
    }

    /* 
     * If fn_check() failed for any peak, mark them all bad and go to
     * the next peak.
     */
    if(is_bad){
      mpCopyFromWorking(mp_fit, i, ERROR);
      continue;
    }

    /* Add working peaks back to image and copy back to current peak. */
    for(j=0;j<mp_fit->n_channels;j++){
      fit_data = mp_fit->fit_data[j];
      fit_data->fn_add_peak(fit_data);
      fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
    }
  }
    
  /*
   * 2. Calculate peak errors.
   */
  for(i=0;i<mp_fit->nfit;i++){
    
    /* Skip ahead if this peak is not RUNNING. */
    if(mp_fit->fit_data[0]->fit[i].status != RUNNING){
      continue;
    }

    /*  Calculate errors for the working peaks. */
    is_bad = 0;
    is_converged = 1;
    for(j=0;j<mp_fit->n_channels;j++){
      fit_data = mp_fit->fit_data[j];
      fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);
      if(mFitCalcErr(fit_data)){
	is_bad = 1;
	if(VERBOSE){
	  printf(" mFitCalcErr() failed %d\n", i);
	}
      }
      if(fit_data->working_peak->status != CONVERGED){
	is_converged = 0;
      }
      fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
    }

    /* If one peak has not converged then mark them all as not converged. */
    if(is_converged != 1){
      for(j=0;j<mp_fit->n_channels;j++){
	mp_fit->fit_data[j]->fit[i].status = RUNNING;
      }
    }

    /* 
     * If one peak has an error, mark them all as error and subtract them
     * out of the fit image.
     */
    if(is_bad){
      for(j=0;j<mp_fit->n_channels;j++){
	fit_data = mp_fit->fit_data[j];

	/* Subtract the peak out of the image. */
	fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);
	fit_data->fn_subtract_peak(fit_data);

	/* Set status to ERROR. */
	fit_data->fit[i].status = ERROR;
      }
    }    
  }
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

  if(VERBOSE){
    printf("mpNP %d\n", n_peaks);
  }
  
  for(i=0;i<mp_fit->n_channels;i++){
    j = i*mp_fit->nfit*NPEAKPAR;
    cfNewPeaks(mp_fit->fit_data[i], &peak_params[j], n_peaks);
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
  int i;
  double tmp;
  fitData *fit_data;
  
  for(i=0;i<mp_fit->n_channels;i++){
    fit_data = mp_fit->fit_data[i];
    
    /* Reset peak (it was changed by fn_update()), increase lambda. */
    tmp = fit_data->working_peak->lambda;
    fit_data->fn_copy_peak(&fit_data->fit[index], fit_data->working_peak);
    fit_data->working_peak->lambda = tmp * LAMBDAUP;

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
 *
 * The overall size is the number of channels times the spline
 * size in z.
 *
 * This cannot be called before mpInitializeChannel() as it needs
 * to know the spline size.
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
 * mpUpdate()
 *
 * Calculate weighted delta and update each channel.
 *
 * mp_fit->heights should be all 1.0 for fixed (relative) heights.
 *
 * Note this assumes that Spliner is using the following convention:
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = ZCENTER;
 *  delta[4] = BACKGROUND;
 */
void mpUpdate(mpFit *mp_fit)
{
  int i,nc,zi;
  double delta,p_ave,p_total,t,xoff,yoff;
  double *params_ch0,*heights;
  peakData *peak;
  fitData *fit_data_ch0;

  heights = mp_fit->heights;
  fit_data_ch0 = mp_fit->fit_data[0];
  params_ch0 = fit_data_ch0->working_peak->params;
  xoff = fit_data_ch0->xoff;
  yoff = fit_data_ch0->yoff;
  
  nc = mp_fit->n_channels;
  zi = ((splinePeak *)fit_data_ch0->working_peak->peak_model)->zi;
  
  /*
   * X parameters depends on the mapping.
   *
   * Note: The meaning of x and y is transposed here compared to in the
   *       mapping. This is also true for the y parameter below.
   */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(VERBOSE){
      printf(" x %d %.3e %.3e", i, heights[i], mp_fit->yt_Nto0[i*3+1]);
      printf(" %.3e %.3e %.3e\n", mp_fit->w_jacobian[i][2], mp_fit->yt_Nto0[i*3+2], mp_fit->w_jacobian[i][1]);
    }
    delta = mp_fit->yt_Nto0[i*3+1] * mp_fit->w_jacobian[i][2];
    delta += mp_fit->yt_Nto0[i*3+2] * mp_fit->w_jacobian[i][1];
    p_ave += delta * mp_fit->w_x[zi*nc+i] * heights[i];
    p_total += mp_fit->w_x[zi*nc+i] * heights[i];
  }
  delta = p_ave/p_total;
  mFitUpdateParam(fit_data_ch0->working_peak, delta, XCENTER);

  /* Y parameters also depend on the mapping. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(VERBOSE){
      printf(" y %d %.3e %.3e", i, heights[i], mp_fit->xt_Nto0[i*3+1]);
      printf(" %.3e %.3e %.3e\n", mp_fit->w_jacobian[i][2], mp_fit->xt_Nto0[i*3+2], mp_fit->w_jacobian[i][1]);
    }
    delta = mp_fit->xt_Nto0[i*3+1] * mp_fit->w_jacobian[i][2];
    delta += mp_fit->xt_Nto0[i*3+2] * mp_fit->w_jacobian[i][1];
    p_ave += delta * mp_fit->w_y[zi*nc+i] * heights[i];
    p_total += mp_fit->w_y[zi*nc+i] * heights[i];
  }
  delta = p_ave/p_total;
  mFitUpdateParam(fit_data_ch0->working_peak, delta, YCENTER);  

  /* 
   * Use mapping to update peak locations in the remaining channels.
   * 
   * Note: The meaning of x and y is transposed here compared to in the
   *       mapping.
   *
   * Note: Spliner uses the upper left corner as 0,0 so we need to adjust
   *       to the center, transform, then adjust back. This is particularly
   *       important if one channel is inverted relative to another.
   */
  for(i=1;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    
    t = mp_fit->yt_0toN[i*3];
    t += mp_fit->yt_0toN[i*3+1] * (params_ch0[YCENTER]+yoff);
    t += mp_fit->yt_0toN[i*3+2] * (params_ch0[XCENTER]+xoff);
    peak->params[XCENTER] = t-xoff;

    t = mp_fit->xt_0toN[i*3];
    t += mp_fit->xt_0toN[i*3+1] * (params_ch0[YCENTER]+yoff);
    t += mp_fit->xt_0toN[i*3+2] * (params_ch0[XCENTER]+xoff);
    peak->params[YCENTER] = t-yoff;
  }

  /* Update peak (integer) location with hysteresis. */
  for(i=0;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    if(fabs(peak->params[XCENTER] - (double)peak->xi - 0.5) > HYSTERESIS){
      peak->xi = (int)peak->params[XCENTER];
    }
    if(fabs(peak->params[YCENTER] - (double)peak->yi - 0.5) > HYSTERESIS){
      peak->yi = (int)peak->params[YCENTER];
    }
  }

  /* Z parameter is a simple weighted average. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    p_ave += mp_fit->w_jacobian[i][3] * mp_fit->w_z[zi*nc+i] * heights[i];
    p_total += mp_fit->w_z[zi*nc+i] * heights[i];
  }
  delta = p_ave/p_total;

  for(i=0;i<nc;i++){
    peak = mp_fit->fit_data[i]->working_peak;
    mFitUpdateParam(peak, delta, ZCENTER);

    /* Update (integer) z position. */
    ((splinePeak *)peak->peak_model)->zi = (int)(peak->params[ZCENTER]);
  }

  /* Backgrounds float independently. */
  for(i=0;i<nc;i++){
    mFitUpdateParam(mp_fit->fit_data[i]->working_peak, mp_fit->w_jacobian[i][4], BACKGROUND);
  }
}


/*
 * mpUpdateFixed()
 *
 * Calculate weighted delta and update each channel for fitting
 * with peak heights fixed relative to each other. This does not
 * change mp_fit->heights;
 *
 * Note: This allows negative heights, which will get removed by fn_check().
 *
 * Note: This assumes that Spliner is using the following convention:
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = ZCENTER;
 *  delta[4] = BACKGROUND;
 */
void mpUpdateFixed(mpFit *mp_fit)
{
  int i,nc,zi;
  double delta, p_ave, p_total;
  fitData *fit_data_ch0;

  fit_data_ch0 = mp_fit->fit_data[0];
  nc = mp_fit->n_channels;
  zi = ((splinePeak *)fit_data_ch0->working_peak->peak_model)->zi;

  /* Height, this is a simple weighted average. */
  p_ave = 0.0;
  p_total = 0.0;
  for(i=0;i<nc;i++){
    if(VERBOSE){
      printf(" h %d %.3e\n", i, mp_fit->w_jacobian[i][0]);
    }
    p_ave += mp_fit->w_jacobian[i][0] * mp_fit->w_h[zi*nc+i];
    p_total += mp_fit->w_h[zi*nc+i];
  }
  delta = p_ave/p_total;
  
  mFitUpdateParam(fit_data_ch0->working_peak, delta, HEIGHT);
  for(i=1;i<nc;i++){
    mp_fit->fit_data[i]->working_peak->params[HEIGHT] = fit_data_ch0->working_peak->params[HEIGHT];
  }

  mpUpdate(mp_fit);
}


/*
 * mpUpdateIndependent()
 *
 * Calculate weighted delta and update each channel for fitting
 * with independently adjustable peak heights.
 *
 * Note: This assumes that Spliner is using the following convention:
 *  delta[0] = HEIGHT;
 *  delta[1] = XCENTER;
 *  delta[2] = YCENTER;
 *  delta[3] = ZCENTER;
 *  delta[4] = BACKGROUND;
 */
void mpUpdateIndependent(mpFit *mp_fit)
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
}
