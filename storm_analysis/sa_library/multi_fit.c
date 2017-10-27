/*
 * Core routines that are common to all of the fitters.
 *
 * Hazen 10/17
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "multi_fit.h"

/* LAPACK Functions */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);

/*
 * mFitCalcErr()
 *
 * Calculate the fit error of working_peak. Technically this is actually the 
 * total error in the pixels that are covered by the peak. When peaks overlap 
 * substantially they will have similar errors.
 *
 * If the difference between the new and the old error is sufficiently small 
 * this will also mark the peak as converged.
 *
 * fit_data - pointer to a fitData structure.
 *
 * Returns 0 if there were no errors.
 */
int mFitCalcErr(fitData *fit_data)
{
  int j,k,l,m;
  double err,fi,xi;
  peakData *peak;
  
  peak = fit_data->working_peak;

  if(peak->status != RUNNING){
    return 0;
  }
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  err = 0.0;
  for(j=0;j<peak->size_y;j++){
    for(k=0;k<peak->size_x;k++){
      m = (j * fit_data->image_size_x) + k + l;
      fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
      if(fi <= 0.0){
	/*
	 * This can happen because the fit background can be negative. I
	 * don't think it is a problem that merits crashing everything.
	 */
	if(TESTING){
	  printf(" Negative f detected!\n");
	  printf("  index %d\n", peak->index);
	  printf("      f %.3f\n", fi);
	  printf("    fit %.3f\n", fit_data->f_data[m]);
	  printf("     bg %.3f\n", fit_data->bg_data[m]);
	  printf("   cnts %d\n\n", fit_data->bg_counts[m]);
	}
	fit_data->n_neg_fi++;
	return 1;
      }
      xi = fit_data->x_data[m];
      /* 
       * This should not happen as the expectation is that negative image
       * values are eliminated upstream of this step.
       *
       * It would probably be better if we just crashed in this situation
       * as it indicates that the image was not correctly processed?
       */
      if(TESTING){
	if(xi <= 0.0){
	  printf(" Negative x detected!\n");
	  printf("   xi %.3f\n\n", xi);
	  err = peak->error;
	  j = peak->size_y + 1;
	  k = peak->size_x + 1;
	}
      }
      err += 2*(fi-xi)-2*xi*log(fi/xi);
      if(TESTING){
	/*
	 * FIXME: Should also test for +- infinity?
	 */
	if (isnan(err)){
	  printf(" NAN error detected!\n");
	  printf("  index %d\n", peak->index);
	  printf("     fi %.3f\n", fi);
	  printf("     xi %.3f\n\n", xi);
	  j = peak->size_y + 1;
	  k = peak->size_x + 1;
	}
      }
    }
  }
  peak->error_old = peak->error;
  peak->error = err;
  if (VERBOSE){
    printf("error: %d %.2f %.2f %.2e\n", peak->index, peak->error_old, peak->error, fit_data->tolerance);
  }
  if(((fabs(err - peak->error_old)/err) < fit_data->tolerance) && (peak->status != ERROR)){
    peak->status = CONVERGED;
  }

  return 0;
}


/*
 * mFitCheck()
 *
 * Check that the parameters of working_peak are still valid.
 *
 * fit_data - pointer to a fitData structure.
 *
 * Return 0 if okay.
 */
int mFitCheck(fitData *fit_data)
{
  int xi,yi;
  peakData *peak;

  peak = fit_data->working_peak;  
  
  /*
   * Check that the peak hasn't moved to close to the 
   * edge of the image. Flag the peak as bad if it has.
   */
  xi = peak->xi;
  yi = peak->yi;
  if((xi < 0)||(xi >= (fit_data->image_size_x - peak->size_x))||(yi < 0)||(yi >= (fit_data->image_size_y - peak->size_y))){
    fit_data->n_margin++;
    if(TESTING){
      printf("object outside margins, %d, %d, %d\n", peak->index, xi, yi);
    }
    return 1;
  }
  
  /* 
   * Check for negative height. 
   */
  if(peak->params[HEIGHT] < 0.0){
    fit_data->n_neg_height++;
    if(TESTING){
      printf("negative height, %d, %.3f\n", peak->index, peak->params[HEIGHT]);
    }
    return 1;
  }

  return 0;
}

/*
 * mfitCleanup()
 *
 * Frees the fitData structure.
 *
 * fit_data - pointer to a fitData structure.
 */
void mFitCleanup(fitData *fit_data)
{
  free(fit_data->working_peak);
  free(fit_data->bg_counts);
  free(fit_data->bg_data);
  free(fit_data->f_data);
  free(fit_data->scmos_term);
  free(fit_data->x_data);
  free(fit_data);
}


/*
 * mFitCopyPeak()
 *
 * Copy the contents of a peakData structure into another peakData 
 * structure. Note that the peak_model pointer is not copied as doing
 * this is the responsibility of the particular instantiation of the
 * the fitter.
 *
 * original - pointer to a peakData structure.
 * copy - pointer to a peakData structure.
 */
void mFitCopyPeak(peakData *original, peakData *copy)
{
  int i;
  
  copy->index = original->index;
  copy->status = original->status;
  copy->xi = original->xi;
  copy->yi = original->yi;
  
  copy->size_x = original->size_x;
  copy->size_y = original->size_y;
  
  copy->error = original->error;
  copy->error_old = original->error_old;

  copy->lambda = original->lambda;

  for(i=0;i<NFITTING;i++){
    copy->sign[i] = original->sign[i];
    copy->clamp[i] = original->clamp[i];
    copy->params[i] = original->params[i];
  }
}


/*
 * mFitGetFitImage()
 *
 * Return an image created from the current best fit peaks.
 */
void mFitGetFitImage(fitData *fit_data, double *fit_image)
{
  int i;
  
  for(i=0;i<(fit_data->image_size_x * fit_data->image_size_y);i++){
    fit_image[i] = fit_data->f_data[i];
  }
}


/*
 * mFitGetResidual(residual).
 *
 * Returns image - fit.
 *
 * fit_data - Pointer to a fitData structure.
 * residual - Pre-allocated space to store the residual values.
 *            This should be square & the same size as the image.
 */
void mFitGetResidual(fitData *fit_data, double *residual)
{
  int i;

  //calcFit(fit_data);
  for(i=0;i<(fit_data->image_size_x * fit_data->image_size_y);i++){
    residual[i] = fit_data->x_data[i] - fit_data->f_data[i];
  }
}


/*
 * getResults(params)
 *
 * Return the current fitting results.
 *
 * fit_data - Pointer to a fitData structure.
 * peak_params - pre-allocated space for storing the peak fitting parameters.
 *            1. height
 *            2. x-center
 *            3. x-sigma
 *            4. y-center
 *            5. y-sigma
 *            6. background
 *            7. z-center
 *            8. status
 *            9. fit error
 *
 *   This array should be n * number of peaks passed in to
 *   initialize the fitting.
 *
 */
void mFitGetResults(fitData *fit_data, double *peak_params)
{
  int i;
  peakData *peak;

  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];

    if(peak->status != ERROR){
      peak_params[i*NPEAKPAR+XWIDTH] = sqrt(1.0/(2.0*peak->params[XWIDTH]));
      peak_params[i*NPEAKPAR+YWIDTH] = sqrt(1.0/(2.0*peak->params[YWIDTH]));
    }
    else{
      peak_params[i*NPEAKPAR+XWIDTH] = 1.0;
      peak_params[i*NPEAKPAR+YWIDTH] = 1.0;
    }
    peak_params[i*NPEAKPAR+HEIGHT]     = peak->params[HEIGHT];
    peak_params[i*NPEAKPAR+XCENTER]    = peak->params[XCENTER] + fit_data->xoff;
    peak_params[i*NPEAKPAR+YCENTER]    = peak->params[YCENTER] + fit_data->yoff;
    peak_params[i*NPEAKPAR+BACKGROUND] = peak->params[BACKGROUND];
    peak_params[i*NPEAKPAR+ZCENTER]    = peak->params[ZCENTER] + fit_data->zoff;

    peak_params[i*NPEAKPAR+STATUS] = (double)peak->status;
    peak_params[i*NPEAKPAR+IERROR] = peak->error;
  }
}


/*
 * mFitGetUnconverged()
 *
 * Return the number of fits that have not yet converged.
 *
 * fit_data - Pointer to a fitData structure.
 */
int mFitGetUnconverged(fitData *fit_data)
{
  int i,count;

  count = 0;
  for(i=0;i<fit_data->nfit;i++){
    if(fit_data->fit[i].status==RUNNING){
      count++;
    }
  }

  return count;
}


/*
 * mFitInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * clamp - The starting clamp values for each peak.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* mFitInitialize(double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y)
{
  int i;
  fitData* fit_data;
  
  /* Initialize fitData structure. */
  fit_data = (fitData*)malloc(sizeof(fitData));

  fit_data->n_dposv = 0;
  fit_data->n_iterations = 0;
  fit_data->n_margin = 0;
  fit_data->n_neg_fi = 0;
  fit_data->n_neg_height = 0;
  fit_data->n_neg_width = 0;
  fit_data->n_non_decr = 0;
  
  fit_data->image_size_x = im_size_x;
  fit_data->image_size_y = im_size_y;
  fit_data->tolerance = tol;
  fit_data->fit = NULL;

  fit_data->xoff = 0.0;
  fit_data->yoff = 0.0;
  fit_data->zoff = 0.0;

  /* Copy sCMOS calibration data. */
  fit_data->scmos_term = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  for(i=0;i<(im_size_x*im_size_y);i++){
    fit_data->scmos_term[i] = scmos_calibration[i];
  }

  /* Copy starting clamp values. */
  for(i=0;i<NFITTING;i++){
    fit_data->clamp_start[i] = clamp[i];
  }

  /* Allocate space for image, fit and background arrays. */
  fit_data->bg_counts = (int *)malloc(sizeof(int)*im_size_x*im_size_y);
  fit_data->bg_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->f_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->x_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);

  /* Allocate space for the working peak. */
  fit_data->working_peak = (peakData *)malloc(sizeof(peakData));
    
  return fit_data;
}


/*
 * mFitIterateLM
 *
 * Perform a single iteration of fitting update for each peaks.
 *
 * fit_data - Pointer to a fitData structure.
 */
void mFitIterateLM(fitData *fit_data)
{
  int i,j,k,l,m;
  int info;
  int n_add;

  double tmp;
  
  double jacobian[NFITTING];           /* Jacobian */
  double w_jacobian[NFITTING];         /* Working copy of the Jacobian. */
  double hessian[NFITTING*NFITTING];   /* Hessian */
  double w_hessian[NFITTING*NFITTING]; /* Working copy of the Hessian. */

  if(VERBOSE){
    printf("mFILM\n");
  }
  
  for(i=0;i<fit_data->nfit;i++){

    if(VERBOSE){
      printf("mFILM %d\n", i);
    }
    
    /* 
     * This is for debugging, to make sure that we not adding more times than
     * we are subtracting. 
     */
    n_add = 1;

    /* Skip ahead if this peak is not RUNNING. */
    if(fit_data->fit[i].status != RUNNING){
      continue;
    }

    /* Copy current peak into working peak. */
    fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);

    /* Calculate Jacobian and Hessian. This is expected to use 'working_peak'. */
    fit_data->fn_calc_JH(fit_data, jacobian, hessian);
    
    /* Subtract current peak out of image. This is expected to use 'working_peak'. */
    fit_data->fn_subtract_peak(fit_data);
    n_add--;

    for(j=0;j<=MAXCYCLES;j++){

      if(VERBOSE){
	printf("  cycle %d %d\n", j, n_add);
      }
      
      /* Update total fitting iterations counter. */
      fit_data->n_iterations++;

      /*
       * Reset status flag. We only started this loop if the peak was RUNNING.
       * However the status might have been changed to error due to a Cholesky
       * solver issue or invalid peak parameters in a previous iteration of
       * this loop.
       */
      fit_data->working_peak->status = RUNNING;

      /* Copy Jacobian and Hessian. */
      for(k=0;k<fit_data->jac_size;k++){
	w_jacobian[k] = jacobian[k];
	m = k*fit_data->jac_size;
	for(l=0;l<fit_data->jac_size;l++){
	  if (k == l){
	    w_hessian[m+l] = (1.0 + fit_data->working_peak->lambda) * hessian[m+l];
	  }
	  else{
	    w_hessian[m+l] = hessian[m+l];
	  }
	}
      }
      
      /* 
       * Solve for update. Note that this also changes w_jacobian
       * which is one of the reasons why we made a copy.
       */
      info = mFitSolve(w_hessian, w_jacobian, fit_data->jac_size);
    
      if(info!=0){
	if(VERBOSE){
	  printf(" mFitSolve() failed %d\n", info);
	}
	fit_data->n_dposv++;
	fit_data->working_peak->status = ERROR;
	
	/* If the solver failed, try again with a larger lambda. */
        fit_data->working_peak->lambda = fit_data->working_peak->lambda * 2.0;
	continue;
      }
      
      /* Update 'working_peak'. mFitSolve returns the update in w_jacobian. */
      fit_data->fn_update(fit_data, w_jacobian);

      /* 
       * Check that it is still in the image, etc.. The fn_check function
       * should return 0 if everything is okay.
       */
      if(fit_data->fn_check(fit_data)){
	if(VERBOSE){
	  printf(" fn_check() failed\n");
	}
	/* 
	 * Try again with a larger lambda. We need to reset the 
	 * peak state because fn_update() changed it.
	 */
	tmp = fit_data->working_peak->lambda;
	fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);
	fit_data->working_peak->lambda = 2.0 * tmp;

	/* Set status to ERROR in case this is the last iteration. */
	fit_data->working_peak->status = ERROR;
	
	continue;	
      }
      
      /* Add peak 'working_peak' back to fit image. */
      fit_data->fn_add_peak(fit_data);
      n_add++;

      /* 
       * Calculate error for 'working_peak' with the new parameters. This
       * will also check if the fit has converged. It will return 0 if
       * everything is okay (the fit image has no negative values).
       */
      if(mFitCalcErr(fit_data)){
	if(VERBOSE){
	  printf(" mFitCalcErr() failed\n");
	}
	/* Subtract 'working_peak' from the fit image. */
	fit_data->fn_subtract_peak(fit_data);
	n_add--;
	 
	/* 
	 * Try again with a larger lambda. We need to reset the peak 
	 * state because fn_update() and fn_add_peak() changed it.
	 */
	tmp = fit_data->working_peak->lambda;
	fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);
	fit_data->working_peak->lambda = 2.0 * tmp;

	/* Set status to ERROR in case this is the last iteration. */
	fit_data->working_peak->status = ERROR;
	
	continue;	
      }

      /* 
       * Break if not still running, mFitCalcErr may have decided that the
       * peak had converged so we don't need to do anything more.
       */
      if(fit_data->working_peak->status != RUNNING){
	break;
      }
      
      /* Check whether the error improved. */
      if(fit_data->working_peak->error > fit_data->working_peak->error_old){

	/* 
	 * If we have reached the maximum number of iterations, then 
	 * the peak stays where it is and we hope for the best.
	 */
	if(j<MAXCYCLES){
	  fit_data->n_non_decr++;
	  
	  /* Subtract 'working_peak' from the fit image. */
	  fit_data->fn_subtract_peak(fit_data);
	  n_add--;

	  /* 
	   * Try again with a larger lambda. We need to reset the 
	   * peak state because fn_update() changed it.
	   */
	  tmp = fit_data->working_peak->lambda;
	  fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);
	  fit_data->working_peak->lambda = 2.0 * tmp;
	}

	if(TESTING){
	  if(j==MAXCYCLES){
	    printf("Reached max cycles with no improvement in peak error for %d\n", i);
	  }
	}
      }
      else{
	/* Decrease lambda and exit the for loop. */
	fit_data->working_peak->lambda = 0.5 * fit_data->working_peak->lambda;
	break;
      }
    }

    /* We expect n_add to be 1 if there were no errors, 0 otherwise. */
    if(TESTING){
      if(fit_data->working_peak->status == ERROR){
	if(n_add != 0){
	  printf("Problem detected in peak addition / subtraction logic, status == ERROR, counts = %d\n", n_add);
	}
      }
      else{
	if(n_add != 1){
	  printf("Problem detected in peak addition / subtraction logic, status != ERROR, counts = %d\n", n_add);
	}
      }
    }
    
    /* Copy updated working peak back into current peak. */
    fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
  }
}


/*
 * mFitIterateOriginal
 *
 * Perform a single iteration of fitting update for each peak
 * using the original 3D-DAOSTORM algorithm.
 *
 * fit_data - Pointer to a fitData structure.
 */
void mFitIterateOriginal(fitData *fit_data)
{
  int i;
  int info;
  
  double jacobian[NFITTING];           /* Jacobian */
  double hessian[NFITTING*NFITTING];   /* Hessian */

  if(VERBOSE){
    printf("mFIO\n");
  }

  if(!USECLAMP){
    printf("Warning! mFitIterateOriginal() without clamping. Mistake?\n");
  }
  
  /*
   * 1. Calculate updated peaks.
   */
  for(i=0;i<fit_data->nfit;i++){

    if(VERBOSE){
      printf("mFIO %d\n", i);
    }

    /* Skip ahead if this peak is not RUNNING. */
    if(fit_data->fit[i].status != RUNNING){
      continue;
    }

    /* Copy current peak into working peak. */
    fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);

    /* Calculate Jacobian and Hessian. This is expected to use 'working_peak'. */
    fit_data->fn_calc_JH(fit_data, jacobian, hessian);
    
    /* Subtract current peak out of image. This is expected to use 'working_peak'. */
    fit_data->fn_subtract_peak(fit_data);
    
    /* Update total fitting iterations counter. */
    fit_data->n_iterations++;

    /* 
     * Solve for update. Note that this also changes jacobian.
     */
    info = mFitSolve(hessian, jacobian, fit_data->jac_size);

    /* If the solver failed, drop this peak from the analysis. */
    if(info!=0){
      if(VERBOSE){
	printf(" mFitSolve() failed %d\n", info);
      }
      fit_data->n_dposv++;
      fit_data->working_peak->status = ERROR;
      fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
      continue;
    }
      
    /* Update 'working_peak'. mFitSolve returns the update in w_jacobian. */
    fit_data->fn_update(fit_data, jacobian);

    /* 
     * Check that it is still in the image, etc.. The fn_check function
     * should return 0 if everything is okay.
     */
    if(fit_data->fn_check(fit_data)){
      if(VERBOSE){
	printf(" fn_check() failed\n");
      }
      
      /* Set status to ERROR and drop this peak from the analysis. */
      fit_data->working_peak->status = ERROR;
      fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
      continue;	
    }
      
    /* Add peak 'working_peak' back to fit image. */
    fit_data->fn_add_peak(fit_data);
        
    /* Copy updated working peak back into current peak. */
    fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
  }


  /*
   * 2) Calculate peak errors.
   */
  for(i=0;i<fit_data->nfit;i++){

    if(VERBOSE){
      printf("mFIO %d\n", i);
    }
    
    /* Skip ahead if this peak is not RUNNING. */
    if(fit_data->fit[i].status != RUNNING){
      continue;
    }
  
    /* 
     * Calculate error for 'working_peak' with the new parameters. This
     * will also check if the fit has converged. It will return 0 if
     * everything is okay (the fit image has no negative values).
     *
     * Note: The replicates the logic flaw in the original 3D-DAOSTORM
     *       algorithm, if the peak had negative fi it was not discarded
     *       but remained present in a somewhat ambiguous state.
     */
    fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);
    if(mFitCalcErr(fit_data)){
      if(VERBOSE){
	printf(" mFitCalcErr() failed\n");
      }
    }
    fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
  }
}


/*
 * mFitNewImage
 *
 * Copy in a new image to fit.
 *
 * fit_data - Pointer to a fitData structure.
 * new_image - Pointer to the image data of size image_size_x by image_size_y.
 */
void mFitNewImage(fitData *fit_data, double *new_image)
{
  int i;

  if(VERBOSE){
    printf("mFNI\n");
  }

  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->x_data[i] = new_image[i];
  }  
}


/*
 * mFitNewPeaks
 *
 * fit_data - Pointer to a fitData structure.
 */
void mFitNewPeaks(fitData *fit_data, double *peak_params, int n_peaks)
{
  int i,j,bg_warning;
  peakData *peak;

  if(VERBOSE){
    printf("mFNP %d\n", n_peaks);
  }
  
  /*
   * Reset fitting arrays.
   */
  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->bg_counts[i] = 0;
    fit_data->bg_data[i] = 0.0;
    fit_data->f_data[i] = 0.0;
  }

  /*
   * Initialize peaks (localizations).
   */
  bg_warning = 0;
  fit_data->nfit = n_peaks;
  fit_data->fit = (peakData *)malloc(sizeof(peakData)*n_peaks);
  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    peak->index = i;

    /* Initial status. */
    peak->status = (int)(peak_params[i*NPEAKPAR+STATUS]);
    if(peak->status==RUNNING){
      peak->error = 0.0;
      peak->error_old = 0.0;
    }
    else {
      peak->error = peak_params[i*NPEAKPAR+IERROR];
      peak->error_old = peak->error;
    }

    /* Initial lambda value. */
    peak->lambda = LAMBDASTART;
    
    /* Initial clamp values. */
    for(j=0;j<NFITTING;j++){
      peak->clamp[j] = fit_data->clamp_start[j];
      peak->sign[j] = 0;
    }

    /* Height and background clamp values are relative. */
    peak->clamp[HEIGHT] = fit_data->clamp_start[HEIGHT]*peak_params[i*NPEAKPAR+HEIGHT];
    peak->clamp[BACKGROUND] = fit_data->clamp_start[BACKGROUND]*peak_params[i*NPEAKPAR+BACKGROUND];

    /* Print a warning if the background is approximately zero and we are using the clamp in fitting. */
    if(USECLAMP){
      if ((peak->clamp[BACKGROUND] < 1.0e-3) && !bg_warning){
	printf("Warning! Background clamp is zero due to peak with zero background!\n");
	bg_warning = 1;
      }
    }
  }
}

  
/*
 * mFitSolve
 *
 * Solve for update vector given jacobian and hessian.
 */
int mFitSolve(double *hessian, double *jacobian, int p_size)
{
  int i,j;
  
  // Lapack
  int n, nrhs = 1, lda, ldb, info;

  n = p_size;
  lda = p_size;
  ldb = p_size;

  // Use Lapack to solve AX=B to calculate update vector
  dposv_("Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info);

  if(VERBOSE){
    if(info!=0){
      printf(" dposv_ failed with %d\n", info);
      for(i=0;i<p_size;i++){
	printf("%.3f\t", jacobian[i]);
      }
      printf("\n\n");
      for(i=0;i<p_size;i++){
	for(j=0;j<p_size;j++){
	  printf("%.3f\t", hessian[i*p_size+j]);
	}
	printf("\n");
      }
      printf("\n");      
    }
  }
  
  return info;
}


/*
 * mFitUpdateParam
 *
 * Update peak parameter based on delta.
 */
void mFitUpdateParam(peakData *peak, double delta, int i)
{
  if(VERBOSE){
    printf("mFUP %d : %d %.3e %.3e %.3e\n", peak->index, i, peak->params[i], delta, peak->clamp[i]);
  }

  /* With clamping. */
  if(USECLAMP){
    if (delta != 0.0){
      
      // update sign & clamp if the solution appears to be oscillating.
      if (peak->sign[i] != 0){
	if ((peak->sign[i] == 1) && (delta < 0.0)){
	  peak->clamp[i] *= 0.5;
	}
	else if ((peak->sign[i] == -1) && (delta > 0.0)){
	  peak->clamp[i] *= 0.5;
	}
      }
      if (delta > 0.0){
	peak->sign[i] = 1;
      }
      else {
	peak->sign[i] = -1;
      }
      
      peak->params[i] -= delta/(1.0 + fabs(delta)/peak->clamp[i]);
    }
  }
  /* No clamping. */
  else{
    peak->params[i] -= delta;
  }
}


/*
 * The MIT License
 *
 * Copyright (c) 2017 Zhuang Lab, Harvard University
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
