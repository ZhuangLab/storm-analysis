/*
 * Core routines that are common to all of the fitters.
 *
 * Hazen 10/17
 */

/* Include */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "multi_fit.h"

/* LAPACK Functions */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);

/*
 * mFitAddPeak()
 *
 * Add working peak to the foreground and background data arrays. This
 * also adds the sCMOS variance term.
 *
 * fit_data - pointer to a fitData structure.
 */
void mFitAddPeak(fitData *fit_data)
{
  int j,k,l,m,n,o;
  double bg,mag;
  peakData *peak;

  peak = fit_data->working_peak;

  peak->added++;

  if(TESTING){
    if(peak->added != 1){
      printf("Peak count error detected in mFitAddPeak()! %d\n", peak->added);
      exit(EXIT_FAILURE);
    }
  }

  /* Add peak to the foreground and background arrays. */
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  mag = peak->params[HEIGHT];
  for(j=0;j<peak->size_y;j++){
    m = j*peak->size_x;
    n = j*fit_data->image_size_x;
    for(k=0;k<peak->size_x;k++){
      o = n + k + l;
      fit_data->f_data[o] += mag * peak->psf[m+k] * fit_data->rqe[o];
      fit_data->bg_counts[o] += 1;
      fit_data->bg_data[o] += bg + fit_data->scmos_term[o];
    }
  }
}

/*
 * mFitAnscombe()
 *
 * Calculate the Anscombe transform.
 *
 * x - Signal / function value.
 * var - (Gaussian) variance.
 */
double mFitAnscombe(double x)
{
  if(TESTING){
    if(x < (-0.375 + 1.0e-6)){
      printf(" Negative value in Anscombe transform! %.3f\n\n", x);
      exit(EXIT_FAILURE);
    }
  }
  return sqrt(x + 0.375);
}

/*
 * mFitCalcErr()
 *
 * Calculate the fit error of a peak. Technically this is actually the total error 
 * in the pixels that are covered by the peak. When peaks overlap substantially they 
 * will have similar errors.
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

  if(VERBOSE){
    printf("mFCE, xi - %d, yi - %d, sx - %d, sy - %d\n", peak->xi, peak->yi, peak->size_x, peak->size_y);
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
       */
      if(TESTING){
	if(xi <= 0.0){
	  printf(" Negative x detected!\n");
	  printf("   xi %.3f\n\n", xi);
	  err = peak->error;
	  j = peak->size_y + 1;
	  k = peak->size_x + 1;
	  exit(EXIT_FAILURE);
	}
      }
      err += 2.0*((fi-xi)-xi*log(fi/xi));
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
	  exit(EXIT_FAILURE);
	}
      }
    }
  }
  peak->error = err;
  
  return 0;
}


/*
 * mFitCalcErrALS()
 *
 * The Anscombe least-squares version of the error function.
 *
 * fit_data - pointer to a fitData structure.
 *
 * Returns 0 if there were no errors.
 */
int mFitCalcErrALS(fitData *fit_data)
{
  int j,k,l,m;
  double err,di,fi;
  peakData *peak;

  peak = fit_data->working_peak;

  if(peak->status != RUNNING){
    return 0;
  }

  if(VERBOSE){
    printf("mFCEL, xi - %d, yi - %d, sx - %d, sy - %d\n", peak->xi, peak->yi, peak->size_x, peak->size_y);
  }

  l = peak->yi * fit_data->image_size_x + peak->xi;
  err = 0.0;
  for(j=0;j<peak->size_y;j++){
    for(k=0;k<peak->size_x;k++){
      m = (j * fit_data->image_size_x) + k + l;
      
      /* The sCMOS variance term is already added into bg_data. */
      fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);

      /* 
       * Need to catch this because the fit background can be 
       * negative. In this case just use 0.0 as the value fi.
       */
      if(fi < (-0.375 + 1.0e-6)){
	if(TESTING){
	  printf(" Negative f detected!\n");
	  printf("  index %d\n", peak->index);
	  printf("      f %.3f\n", fi);
	  printf("    fit %.3f\n", fit_data->f_data[m]);
	  printf("     bg %.3f\n", fit_data->bg_data[m]);
	  printf("   cnts %d\n\n", fit_data->bg_counts[m]);
	}
	fit_data->n_neg_fi++;
	di = fit_data->a_data[m];
      }
      else{
	di = mFitAnscombe(fi) - fit_data->a_data[m];
      }
      err += di*di;
    }
  }
  peak->error = err;
  
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
   * Check that the peak hasn't moved to close to the edge of the image.
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
  int i;
  
  /* 
   * Free individual peaks. Freeing of the fitter specialized parts
   * are handled by the fitter.
   */
  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->max_nfit;i++){
      free(fit_data->fit[i].psf);
    }
    free(fit_data->fit);
  }

  free(fit_data->working_peak->psf);
  free(fit_data->working_peak);

  free(fit_data->a_data);
  free(fit_data->bg_counts);
  free(fit_data->bg_data);
  free(fit_data->bg_estimate);
  free(fit_data->f_data);
  free(fit_data->rqe);
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

  copy->added = original->added;
  copy->index = original->index;
  copy->iterations = original->iterations;
  copy->status = original->status;
  copy->xi = original->xi;
  copy->yi = original->yi;
  
  copy->size_x = original->size_x;
  copy->size_y = original->size_y;
  
  copy->error = original->error;

  copy->lambda = original->lambda;

  for(i=0;i<NFITTING;i++){
    copy->sign[i] = original->sign[i];
    copy->clamp[i] = original->clamp[i];
    copy->params[i] = original->params[i];
  }

  for(i=0;i<(original->size_x*original->size_y);i++){
    copy->psf[i] = original->psf[i];
  }
}


/*
 * mFitEstimatePeakHeight()
 *
 * Calculate the area under the peak of unit height and compare this to
 * the area under (image - current fit - estimated background) x peak.
 *
 * We are minimizing : fi * (h*fi - xi)^2
 * where fi is the peak shape at pixel i, h is the height and xi is
 * is the data (minus the current fit & estimated background.
 *
 * Taking the derivative with respect to h gives fi*fi*xi/fi*fi*fi as the
 * value for h that will minimize the above.
 */
void mFitEstimatePeakHeight(fitData *fit_data)
{
  int j,k,l,m,n,o;
  double sp,sx,t1;
  peakData *peak;

  peak = fit_data->working_peak;
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  sp = 0.0;  /* This is fi*fi*fi. */
  sx = 0.0;  /* This is fi*fi*xi. */
  for(j=0;j<peak->size_y;j++){
    m = j*peak->size_x;
    n = j*fit_data->image_size_x;
    for(k=0;k<peak->size_x;k++){
      o = n + k + l;
      t1 = peak->psf[m+k];
      sp += t1*t1*t1;
      sx += t1*t1*(fit_data->x_data[o] - fit_data->f_data[o] - fit_data->bg_estimate[o]);
    }
  }
  peak->params[HEIGHT] = sx/sp;
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
 * mFitGetNError()
 *
 * Return the number of fits that are in the error state.
 *
 * fit_data - Pointer to a fitData structure.
 */
int mFitGetNError(fitData *fit_data)
{
  int i,count;

  count = 0;
  for(i=0;i<fit_data->nfit;i++){
    if(fit_data->fit[i].status==ERROR){
      count++;
    }
  }

  return count;
}


/*
 * mFitGetPeakPropertyDouble()
 *
 * Return requested peak property (double).
 *
 * fit_data - Pointer to a fitData structure.
 * values - Pre-allocated storage for the results.
 * what - Which property to get.
 */
void mFitGetPeakPropertyDouble(fitData *fit_data, double *values, char *what)
{
  int i;
  
  if (!strcmp(what, "background")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].params[BACKGROUND];
    }
  }
  else if (!strcmp(what, "bg_sum")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = mFitPeakBgSum(fit_data, &fit_data->fit[i]);
    }
  }
  else if (!strcmp(what, "error")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].error;
    }
  }
  else if (!strcmp(what, "fg_sum")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = mFitPeakFgSum(fit_data, &fit_data->fit[i]);
    }
  }
  else if (!strcmp(what, "height")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].params[HEIGHT];
    }
  }
  else if (!strcmp(what, "sum")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = mFitPeakSum(&fit_data->fit[i]);
    }
  }
  else if (!strcmp(what, "x")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].params[XCENTER] + fit_data->xoff;
    }
  }
  else if (!strcmp(what, "xsigma")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = sqrt(1.0/(2.0*fit_data->fit[i].params[XWIDTH]));
    }
  }
  else if (!strcmp(what, "y")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].params[YCENTER] + fit_data->yoff;
    }
  }
  else if (!strcmp(what, "ysigma")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = sqrt(1.0/(2.0*fit_data->fit[i].params[YWIDTH]));
    }
  }
  else if (!strcmp(what, "z")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].params[ZCENTER] + fit_data->zoff;
    }
  }
  else{
    printf("Unrecognized parameter '%s'!\n", what);
  }
}


/*
 * mFitGetPeakPropertyInt()
 *
 * Return requested peak values as integers.
 *
 * fit_data - Pointer to a fitData structure.
 * values - Pre-allocated storage for the results.
 * what - Which property to get.
 */
void mFitGetPeakPropertyInt(fitData *fit_data, int32_t *values, char *what)
{
  int i;

  if (!strcmp(what, "iterations")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].iterations;
    }
  }
  else if (!strcmp(what, "status")){
    for(i=0;i<fit_data->nfit;i++){
      values[i] = fit_data->fit[i].status;
    }
  }
  else{
    printf("Unrecognized parameter '%s'!\n", what);
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
 * rqe - Pixel relative quantum efficiency.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * clamp - The starting clamp values for each peak.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* mFitInitialize(double *rqe, double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y)
{
  int i;
  fitData* fit_data;
  
  /* Initialize fitData structure. */
  fit_data = (fitData*)malloc(sizeof(fitData));

  fit_data->n_dposv = 0;
  fit_data->n_iterations = 0;
  fit_data->n_lost = 0;
  fit_data->n_margin = 0;
  fit_data->n_neg_fi = 0;
  fit_data->n_neg_height = 0;
  fit_data->n_neg_width = 0;
  fit_data->n_non_converged = 0;
  fit_data->n_non_decr = 0;

  fit_data->max_nfit = 0;
  fit_data->nfit = 0;
  fit_data->image_size_x = im_size_x;
  fit_data->image_size_y = im_size_y;
  fit_data->tolerance = tol;

  /* 
   * The default behavior is to immediately ERROR out peaks that start 
   * with a negative height. However this is a problem for multi-plane
   * analysis where some peaks in a group could have negative heights 
   * due to noise and large z values, so multi-plane sets this to a 
   * small positive value.
   */
  fit_data->minimum_height = -1.0;
  
  fit_data->xoff = 0.0;
  fit_data->yoff = 0.0;
  fit_data->zoff = 0.0;

  /* Copy RQE data. */
  fit_data->rqe = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  printf("RQE[0] %.3f\n", rqe[0]);
  for(i=0;i<(im_size_x*im_size_y);i++){
    fit_data->rqe[i] = rqe[i];
  }

  /* Copy sCMOS calibration data. */
  fit_data->scmos_term = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  printf("sCMOS[0] %.3f\n", scmos_calibration[0]);
  for(i=0;i<(im_size_x*im_size_y);i++){
    fit_data->scmos_term[i] = scmos_calibration[i];
  }

  /* Copy starting clamp values. */
  for(i=0;i<NFITTING;i++){
    fit_data->clamp_start[i] = clamp[i];
  }

  /* Allocate space for image, fit and background arrays. */
  fit_data->a_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->bg_counts = (int *)malloc(sizeof(int)*im_size_x*im_size_y);
  fit_data->bg_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->bg_estimate = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->f_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->x_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);

  /* Allocate space for the working peak. */
  fit_data->working_peak = (peakData *)malloc(sizeof(peakData));
  fit_data->working_peak->psf = NULL;
  fit_data->working_peak->peak_model = NULL;

  fit_data->fit = NULL;
  fit_data->fit_model = NULL;

  fit_data->fn_alloc_peaks = NULL;
  fit_data->fn_calc_JH = NULL;
  fit_data->fn_calc_peak_shape = NULL;
  fit_data->fn_check = NULL;
  fit_data->fn_copy_peak = NULL;
  fit_data->fn_error_fn = mFitCalcErr;
  fit_data->fn_free_peaks = NULL;
  fit_data->fn_update = NULL;
  
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

  double starting_error;               /* Initial error value for the peak. */
  
  double jacobian[NFITTING];           /* Jacobian */
  double w_jacobian[NFITTING];         /* Working copy of the Jacobian. */
  double hessian[NFITTING*NFITTING];   /* Hessian */
  double w_hessian[NFITTING*NFITTING]; /* Working copy of the Hessian. */

  if(VERBOSE){
    printf("mFILM\n");
  }

  for(i=0;i<fit_data->nfit;i++){

    if(VERBOSE){
      printf("\nmFILM peak - %d\n", i);
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

    /* 
     * Calculate initial error.
     *
     * Why? This might have changed from the previous cycle because the peak
     * background value could be shifted by neighboring peaks, creating a 
     * situation where it is impossible to improve on the old error value.
     */
    /* 
     * FIXME: We're ignoring the return value here and just assuming that
     *        the error function will work. Maybe okay because when we 
     *        calculate the updated error we check the return value? Not
     *        sure what would be the correct thing to do anyway.
     */
    fit_data->fn_error_fn(fit_data);
    starting_error = fit_data->working_peak->error;

    /* Calculate Jacobian and Hessian. This is expected to use 'working_peak'. */
    fit_data->fn_calc_JH(fit_data, jacobian, hessian);
    
    /* Subtract working peak out of image. */
    mFitSubtractPeak(fit_data);
    n_add--;

    j = 0;
    while(1){
      j++;

      if(VERBOSE){
	printf("  cycle %d %d\n", j, n_add);
      }

      /* Check if we are stuck on this peak, error it out if we are. */
      if(fit_data->working_peak->lambda > LAMBDAMAX){
	fit_data->n_lost++;
	fit_data->working_peak->status = ERROR;
	break;
      }
	
      /* Update peak iterations counter. */
      fit_data->working_peak->iterations++;
      
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
        fit_data->working_peak->lambda = fit_data->working_peak->lambda * LAMBDAUP;
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
	mFitResetPeak(fit_data, i);

	/* Set status to ERROR in case this is the last iteration. */
	fit_data->working_peak->status = ERROR;
	
	continue;	
      }
      
      /* Add peak 'working_peak' back to fit image. */
      fit_data->fn_calc_peak_shape(fit_data);
      mFitAddPeak(fit_data);
      n_add++;

      /* 
       * Calculate error for 'working_peak' with the new parameters. This
       * will also check if the fit has converged. It will return 0 if
       * everything is okay (the fit image has no negative values).
       */
      if(fit_data->fn_error_fn(fit_data)){
	if(VERBOSE){
	  printf(" Error calculation failed\n");
	}
	/* Subtract 'working_peak' from the fit image. */
	mFitSubtractPeak(fit_data);
	n_add--;
	 
	/* 
	 * Try again with a larger lambda. We need to reset the peak 
	 * state because fn_update() and fn_add_peak() changed it.
	 */
	mFitResetPeak(fit_data, i);

	/* Set status to ERROR in case this is the last iteration. */
	fit_data->working_peak->status = ERROR;
	
	continue;	
      }
      
      /* Check whether the error improved. */
      if(fit_data->working_peak->error > starting_error){
	if(VERBOSE){
	  printf("    increasing error %.6e %.6e %.6e\n", fit_data->working_peak->error, starting_error, fit_data->working_peak->lambda);
	}

	/* 
	 * Check for error convergence. 
	 *
	 * Usually this will happen because the lambda term has gotten so 
	 * large that the peak will barely move in the update.
	 */
      	if (((fit_data->working_peak->error - starting_error)/starting_error) < fit_data->tolerance){
	  fit_data->working_peak->status = CONVERGED;
	  break;
	}
	else{
	  fit_data->n_non_decr++;
	  
	  /* Subtract 'working_peak' from the fit image. */
	  mFitSubtractPeak(fit_data);
	  n_add--;

	  /* 
	   * Try again with a larger lambda. We need to reset the 
	   * peak state because fn_update() changed it.
	   */
	  mFitResetPeak(fit_data, i);
	}
      }
      else{
	
	if(VERBOSE){
	  printf("    decreasing error %.6e %.6e %.6e\n", fit_data->working_peak->error, starting_error, fit_data->working_peak->lambda);
	}

	/* Check for error convergence. */
      	if (((starting_error - fit_data->working_peak->error)/starting_error) < fit_data->tolerance){
	  fit_data->working_peak->status = CONVERGED;
	}
	
	/* Decrease lambda and exit the while loop. */
	else if(fit_data->working_peak->lambda > LAMBDAMIN){
	  fit_data->working_peak->lambda = LAMBDADOWN * fit_data->working_peak->lambda;
	}
	break;
      }
    }

    /* We expect n_add to be 1 if there were no errors, 0 otherwise. */
    if(TESTING){
      if(fit_data->working_peak->status == ERROR){
	if(n_add != 0){
	  printf("Problem detected in peak addition / subtraction logic, status == ERROR, counts = %d\n", n_add);
	  exit(EXIT_FAILURE);
	}
      }
      else{
	if(n_add != 1){
	  printf("Problem detected in peak addition / subtraction logic, status != ERROR, counts = %d\n", n_add);
	  exit(EXIT_FAILURE);
	}
      }
    }
    
    /* Copy updated working peak back into current peak. */
    fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
  }

  /* Recenter peaks as necessary. */
  mFitRecenterPeaks(fit_data);
}


/*
 * mFitIterateOriginal
 *
 * FIXME: Fix bit-rot or remove this fitting approach.
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
    mFitSubtractPeak(fit_data);
    
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
    fit_data->fn_calc_peak_shape(fit_data);
    mFitAddPeak(fit_data);
    
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
     * Note: This replicates the logic flaw in the original 3D-DAOSTORM
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
 * mFitNewBackground
 *
 * Copy in a new estimate of the background.
 *
 * fit_data - Pointer to a fitData structure.
 * background - Pointer to the background data of size image_size_x by image_size_y.
 */
void mFitNewBackground(fitData *fit_data, double *background)
{
  int i;

  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->bg_estimate[i] = background[i];
  }  
}


/*
 * mFitNewImage
 *
 * Copy in a new image to fit. We also reset everything that is fitting related
 * as this call indicates the start of a new cycle of fitting.
 *
 * fit_data - Pointer to a fitData structure.
 * new_image - Pointer to the image data of size image_size_x by image_size_y.
 * use_als - Doing Anscombe least-squares fitting.
 */
void mFitNewImage(fitData *fit_data, double *new_image, int use_als)
{
  int i;

  if(VERBOSE){
    printf("mFNI\n");
  }

  /* Copy the image & add scmos term (variance / gain * gain). */
  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->x_data[i] = new_image[i] + fit_data->scmos_term[i];
  }
  
  if (use_als != 0){
    if (VERBOSE){
      printf("Using ALS fit\n");
    }
    for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
      fit_data->a_data[i] = mFitAnscombe(fit_data->x_data[i]);
    }
  }

  /* Reset fitting arrays. */
  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->bg_counts[i] = 0;
    fit_data->bg_data[i] = 0.0;
    fit_data->f_data[i] = 0.0;
  }

  fit_data->nfit = 0;
}


/*
 * mFitNewPeaks
 *
 * These are new peaks to add to our current list of peaks.
 * 
 * fit_data - Pointer to a fitData structure.
 * n_peaks - The number of peaks.
 */
void mFitNewPeaks(fitData *fit_data, int n_peaks)
{
  int i,j,n_alloc,start,stop;
  peakData *peak,*new_peaks;

  if(VERBOSE){
    printf("mFNP %d\n", n_peaks);
  }
  
  /* 1. Check if we need more storage. */
  if ((fit_data->nfit + n_peaks) > fit_data->max_nfit){

    n_alloc = INCNPEAKS*((fit_data->nfit + n_peaks)/INCNPEAKS + 1);

    new_peaks = (peakData *)malloc(sizeof(peakData)*n_alloc);
    for(j=0;j<n_alloc;j++){
      new_peaks[j].psf = NULL;
    }
    fit_data->fn_alloc_peaks(new_peaks, n_alloc);

    i = 0;
    for(j=0;j<fit_data->nfit;j++){
      if(fit_data->fit[j].status != ERROR){
	fit_data->fn_copy_peak(&fit_data->fit[j], &new_peaks[i]);
	i += 1;

	/* Check that this peak is in the image. */
	if(TESTING){
	  if(fit_data->fit[j].added == 0){
	    printf("Peak %d is not in the image.\n", j);
	    exit(EXIT_FAILURE);
	  }
	}
      }
      else{
	/* Check that this peak is not in the image. */
	if(TESTING){
	  if(fit_data->fit[j].added != 0){
	    printf("Peak %d is in error state, but still in the image.\n", j);
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }
    
    /* Free old peak storage (if necessary). */
    if(fit_data->fit != NULL){
      fit_data->fn_free_peaks(fit_data->fit, fit_data->max_nfit);
      for(j=0;j<fit_data->max_nfit;j++){
	free(fit_data->fit[j].psf);
      }
      free(fit_data->fit);
    }

    /* Point to new peak storage and update counters. */
    fit_data->fit = new_peaks;
    fit_data->max_nfit = n_alloc;
    fit_data->nfit = i;
  }
  
  /* 2. Generic peak initialization. */
  start = fit_data->nfit;
  stop = fit_data->nfit + n_peaks;
  for(i=start;i<stop;i++){
    peak = &fit_data->fit[i];
    peak->added = 0;
    peak->index = i;
    peak->iterations = 0;

    /* Initial status. */
    peak->status = RUNNING;

    /* Initial error values. */
    peak->error = 0.0;

    /* Initial lambda value. */
    peak->lambda = LAMBDASTART;
  }
  
  /* 3. Caller must update the value of fit_data->nfit! */
}


/*
 * mFitPeakBgSum
 *
 * Return the background sum for the purposes of peak significance
 * calculations.
 *
 * This subtracts the fit shape off of the image. Then it multiplies
 * what is left by the PSF squared. The result is the number of
 * photo-electrons in the background if the background was filtered
 * by the PSF that describes the peak. Assuming Poisson statistics this
 * is also the variance. It is slightly complicated by the fact that
 * the PSF is not normalized.
 *
 * This is the same algorithm that is used in the Python peak finder
 * module for image segmentation.
 */
double mFitPeakBgSum(fitData *fit_data, peakData *peak)
{
  int j,k,l,m,n,o;
  double bg,bg_sum,fg,mag,psf,psf_sum;

  bg_sum = 0.0;
  psf_sum = 0.0;
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  mag = peak->params[HEIGHT];
  for(j=0;j<peak->size_y;j++){
    m = j*peak->size_x;
    n = j*fit_data->image_size_x;
    for(k=0;k<peak->size_x;k++){
      o = n + k + l;

      psf = peak->psf[m+k];
      psf_sum += psf;
      fg = mag * psf;

      /*
       * The background is the sum of the residual after subtracting the
       * foreground fit, and then convolved with the PSF.
       */
      bg = (fit_data->x_data[o] - fit_data->scmos_term[o]) - fg;
      bg_sum += bg*psf*psf;
    }
  }

  if((bg_sum > 0.0)&&(psf_sum > 0.0)){
    
    /* Correct bg sum calculation for the PSF not being normalized. */
    bg_sum = bg_sum/(psf_sum*psf_sum);
  }
  else{
    bg_sum = 1.0;
    if(TESTING){
      printf("0 or negative background sum or PSF sum in peak background sum calculation.\n");
    }
  }
  
  return bg_sum;
}


/*
 * mFitPeakFgSum
 *
 * Return the foreground sum for the purposes of peak significance
 * calculations.
 *
 * This convolves the fit peak shape with the PSF, mirroring the
 * calculation that is done to determine the background sum. It is
 * also slightly complicated by the fact that the PSF is not 
 * normalized.
 */
double mFitPeakFgSum(fitData *fit_data, peakData *peak)
{
  int j,k,l;
  double fg_sum,psf,psf_sum;

  fg_sum = 0.0;
  psf_sum = 0.0;
  
  for(j=0;j<peak->size_y;j++){
    l = j*peak->size_x;
    for(k=0;k<peak->size_x;k++){
      psf = peak->psf[l+k];
      psf_sum += psf;
      fg_sum += psf*psf;
    }
  }
  fg_sum = peak->params[HEIGHT]*fg_sum;

  if(psf_sum > 0.0){
    
    /* Correct foreground sum calculation for the PSF not being normalized. */
    fg_sum = fg_sum/psf_sum;
  }
  else{
    fg_sum = 1.0;
    if(TESTING){
      printf("0 or negative psf sum in peak foreground sum calculation.\n");
    }
  }
  
  return fg_sum;
}


/*
 * mFitPeakSum()
 *
 * Return the integrated intensity of the requested peak.
 *
 * peak - the peak to sum.
 */
double mFitPeakSum(peakData *peak)
{
  int i;
  double sum;

  sum = 0.0;
  for(i=0;i<(peak->size_y*peak->size_x);i++){
    sum += peak->psf[i];
  }
  sum = sum * peak->params[HEIGHT];
  
  return sum;
}


/*
 * mFitRecenterPeaks()
 *
 * Recenter the ROIs of any RUNNING peaks that are not centered. This is
 * called after LM fitting to move the ROIs of peaks whose center has
 * moved too far from the ROI center. We don't force the center to be
 * exactly inside the pixel (i.e. between 0.0 - 1.0) as this causes
 * problems for peaks that are near the pixel boundary.
 */
void mFitRecenterPeaks(fitData *fit_data)
{
  int i;
  double dx,dy;
  peakData *peak;

  for(i=0;i<fit_data->nfit;i++){

    if(fit_data->fit[i].status != RUNNING){
      continue;
    }

    peak = &fit_data->fit[i];

    dx = peak->params[XCENTER] - peak->xi;
    dy = peak->params[YCENTER] - peak->yi;
    
    if((dx<-0.25)||(dx>1.25)||(dy<-0.25)||(dy>1.25)){

      /* Copy into working peak. */
      fit_data->fn_copy_peak(peak, fit_data->working_peak);

      /* Subtract peak at current ROI. */
      mFitSubtractPeak(fit_data);

      /* Move the ROI. */
      fit_data->working_peak->xi = (int)floor(peak->params[XCENTER]);
      fit_data->working_peak->yi = (int)floor(peak->params[YCENTER]);

      /* 
       * Check that the ROI is still in the image. Mark it for removal if
       * it is not. 
       */
      if(mFitCheck(fit_data)){
	fit_data->n_lost += 1;
	fit_data->working_peak->status = ERROR;
      }
      else{
	fit_data->fn_calc_peak_shape(fit_data);
	mFitAddPeak(fit_data);
      }

      /* Copy working peak back. */
      fit_data->fn_copy_peak(fit_data->working_peak, peak);      
    }
  }
}

/*
 * mFitRemoveErrorPeaks()
 *
 * This removes all the peaks that are in the error state from the
 * peaks array. If none of the peaks are in the error state this
 * this is basically a NOP. If some are then there is a bit of
 * copying as it overwrites the error peaks with good peaks, shortening
 * the list of peaks in the process.
 *
 * The idea is that this is called after fitting and before you start
 * filtering out adjacent peaks, etc. in Python.
 */
void mFitRemoveErrorPeaks(fitData *fit_data)
{
  int i,j;
  
  i = 0;
  for(j=0;j<fit_data->nfit;j++){
    if(fit_data->fit[j].status != ERROR){
      if(j!=i){
	fit_data->fn_copy_peak(&fit_data->fit[j], &fit_data->fit[i]);
      }
      i += 1;

      /* Check that this peak is in the image. */
      if(TESTING){
	if(fit_data->fit[j].added == 0){
	  printf("Peak %d is not in the image.\n", j);
	  exit(EXIT_FAILURE);
	}
      }
    }
    else{
      /* Check that this peak is not in the image. */
      if(TESTING){
	if(fit_data->fit[j].added > 0){
	  printf("Peak %d is in error state, but still in the image (%d).\n", j, fit_data->fit[j].added);
	  exit(EXIT_FAILURE);
	}
      }
    }
  }
  fit_data->nfit = i;
}


/*
 * mFitRemoveRunningPeaks()
 *
 * This removes all the peaks that are in the running state from the
 * peaks array. If none of the peaks are in the running state this
 * this is basically a NOP. If some are then there is a bit of
 * copying as it overwrites the running peaks with converged peaks, 
 * shortening the list of peaks in the process.
 *
 * The idea is that this is called at the end of the analysis before
 * you save the results.
 */
void mFitRemoveRunningPeaks(fitData *fit_data)
{
  int i,j;
  
  i = 0;
  for(j=0;j<fit_data->nfit;j++){
    if(fit_data->fit[j].status != RUNNING){
      if(j!=i){
	fit_data->fn_copy_peak(&fit_data->fit[j], &fit_data->fit[i]);
      }
      i += 1;

      /* 
       * Check that this peak is not in the ERROR state. Peaks that
       * are in the error state should have been removed before this
       * function was called.
       */
      if(TESTING){
	if(fit_data->fit[j].status == ERROR){
	  printf("Peak %d is in the error state!\n", j);
	  exit(EXIT_FAILURE);
	}
      }
    }
    else{
      fit_data->n_non_converged++;
    }
  }
  fit_data->nfit = i;
}


/*
 * mFitResetClampValues()
 *
 * Resets the clamp values for all of the peaks.
 */
void mFitResetClampValues(fitData *fit_data)
{
  int i,j;
  int bg_warning;
  peakData *peak;

  bg_warning = 0;
  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];

    /* Initial clamp values. */
    for(j=0;j<NFITTING;j++){
      peak->clamp[j] = fit_data->clamp_start[j];
      peak->sign[j] = 0;
    }

    /* Height and background clamp values are relative. */
    peak->clamp[HEIGHT] = fit_data->clamp_start[HEIGHT]*peak->params[HEIGHT];
    peak->clamp[BACKGROUND] = fit_data->clamp_start[BACKGROUND]*peak->params[BACKGROUND];
    
    /*
     * Print a warning if the background is approximately zero and we are using the 
     * clamp in fitting. This is more likely than height being zero..
     */
    if ((peak->clamp[BACKGROUND] < 1.0e-3) && !bg_warning){
      printf("Warning! Background clamp is zero due to peak with zero background!\n");
      bg_warning = 1;
    }
  }
}


/*
 * mFitResetPeak()
 *
 * This is used during fitting to restore the working peak to it's previous
 * state, but with increased lambda.
 */
void mFitResetPeak(fitData *fit_data, int index)
{
  int tmp_added;
  double tmp_lambda;

  tmp_added = fit_data->working_peak->added;
  tmp_lambda = fit_data->working_peak->lambda;
  fit_data->fn_copy_peak(&fit_data->fit[index], fit_data->working_peak);
  fit_data->working_peak->added = tmp_added;
  fit_data->working_peak->lambda = tmp_lambda * LAMBDAUP;
}


/*
 * mFitSetPeakStatus()
 *
 * Set the status of the peaks.
 *
 * fit_data - Pointer to a fitData structure.
 * status - Integer array containing new status values.
 */
void mFitSetPeakStatus(fitData *fit_data, int32_t *status)
{
  int i;

  for(i=0;i<fit_data->nfit;i++){   
    /* 
     * If we marked this peak as ERROR and it was in the image we 
     * need to subtract it from the image.
     */
    if(status[i] == ERROR){
      if(VERBOSE){
	printf(" mFSPS %d %d\n", i, status[i]);
      }
      if(fit_data->fit[i].added > 0){
	fit_data->fn_copy_peak(&fit_data->fit[i], fit_data->working_peak);
	mFitSubtractPeak(fit_data);
	fit_data->fn_copy_peak(fit_data->working_peak, &fit_data->fit[i]);
      }
    }

    fit_data->fit[i].status = status[i];
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
 * mFitSubtractPeak()
 *
 * Subtract working peak out of the current fit, basically 
 * this just undoes addPeak().
 *
 * fit_data - pointer to a fitData structure.
 */
void mFitSubtractPeak(fitData *fit_data)
{
  int j,k,l,m,n,o;
  double bg,mag;
  peakData *peak;

  peak = fit_data->working_peak;

  peak->added--;

  if(TESTING){
    if(peak->added != 0){
      printf("Peak count error detected in mFitSubtractPeak()! %d\n", peak->added);
      exit(EXIT_FAILURE);
    }
  }
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  mag = peak->params[HEIGHT];
  for(j=0;j<peak->size_y;j++){
    m = j*peak->size_x;
    n = j*fit_data->image_size_x;
    for(k=0;k<peak->size_x;k++){
      o = n + k + l;
      fit_data->f_data[o] -= mag * peak->psf[m+k] * fit_data->rqe[o];
      fit_data->bg_counts[o] -= 1;
      fit_data->bg_data[o] -= (bg + fit_data->scmos_term[o]);
    }
  }
}


/*
 * mFitUpdate()
 *
 * Updates working_peak (integer) center. This is called in the 
 * LM loop. It moves the pixel location of the peak.
 *
 * peak - pointer to a peakData structure.
 */
void mFitUpdate(peakData *peak)
{
  double dx,dy;
  
  dx = peak->params[XCENTER] - peak->xi;
  if((dx<-0.5)||(dx>1.5)){
    peak->xi = (int)floor(peak->params[XCENTER]);
  }

  dy = peak->params[YCENTER] - peak->yi;
  if((dy<-0.5)||(dy>1.5)){
    peak->yi = (int)floor(peak->params[YCENTER]);
  }
}


/*
 * mFitUpdateParam
 *
 * Update peak parameter based on delta.
 */
void mFitUpdateParam(peakData *peak, double delta, int i)
{
  if(VERBOSE){
    if(USECLAMP){
      printf("mFUP %d : %d %.3e %.3e %.3e\n", peak->index, i, peak->params[i], delta, peak->clamp[i]);
    }
    else{
      printf("mFUP %d : %d %.3e %.3e\n", peak->index, i, peak->params[i], delta);
    }
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
