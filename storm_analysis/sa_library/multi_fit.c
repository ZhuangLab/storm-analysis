/*
 * Core routines that are common to all of the fitters.
 *
 * Hazen 10/16
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "multi_fit.h"


/*
 * mFitCalcErr()
 *
 * Calculate the fit error of the peak. Technically this is
 * actually the total error in the pixels that are covered
 * by the peak. When peaks overlap substantially they will
 * have similar errors.
 *
 * If the difference between the new and the old error is
 * sufficiently small this will also mark the peak as
 * converged.
 *
 * fit_data - pointer to a fitData structure.
 * peak_data - pointer to a peakData structure.
 */
void mFitCalcErr(fitData *fit_data, peakData *peak)
{
  int j,k,l,m;
  double err,fi,xi;

  if(peak->status == RUNNING){
    l = peak->yi * fit_data->image_size_x + peak->xi;
    err = 0.0;
    for(j=0;j<peak->size_y;j++){
      for(k=0;k<peak->size_x;k++){
	m = (j * fit_data->image_size_x) + k + l;
	fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
	if(fi <= 0.0){
	  if(TESTING){
	    printf(" Negative f detected! %.3f %d\n", fi, peak->index);
	  }
	  peak->status = ERROR;
	  fit_data->n_neg_fi++;
	  j = peak->size_y + 1;
	  k = peak->size_x + 1;
	}
	xi = fit_data->x_data[m];
	/* 
	 * This should not happen as the expectation is that negative image
	 * values are eliminated upstream of this step.
	 */
	if(TESTING){
	  if(xi <= 0.0){
	    printf(" Negative x detected! %.3f %d\n", xi, m);
	  }
	}
	err += 2*(fi-xi)-2*xi*log(fi/xi);
      }
    }
    peak->error_old = peak->error;
    peak->error = err;
    if (VERBOSE){
      printf("error: %d %f %f %f\n", peak->index, peak->error_old, peak->error, fit_data->tolerance);
    }
    if(((fabs(err - peak->error_old)/err) < fit_data->tolerance) && (peak->status != ERROR)){
      peak->status = CONVERGED;
    }
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
  for(i=0;i<7;i++){
    fit_data->clamp_start[i] = clamp[i];
  }

  /* Allocate space for image, fit and background arrays. */
  fit_data->bg_counts = (int *)malloc(sizeof(int)*im_size_x*im_size_y);
  fit_data->bg_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->f_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);
  fit_data->x_data = (double *)malloc(sizeof(double)*im_size_x*im_size_y);

  return fit_data;
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

  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->x_data[i] = new_image[i];
  }  
}


/*
 * mFitNewPeaks
 *
 * fit_data - Pointer to a fitData structure.
 */
void mFitNewPeaks(fitData *fit_data)
{
  int i;

  /*
   * Reset diagnostics.
   */
  fit_data->n_dposv = 0;
  fit_data->n_margin = 0;
  fit_data->n_neg_fi = 0;
  fit_data->n_neg_height = 0;
  fit_data->n_neg_width = 0;
 
  /*
   * Reset fitting arrays.
   */
  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->bg_counts[i] = 0;
    fit_data->bg_data[i] = 0.0;
    fit_data->f_data[i] = 0.0;
  }
}

  
/*
 * mFitUpdateParams
 *
 * Update peak parameters based on delta.
 */
void mFitUpdateParams(peakData *peak, double *delta)
{
  int i;

  if(VERBOSE){
    printf("%d : ", peak->index);  
  }
  for(i=0;i<NFITTING;i++){

    if(VERBOSE){
      printf("%.3e %.3f | ", delta[i], peak->clamp[i]);
    }
    
    if (delta[i] != 0.0){
      
      // update sign & clamp if the solution appears to be oscillating.
      if (peak->sign[i] != 0){
	if ((peak->sign[i] == 1) && (delta[i] < 0.0)){
	  peak->clamp[i] *= 0.5;
	}
	else if ((peak->sign[i] == -1) && (delta[i] > 0.0)){
	  peak->clamp[i] *= 0.5;
	}
      }
      if (delta[i] > 0.0){
	peak->sign[i] = 1;
      }
      else {
	peak->sign[i] = -1;
      }

      peak->params[i] -= delta[i]/(1.0 + fabs(delta[i])/peak->clamp[i]);
    }
  }
  if(VERBOSE){
    printf("\n");
  }
}


/*
 * The MIT License
 *
 * Copyright (c) 2013 Zhuang Lab, Harvard University
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
