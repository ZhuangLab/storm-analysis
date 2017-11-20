/*
 * Fit multiple, possible overlapping, cubic splines simultaneously to
 * image data.
 *
 * Hazen 01/14
 */

/* Include */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cubic_fit.h"


/*
 * cfAddPeak()
 *
 * Calculate peak shape and add the working peak to the 
 * foreground and background data arrays.
 *
 * fit_data - pointer to a fitData structure.
 */
void cfAddPeak(fitData *fit_data)
{
  int j,k,l,m,psx,psy;
  double bg,height;
  peakData *peak;
  splinePeak *spline_peak;

  peak = fit_data->working_peak;
  spline_peak = (splinePeak *)peak->peak_model;

  peak->added++;
  
  if(TESTING){
    if(peak->added != 1){
      printf("Peak count error detected in cfAddPeak()! %d\n", peak->added);
      exit(EXIT_FAILURE);
    }
  }  
  
  psx = peak->size_x;
  psy = peak->size_y;
  
  /* Add peak to the foreground and background arrays. */
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  height = peak->params[HEIGHT];
  for (j=0;j<psy;j++){
    for (k=0;k<psx;k++){
      m = j * fit_data->image_size_x + k + l;
      fit_data->f_data[m] += height*spline_peak->peak_values[j*psx+k];
      fit_data->bg_counts[m] += 1;
      fit_data->bg_data[m] += bg + fit_data->scmos_term[m];
    }
  }
}


/*
 * cfAllocPeaks()
 *
 * Allocate storage for cfPeaks. Note that this does not allocate
 * space for the peak_values element, which is done in cfNewPeaks().
 */
struct peakData *cfAllocPeaks(int n_peaks)
{
  int i;
  peakData *new_peaks;

  new_peaks = (peakData *)malloc(sizeof(peakData)*n_peaks);  
  for(i=0;i<n_peaks;i++){
    new_peaks[i].peak_model = (splinePeak *)malloc(sizeof(splinePeak));
  }
  return new_peaks;
}


/* 
 * cfCalcJH2D()
 *
 * Calculate Jacobian and Hessian for the 2D spline.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void cfCalcJH2D(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m,n;
  int x_start, y_start;
  double height,fi,t1,t2,xi;
  double jt[4];
  peakData *peak;
  splinePeak *spline_peak;
  splineFit *spline_fit;

  /* Initializations. */
  peak = fit_data->working_peak;
  spline_peak = (splinePeak *)peak->peak_model;
  spline_fit = (splineFit *)fit_data->fit_model;
  x_start = spline_peak->x_start;
  y_start = spline_peak->y_start;
    
  for(i=0;i<4;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<16;i++){
    hessian[i] = 0.0;
  }

  /* Calculate values x, y, xx, xy, yy, etc. terms for a 2D spline. */
  computeDelta2D(spline_fit->spline_data, spline_peak->y_delta, spline_peak->x_delta);

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
      jt[1] = -0.5*height*dxfAt2D(spline_fit->spline_data,2*j+y_start,2*k+x_start);
      jt[2] = -0.5*height*dyfAt2D(spline_fit->spline_data,2*j+y_start,2*k+x_start);
      jt[3] = 1.0;

      /* Calculate jacobian. */
      t1 = 2.0*(1.0 - xi/fi);
      for(m=0;m<4;m++){
	jacobian[m] += t1*jt[m];
      }
	  
      /* Calculate hessian. */
      t2 = 2.0*xi/(fi*fi);
      for(m=0;m<4;m++){
	for(n=m;n<4;n++){
	  hessian[m*4+n] += t2*jt[m]*jt[n];
	}
      }
    }
  }
}


/* 
 * cfCalcJH3D()
 *
 * Calculate Jacobian and Hessian for the 3D spline.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void cfCalcJH3D(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m,n,zi;
  int x_start, y_start;
  double height,fi,t1,t2,xi;
  double jt[5];
  peakData *peak;
  splinePeak *spline_peak;
  splineFit *spline_fit;

  /* Initializations. */
  peak = fit_data->working_peak;
  spline_peak = (splinePeak *)peak->peak_model;
  spline_fit = (splineFit *)fit_data->fit_model;
  
  x_start = spline_peak->x_start;
  y_start = spline_peak->y_start;
  zi = spline_peak->zi;

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
}


/*
 * cfCalcPeakShape()
 *
 * Calculate peak shape.
 *
 * fit_data - pointer to a fitData structure.
 */
void cfCalcPeakShape(fitData *fit_data)
{
  int j,k,psx,psy,x_start,y_start;
  double xd,yd,zd;
  peakData *peak;
  splinePeak *spline_peak;
  splineFit *spline_fit;

  peak = fit_data->working_peak;
  spline_peak = (splinePeak *)peak->peak_model;
  spline_fit = (splineFit *)fit_data->fit_model;

  psx = peak->size_x;
  psy = peak->size_y;
  
  /* 
   * Calculate peak shape using the cubic_spline library.
   *
   * The spline resolution is 2x that of the analyzed image in x and y, 
   * so there is some fiddling here to evaluate the spline at the correct 
   * location.
   *
   * Basically xd (short for xdelta), yd and zd are all expected to be 
   * in the range 0.0 - 1.0. These values are the same for each unit
   * cell of the spline, so as an optimization we are only calculating 
   * them and the values that depend on them (in the cubic_spline library) 
   * once. This is more complicated for xd and yd as due to the 2x size 
   * difference, we need to figure out whether we want the spline to be
   * evaluated in the first or the second half pixel.
   */
  xd = 2.0*(2.0 - (peak->params[XCENTER] - peak->xi));
  yd = 2.0*(2.0 - (peak->params[YCENTER] - peak->yi));
  zd = (peak->params[ZCENTER] - spline_peak->zi);

  x_start = 0;
  while(xd>1.0){
    x_start += 1;
    xd -= 1.0;
  }

  y_start = 0;
  while(yd>1.0){
    y_start += 1;
    yd -= 1.0;
  }

  /*
   * Save the values for convenience as we'll need them again when 
   * calculating how to update the fit.
   */
  spline_peak->x_start = x_start;
  spline_peak->y_start = y_start;
  spline_peak->x_delta = xd;
  spline_peak->y_delta = yd;
  spline_peak->z_delta = zd;

  /* 
   * FIXME: These functions are also calculating the values necessary for
   *        calculating derivatives. A possible optimization would be to
   *        have a flag so that this does not happen.
   */
  if(spline_fit->fit_type == S3D){
    computeDelta3D(spline_fit->spline_data, zd, yd, xd);
    for(j=0;j<psy;j++){
      for(k=0;k<psx;k++){
	spline_peak->peak_values[j*psx+k] = fAt3D(spline_fit->spline_data,spline_peak->zi,2*j+y_start,2*k+x_start);
      }
    }
  }
  else{
    computeDelta2D(spline_fit->spline_data, yd, xd);
    for(j=0;j<psy;j++){
      for(k=0;k<psx;k++){
	spline_peak->peak_values[j*psx+k] = fAt2D(spline_fit->spline_data,2*j+y_start,2*k+x_start);
      }
    }    
  }
}


/*
 * cfCleanup()
 *
 * Frees the fitData structure.
 *
 * fit_data - pointer to a fitData structure.
 */
void cfCleanup(fitData *fit_data)
{
  int i;
  splinePeak *spline_peak;
  splineFit *spline_fit;

  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->nfit;i++){
      spline_peak = (splinePeak *)(fit_data->fit[i].peak_model);
      free(spline_peak->peak_values);
      free(spline_peak);
    }
    free(fit_data->fit);
  }
  
  spline_peak = (splinePeak *)(fit_data->working_peak->peak_model);
  free(spline_peak->peak_values);
  free(spline_peak);
    
  spline_fit = (splineFit *)fit_data->fit_model;
  splineCleanup(spline_fit->spline_data);

  mFitCleanup(fit_data);
}


/*
 * cfCopyPeak()
 *
 * Copies the contents of peak structure into another peak structure.
 *
 * original - pointer to a peakData structure.
 * copy - pointer to a peakData structure.
 */
void cfCopyPeak(peakData *original, peakData *copy)
{
  int i;
  splinePeak *spline_copy, *spline_original;

  spline_copy = (splinePeak *)copy->peak_model;
  spline_original = (splinePeak *)original->peak_model;

  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(original, copy);

  /* Copy the parts that are specific to Spliner. */
  spline_copy->zi = spline_original->zi;

  spline_copy->x_start = spline_original->x_start;
  spline_copy->y_start = spline_original->y_start;

  spline_copy->x_delta = spline_original->x_delta;
  spline_copy->y_delta = spline_original->y_delta;
  spline_copy->z_delta = spline_original->z_delta;

  for(i=0;i<(copy->size_x*copy->size_y);i++){
    spline_copy->peak_values[i] = spline_original->peak_values[i];
  }
}


/*
 * cfFreePeaks()
 *
 * Frees a peakData array.
 *
 * peaks - Pointer to an array of peakData.
 * n_peaks - The size of the array.
 */
void cfFreePeaks(peakData *peaks, int n_peaks)
{
  int i;
  splinePeak *spline_peak;

  for(i=0;i<n_peaks;i++){
    spline_peak = (splinePeak *)(peaks[i].peak_model);
    free(spline_peak->peak_values);
    free(spline_peak);
  }
  free(peaks);
}


/*
 * cfInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * spline_data - Pointer to a splineData structure.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * clamp - The starting clamp values for each peak.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* cfInitialize(splineData *spline_data, double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y)
{
  int sx,sy;
  fitData* fit_data;

  fit_data = mFitInitialize(scmos_calibration, clamp, tol, im_size_x, im_size_y);

  /*
   * Initialize fit model.
   */
  fit_data->fit_model = (splineFit *)malloc(sizeof(splineFit));
  ((splineFit *)fit_data->fit_model)->fit_type = spline_data->type;
  ((splineFit *)fit_data->fit_model)->spline_size_z = spline_data->zsize;
  ((splineFit *)fit_data->fit_model)->spline_data = spline_data;

  /*
   * Calculate offset to the spline center as XCENTER and YCENTER
   * in spline fitting correspond to the location of the top left
   * corner of the spline. (In the assumed coordinate system where 
   * x = 0, y = 0 is the upper left hand corner of the image).
   *
   * Notes:
   *  1. Spline sizes are 2n + 1. 
   *  2. This assumes that the splines are square.
   */
  sx = spline_data->xsize/2 - 1;
  sy = spline_data->ysize/2 - 1;

  if((sx%2)==1){
    fit_data->xoff = (double)(sx/2) - 1.0;
    fit_data->yoff = (double)(sy/2) - 1.0;
  }
  else{
    fit_data->xoff = (double)(sx/2) - 1.5;
    fit_data->yoff = (double)(sy/2) - 1.5;
  }

  /*
   * Save spline size in pixels in x and y for later convenience.
   */
  ((splineFit *)fit_data->fit_model)->spline_size_x = sx;
  ((splineFit *)fit_data->fit_model)->spline_size_y = sy;
  
  /* Allocate storage for the working peak. */
  fit_data->working_peak->peak_model = (splinePeak *)malloc(sizeof(splinePeak));
  ((splinePeak *)fit_data->working_peak->peak_model)->peak_values = (double *)malloc(sizeof(double)*sx*sy);

  /* Set function pointers. */
  fit_data->fn_add_peak = &cfAddPeak;
  fit_data->fn_alloc_peaks = &cfAllocPeaks;
  fit_data->fn_calc_peak_shape = &cfCalcPeakShape;
  fit_data->fn_check = &mFitCheck;
  fit_data->fn_copy_peak = &cfCopyPeak;
  fit_data->fn_free_peaks = &cfFreePeaks;
  fit_data->fn_subtract_peak = &cfSubtractPeak;  
  
  return fit_data;
}


/*
 * cfInitialize2D()
 *
 * Initializes 2D spline fitting.
 *
 * fit_data - pointer to a fitData structure.
 */
void cfInitialize2D(fitData *fit_data)
{
  fit_data->jac_size = 4;
  
  fit_data->fn_calc_JH = &cfCalcJH2D;
  fit_data->fn_update = &cfUpdate2D;
}


/*
 * cfInitialize3D()
 *
 * Initializes 3D spline fitting.
 *
 * fit_data - pointer to a fitData structure.
 */
void cfInitialize3D(fitData *fit_data)
{
  fit_data->jac_size = 5;
  
  fit_data->fn_calc_JH = &cfCalcJH3D;
  fit_data->fn_update = &cfUpdate3D;
}


/*
 * cfNewPeaks
 *
 * fit_data - Pointer to a fitData structure.
 * peak_params - Input values for the peak parameters.
 * n_peaks - The number of peaks.
 */
void cfNewPeaks(fitData *fit_data, double *peak_params, char *p_type, int n_peaks)
{
  int i,j,k,l,m,n,xc,yc;
  int start,stop;
  double sp,sx,t1;
  peakData *peak;
  splinePeak *spline_peak;

  if(VERBOSE){
    printf("cfNP %d\n", n_peaks);
  }

  /* Generic initializations. */
  mFitNewPeaks(fit_data, n_peaks);

  /* Spliner specific initializations. */
  start = fit_data->nfit;
  stop = fit_data->nfit + n_peaks;
  
  /*
   * 'finder' or 'testing' parameters, these are the peak x,y and z 
   * values as an n_peaks x 3 array.
   */
  if(!strcmp(p_type, "finder") || !strcmp(p_type, "testing")){
    for(i=start;i<stop;i++){
      j = 3*(i-start);
      
      peak = &fit_data->fit[i];
      spline_peak = (splinePeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      
      /*
       * Note: Even though these are the same for every peak (as the spline
       *       does not change size during fitting), they are duplicated
       *       for each peak for the benefit of the mFit functions.
       */
      peak->size_x = ((splineFit *)fit_data->fit_model)->spline_size_x;
      peak->size_y = ((splineFit *)fit_data->fit_model)->spline_size_y;

      /* Allocate space for saving the peak shape. */
      spline_peak->peak_values = (double *)malloc(sizeof(double)*peak->size_x*peak->size_y);

      /* Calculate (integer) peak locations. */
      peak->xi = (int)round(peak->params[XCENTER]);
      peak->yi = (int)round(peak->params[YCENTER]);
      spline_peak->zi = (int)peak->params[ZCENTER];

      /* Estimate background. */
      xc = (int)round(peak_params[j]);
      yc = (int)round(peak_params[j+1]);      
      peak->params[BACKGROUND] = fit_data->bg_estimate[yc * fit_data->image_size_x + xc];

      /* Arbitrary initial value for HEIGHT. */
      peak->params[HEIGHT] = 1.0;      

      /* Copy into working peak. */
      cfCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	cfCopyPeak(fit_data->working_peak, peak);
	continue;
      }

      if(!strcmp(p_type, "finder")){

	/* Calculate peak shape (of working peak). */
	cfCalcPeakShape(fit_data);

	/* 
	 * Estimate height. 
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
	spline_peak = (splinePeak *)fit_data->working_peak->peak_model;

	k = peak->yi * fit_data->image_size_x + peak->xi;
	sp = 0.0;  /* This is fi*fi*fi. */
	sx = 0.0;  /* This is fi*fi*xi. */
	for(l=0;l<peak->size_y;l++){
	  for(m=0;m<peak->size_x;m++){
	    n = l * fit_data->image_size_x + m + k;
	    t1 = spline_peak->peak_values[l*peak->size_x+m];
	    sp += t1*t1*t1;
	    sx += t1*t1*(fit_data->x_data[n] - fit_data->f_data[n] - fit_data->bg_estimate[n]);
	  }
	}
	fit_data->working_peak->params[HEIGHT] = sx/sp;
	
	if(fit_data->working_peak->params[HEIGHT] < fit_data->minimum_height){
	  fit_data->working_peak->params[HEIGHT] = fit_data->minimum_height;
	}

	/* Check that the initial height is positive, error it out if not. */
	if(fit_data->working_peak->params[HEIGHT] <= 0.0){
	  printf("Warning peak %d has negative estimated height!\n", (i-start));
	  fit_data->working_peak->status = ERROR;
	  cfCopyPeak(fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      cfAddPeak(fit_data);

      /* Copy values back from working peak. */
      cfCopyPeak(fit_data->working_peak, peak);      
    }
  }
  /*
   * "pre-specified" parameters, these are the peak x, y, z, background
   * and height as an n_peaks x 5 array.
   */
  else if(!strcmp(p_type, "text") || !strcmp(p_type, "insight3")){
    for(i=start;i<stop;i++){
      j = 5*(i-start);
      peak = &fit_data->fit[i];
      spline_peak = (splinePeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      peak->params[BACKGROUND] = peak_params[j+3];
      peak->params[HEIGHT] = peak_params[j+4];
      
      /*
       * Note: Even though these are the same for every peak (as the spline
       *       does not change size during fitting), they are duplicated
       *       for each peak for the benefit of the mFit functions.
       */
      peak->size_x = ((splineFit *)fit_data->fit_model)->spline_size_x;
      peak->size_y = ((splineFit *)fit_data->fit_model)->spline_size_y;

      /* Allocate space for saving the peak shape. */
      spline_peak->peak_values = (double *)malloc(sizeof(double)*peak->size_x*peak->size_y);

      /* Calculate (integer) peak locations. */
      peak->xi = (int)round(peak->params[XCENTER]);
      peak->yi = (int)round(peak->params[YCENTER]);
      spline_peak->zi = (int)peak->params[ZCENTER];

      /* Copy into working peak. */
      cfCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	cfCopyPeak(fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      cfCalcPeakShape(fit_data);
      cfAddPeak(fit_data);

      /* Copy values back from working peak. */
      cfCopyPeak(fit_data->working_peak, peak);
    }
  }
  else{
    printf("Unrecognized peak type '%s'!\n", p_type);
  }  

  /*
   * Initial error calculation. 
   */
  for(i=start;i<stop;i++){
    peak = &fit_data->fit[i];
    cfCopyPeak(peak, fit_data->working_peak);
    mFitCalcErr(fit_data);
    cfCopyPeak(fit_data->working_peak, peak);
  }

  fit_data->nfit = stop;

  /* Reset the clamp values on all the peaks. */
  if(USECLAMP){
    mFitResetClampValues(fit_data);
  }  
}


/*
 * cfSubtractPeak()
 *
 * Subtract the working peak out of the current fit, basically 
 * this just undoes addPeak().
 *
 * fit_data - pointer to a fitData structure.
 */
void cfSubtractPeak(fitData *fit_data)
{
  int j,k,l,m;
  double bg,height;
  peakData *peak;
  splinePeak *spline_peak;

  peak = fit_data->working_peak;
  spline_peak = (splinePeak *)peak->peak_model;

  peak->added--;

  if(TESTING){
    if(peak->added != 0){
      printf("Peak count error detected in cfSubtractPeak()! %d\n", peak->added);
      exit(EXIT_FAILURE);
    }
  }
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  height = peak->params[HEIGHT];
  for (j=0;j<peak->size_y;j++){
    for (k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      fit_data->f_data[m] -= height*spline_peak->peak_values[j*peak->size_x+k];
      fit_data->bg_counts[m] -= 1;
      fit_data->bg_data[m] -= (bg + fit_data->scmos_term[m]);
    }
  }
}


/*
 * cfUpdate()
 *
 * Updates working_peak location with hysteresis.
 *
 * peak - pointer to a peakData structure.
 */
void cfUpdate(peakData *peak)
{
  /* Update peak (integer) location with hysteresis. */
  if(fabs(peak->params[XCENTER] - (double)peak->xi) > HYSTERESIS){
    peak->xi = (int)round(peak->params[XCENTER]);
  }
  if(fabs(peak->params[YCENTER] - (double)peak->yi) > HYSTERESIS){
    peak->yi = (int)round(peak->params[YCENTER]);
  }
}


/*
 * cfUpdate2D()
 *
 * Update for 2D spline fitting.
 *
 * fit_data - pointer to a fitData structure.
 * delta - the deltas for different parameters.
 */
void cfUpdate2D(fitData *fit_data, double *delta)
{
  peakData *peak;

  peak = fit_data->working_peak;

  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], YCENTER);
  mFitUpdateParam(peak, delta[3], BACKGROUND);

  cfUpdate(peak);
}


/*
 * cfUpdate3D()
 *
 * Update for a 3D spline fitting.
 *
 * fit_data - pointer to a fitData structure.
 * delta - the deltas for different parameters.
 */
void cfUpdate3D(fitData *fit_data, double *delta)
{
  peakData *peak;

  peak = fit_data->working_peak;

  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], YCENTER);
  mFitUpdateParam(peak, delta[3], ZCENTER);
  mFitUpdateParam(peak, delta[4], BACKGROUND);

  cfUpdate(peak);

  /* Keep Z in a fixed range. */
  cfZRangeCheck(fit_data);
}

/*
 * cfZRangeCheck()
 *
 * Keep peak z value inside a specific range and also update zi.
 * 
 * This is a separate function as multiplane also uses it.
 */
void cfZRangeCheck(fitData *fit_data)
{
  double maxz;
  peakData *peak;
  splineFit *spline_fit;

  peak = fit_data->working_peak;
  spline_fit = (splineFit *)fit_data->fit_model;  

  /* Force z value to stay in range. */
  if(peak->params[ZCENTER] < 1.0e-12){
    peak->params[ZCENTER] = 1.0e-12;
  }
  maxz = ((double)spline_fit->spline_size_z) - 1.0e-12;
  if(peak->params[ZCENTER] > maxz){
    peak->params[ZCENTER] = maxz;
  }
  ((splinePeak *)peak->peak_model)->zi = (int)(peak->params[ZCENTER]);
}
