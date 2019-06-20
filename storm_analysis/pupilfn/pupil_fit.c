/*
 * Fit multiple, possible overlapping, pupil functions to 
 * image data.
 *
 * Hazen 10/17
 */

/* Include */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pupil_fit.h"


/*
 * pfitAllocPeaks()
 *
 * Allocate storage for pfitPeaks. Note that this does not allocate
 * space for the psf_c and psf_r elements, which is done in pfitNewPeaks().
 */
void pfitAllocPeaks(peakData *new_peaks, int n_peaks)
{
  int i;

  for(i=0;i<n_peaks;i++){
    new_peaks[i].peak_model = (pupilPeak *)malloc(sizeof(pupilPeak));
    ((pupilPeak *)new_peaks[i].peak_model)->psf_c = NULL;
    ((pupilPeak *)new_peaks[i].peak_model)->psf_r = NULL;
  }
}


/* 
 * pfitCalcJH3D()
 *
 * Calculate Jacobian and Hessian for the pupil function.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void pfitCalcJH3D(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m,n;
  double height,fi,rqei,t1,t2,xi;
  double jt[5];
  double *dx_c,*dx_r,*dy_c,*dy_r,*dz_c,*dz_r,*psf_c,*psf_r;
  peakData *peak;
  pupilPeak *pupil_peak;
  pupilFit *pupil_fit;

  /* Initializations. */
  peak = fit_data->working_peak;
  pupil_peak = (pupilPeak *)peak->peak_model;
  pupil_fit = (pupilFit *)fit_data->fit_model;
  
  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }

  /* 
   * Calculate derivatives in x,y,z. 
   */
  
  /* Translate PF by dx, dy, dz. */
  pfnTranslate(pupil_fit->pupil_data, pupil_peak->dx, pupil_peak->dy, pupil_peak->dz);

  pfnGetPSFdx(pupil_fit->pupil_data, pupil_fit->dx_r, pupil_fit->dx_c);
  pfnGetPSFdy(pupil_fit->pupil_data, pupil_fit->dy_r, pupil_fit->dy_c);
  pfnGetPSFdz(pupil_fit->pupil_data, pupil_fit->dz_r, pupil_fit->dz_c);
  
  /* 
   * Calculate jacobian and hessian. 
   */
  
  /* There are just for convenience. */
  psf_r = pupil_peak->psf_r;
  psf_c = pupil_peak->psf_c;
  dx_r = pupil_fit->dx_r;
  dx_c = pupil_fit->dx_c;
  dy_r = pupil_fit->dy_r;
  dy_c = pupil_fit->dy_c;
  dz_r = pupil_fit->dz_r;
  dz_c = pupil_fit->dz_c;  

  height = peak->params[HEIGHT];
  i = peak->yi * fit_data->image_size_x + peak->xi;
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    l = fit_data->roi_x_index[j]*fit_data->fit_size_y + fit_data->roi_y_index[j];

    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];

    /* Calculate derivatives. */
    jt[0] = rqei*(psf_r[l]*psf_r[l]+psf_c[l]*psf_c[l]);
    jt[1] = rqei*2.0*height*(psf_r[l]*dx_r[l]+psf_c[l]*dx_c[l]);
    jt[2] = rqei*2.0*height*(psf_r[l]*dy_r[l]+psf_c[l]*dy_c[l]);
    jt[3] = rqei*2.0*height*(psf_r[l]*dz_r[l]+psf_c[l]*dz_c[l]);
    jt[4] = rqei;

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


/*
 * pfitCalcPeakShape()
 *
 * Calculate peak shape.
 *
 * fit_data - pointer to a fitData structure.
 */
void pfitCalcPeakShape(fitData *fit_data)
{
  int i,j;
  double *psf_r,*psf_c;
  peakData *peak;
  pupilPeak *pupil_peak;
  pupilFit *pupil_fit;

  peak = fit_data->working_peak;
  pupil_peak = (pupilPeak *)peak->peak_model;
  pupil_fit = (pupilFit *)fit_data->fit_model;
  
  /* 
   * Calculate PSF shape using the pupil_function library.
   */  
  pupil_peak->dx = peak->params[XCENTER] - (double)peak->xi;
  pupil_peak->dy = peak->params[YCENTER] - (double)peak->yi; 
  pupil_peak->dz = peak->params[ZCENTER];

  /* Translate PF by dx, dy, dz. */
  pfnTranslate(pupil_fit->pupil_data, pupil_peak->dx, pupil_peak->dy, pupil_peak->dz);

  /* Get PSF (complex) values, save with the peak. */
  pfnGetPSF(pupil_fit->pupil_data, pupil_peak->psf_r, pupil_peak->psf_c);

  /* 
   * Calculate PSF. 
   * 
   * FIXME: This would be more straightforward if the pupil function library and this
   *        library agreed on which axis was x and which was y..
   */
  psf_r = pupil_peak->psf_r;
  psf_c = pupil_peak->psf_c;
  for(i=0;i<fit_data->roi_n_index;i++){
    j = fit_data->roi_x_index[i]*fit_data->fit_size_y + fit_data->roi_y_index[i];
    peak->psf[i] = psf_r[j]*psf_r[j]+psf_c[j]*psf_c[j];
  }
}


/*
 * pfitCleanup()
 *
 * Frees the fitData structure.
 *
 * fit_data - pointer to a fitData structure.
 */
void pfitCleanup(fitData *fit_data)
{
  int i;
  pupilPeak *pupil_peak;
  pupilFit *pupil_fit;

  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->max_nfit;i++){
      pupil_peak = (pupilPeak *)(fit_data->fit[i].peak_model);
      free(pupil_peak->psf_r);
      free(pupil_peak->psf_c);
      free(pupil_peak);
    }
  }
  
  pupil_peak = (pupilPeak *)(fit_data->working_peak->peak_model);
  free(pupil_peak->psf_r);
  free(pupil_peak->psf_c);
  free(pupil_peak);
    
  pupil_fit = (pupilFit *)fit_data->fit_model;
  free(pupil_fit->dx_r);
  free(pupil_fit->dx_c);
  free(pupil_fit->dy_r);
  free(pupil_fit->dy_c);
  free(pupil_fit->dz_r);
  free(pupil_fit->dz_c);
  pfnCleanup(pupil_fit->pupil_data);
  free(pupil_fit);

  mFitCleanup(fit_data);
}


/*
 * pfitCopyPeak()
 *
 * Copies the contents of peak structure into another peak structure.
 *
 * fit_data - pointer to a fitData structure.
 * original - pointer to a peakData structure.
 * copy - pointer to a peakData structure.
 */
void pfitCopyPeak(fitData *fit_data, peakData *original, peakData *copy)
{
  int n;
  pupilPeak *pupil_copy, *pupil_original;

  pupil_copy = (pupilPeak *)copy->peak_model;
  pupil_original = (pupilPeak *)original->peak_model;

  /* Allocate storage, if necessary. */
  if(copy->psf == NULL){
    n = fit_data->fit_size_x*fit_data->fit_size_y;
    copy->psf = (double *)malloc(sizeof(double)*n);
    pupil_copy->psf_r = (double *)malloc(sizeof(double)*n);
    pupil_copy->psf_c = (double *)malloc(sizeof(double)*n);
  }
  
  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(fit_data, original, copy);

  /* Copy the parts that are specific to Pupilfn. */
  pupil_copy->dx = pupil_original->dx;
  pupil_copy->dy = pupil_original->dy;
  pupil_copy->dz = pupil_original->dz;

  memcpy(pupil_copy->psf_r, pupil_original->psf_r, sizeof(double)*fit_data->roi_n_index);
  memcpy(pupil_copy->psf_c, pupil_original->psf_c, sizeof(double)*fit_data->roi_n_index);
}


/*
 * pfitFreePeaks()
 *
 * Frees a peakData array.
 *
 * peaks - Pointer to an array of peakData.
 * n_peaks - The size of the array.
 */
void pfitFreePeaks(peakData *peaks, int n_peaks)
{
  int i;
  pupilPeak *pupil_peak;

  for(i=0;i<n_peaks;i++){
    pupil_peak = (pupilPeak *)(peaks[i].peak_model);
    free(pupil_peak->psf_c);
    free(pupil_peak->psf_r);
    free(pupil_peak);
  }
}


/*
 * pfitInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * pupil_data - Pointer to a pupilData structure.
 * rqe - Pixel relative quantum efficiency.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* pfitInitialize(pupilData *pupil_data, double *rqe, double *scmos_calibration, double tol, int im_size_x, int im_size_y)
{
  int pupil_size;
  fitData* fit_data;
  pupilFit *pupil_fit;

  fit_data = mFitInitialize(rqe, scmos_calibration, tol, im_size_x, im_size_y);
  fit_data->jac_size = 5;

  pupil_size = pfnGetSize(pupil_data);
  fit_data->xoff = 0.5*((double)pupil_size);
  fit_data->yoff = 0.5*((double)pupil_size);

  mFitInitializeROIIndexing(fit_data, pupil_size);

  /*
   * Initialize fit model.
   */
  fit_data->fit_model = (pupilFit *)malloc(sizeof(pupilFit));
  pupil_fit = (pupilFit *)fit_data->fit_model;
  pupil_fit->pupil_data = pupil_data;
  pupil_fit->pupil_size = pupil_size;
  pupil_fit->dx_r = (double *)malloc(pupil_size*pupil_size*sizeof(double));
  pupil_fit->dx_c = (double *)malloc(pupil_size*pupil_size*sizeof(double));
  pupil_fit->dy_r = (double *)malloc(pupil_size*pupil_size*sizeof(double));
  pupil_fit->dy_c = (double *)malloc(pupil_size*pupil_size*sizeof(double));
  pupil_fit->dz_r = (double *)malloc(pupil_size*pupil_size*sizeof(double));
  pupil_fit->dz_c = (double *)malloc(pupil_size*pupil_size*sizeof(double));

  /* Allocate storage for the working peak. */
  fit_data->working_peak->peak_model = (pupilPeak *)malloc(sizeof(pupilPeak));
  fit_data->working_peak->psf = (double *)malloc(sizeof(double)*pupil_size*pupil_size);
  ((pupilPeak *)fit_data->working_peak->peak_model)->psf_r = (double *)malloc(sizeof(double)*pupil_size*pupil_size);
  ((pupilPeak *)fit_data->working_peak->peak_model)->psf_c = (double *)malloc(sizeof(double)*pupil_size*pupil_size);

  /* Set function pointers. */
  fit_data->fn_alloc_peaks = &pfitAllocPeaks;
  fit_data->fn_calc_JH = &pfitCalcJH3D;
  fit_data->fn_calc_peak_shape = &pfitCalcPeakShape;
  fit_data->fn_check = &mFitCheck;
  fit_data->fn_copy_peak = &pfitCopyPeak;
  fit_data->fn_update = &pfitUpdate3D;
  
  return fit_data;
}


/*
 * pfitNewPeaks
 *
 * fit_data - Pointer to a fitData structure.
 * peak_params - Input values for the peak parameters.
 * n_peaks - The number of peaks.
 */
void pfitNewPeaks(fitData *fit_data, double *peak_params, char *p_type, int n_peaks)
{
  int i,j,n;
  int start,stop;
  peakData *peak;
  pupilPeak *pupil_peak;

  if(VERBOSE){
    printf("cfNP %d\n", n_peaks);
  }

  /*
   * Generic initializations.
   */
  mFitNewPeaks(fit_data, n_peaks);

  /* Pupil function specific initializations. */
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
      pupil_peak = (pupilPeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;

      /* Allocate space for saving the PSF. */
      if(peak->psf == NULL){
	n = fit_data->fit_size_x*fit_data->fit_size_y;
	peak->psf = (double *)malloc(sizeof(double)*n);
	pupil_peak->psf_r = (double *)malloc(sizeof(double)*n);
	pupil_peak->psf_c = (double *)malloc(sizeof(double)*n);
      }

      /* Calculate (integer) peak locations. */
      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);

      /* Arbitrary initial values for BACKGROUND, HEIGHT. */
      peak->params[BACKGROUND] = 1.0;
      peak->params[HEIGHT] = 1.0;

      /* Copy into working peak. */
      pfitCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	pfitCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }

      /* Calculate peak shape (of working peak). */
      pfitCalcPeakShape(fit_data);

      /* Estimate best starting background. */
      mFitEstimatePeakBackground(fit_data);
      
      if(!strcmp(p_type, "finder")){

	/* Estimate best starting height. */
	mFitEstimatePeakHeight(fit_data);
	
	if(fit_data->working_peak->params[HEIGHT] < fit_data->minimum_height){
	  fit_data->working_peak->params[HEIGHT] = fit_data->minimum_height;
	}
	
	/* Check that the initial height is positive, error it out if not. */
	if(fit_data->working_peak->params[HEIGHT] <= 0.0){
	  printf("Warning peak %d has negative estimated height!\n", (i-start));
	  fit_data->working_peak->status = ERROR;
	  pfitCopyPeak(fit_data, fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      pfitCopyPeak(fit_data, fit_data->working_peak, peak);      
    }
  }
  /*
   * "pre-specified" parameters, these are the peak x, y, z, background
   * and height as an n_peaks x 5 array.
   */
  else if(!strcmp(p_type, "text") || !strcmp(p_type, "hdf5")){
    for(i=start;i<stop;i++){
      j = 5*(i-start);
      peak = &fit_data->fit[i];
      pupil_peak = (pupilPeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      peak->params[BACKGROUND] = peak_params[j+3];
      peak->params[HEIGHT] = peak_params[j+4];

      /* Allocate space for saving the PSF. */
      if(peak->psf == NULL){
	n = fit_data->fit_size_x*fit_data->fit_size_y;
	peak->psf = (double *)malloc(sizeof(double)*n);
	pupil_peak->psf_r = (double *)malloc(sizeof(double)*n);
	pupil_peak->psf_c = (double *)malloc(sizeof(double)*n);
      }
      
      /* Calculate (integer) peak locations. */
      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);

      /* Copy into working peak. */
      pfitCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	pfitCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      pfitCalcPeakShape(fit_data);
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      pfitCopyPeak(fit_data, fit_data->working_peak, peak);
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
    pfitCopyPeak(fit_data, peak, fit_data->working_peak);
    fit_data->fn_error_fn(fit_data);
    pfitCopyPeak(fit_data, fit_data->working_peak, peak);
  }

  fit_data->nfit = stop;
}
  

/*
 * pfitSetZRange()
 *
 * Set the fitting range for Z (in microns).
 */
void pfitSetZRange(fitData *fit_data, double min_z, double max_z)
{
  pupilFit *pupil_fit;

  pupil_fit = (pupilFit *)fit_data->fit_model;
  pupil_fit->min_z = min_z;
  pupil_fit->max_z = max_z;
}


/*
 * pfitUpdate3D()
 *
 * Update for a pupil function fitting.
 *
 * fit_data - pointer to a fitData structure.
 * delta - the deltas for different parameters.
 */
void pfitUpdate3D(fitData *fit_data, double *delta)
{
  peakData *peak;

  peak = fit_data->working_peak;

  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], YCENTER);
  mFitUpdateParam(peak, delta[3], ZCENTER);
  mFitUpdateParam(peak, delta[4], BACKGROUND);

  /* Update peak (integer) location with hysteresis. */
  mFitUpdate(peak);

  /* Keep Z in a fixed range. */
  pfitZRangeCheck(fit_data);
}

/*
 * pfitZRangeCheck()
 *
 * Keep peak z value inside a specific range. 
 * 
 * This is a separate function as multiplane also uses it.
 */
void pfitZRangeCheck(fitData *fit_data)
{
  peakData *peak;
  pupilFit *pupil_fit;

  peak = fit_data->working_peak;
  pupil_fit = (pupilFit *)fit_data->fit_model;

  if(peak->params[ZCENTER] < pupil_fit->min_z){
    peak->params[ZCENTER] = pupil_fit->min_z;
  }

  if(peak->params[ZCENTER] > pupil_fit->max_z){
    peak->params[ZCENTER] = pupil_fit->max_z;
  }  
}
