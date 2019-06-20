/*
 * Fit multiple, possible overlapping, FFT PSFs to 
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

#include "fft_fit.h"


/*
 * ftFitAllocPeaks()
 *
 * Allocate storage for psfFFTPeaks. Note that this does not allocate
 * space for the psf element, which is done in ftFitNewPeaks().
 */
void ftFitAllocPeaks(peakData *new_peaks, int n_peaks)
{
  int i;

  for(i=0;i<n_peaks;i++){
    new_peaks[i].peak_model = (psfFFTPeak *)malloc(sizeof(psfFFTPeak));
  }
}


/* 
 * ftFitCalcJH3D()
 *
 * Calculate Jacobian and Hessian for the psf FFT function.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void ftFitCalcJH3D(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m,n;
  double height,fi,rqei,t1,t2,xi;
  double jt[5];
  double *dx,*dy,*dz,*psf;
  peakData *peak;
  psfFFTPeak *psf_fft_peak;
  psfFFTFit *psf_fft_fit;

  /* Initializations. */
  peak = fit_data->working_peak;
  psf_fft_peak = (psfFFTPeak *)peak->peak_model;
  psf_fft_fit = (psfFFTFit *)fit_data->fit_model;
  
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
  pFTTranslate(psf_fft_fit->psf_fft_data, psf_fft_peak->dx, psf_fft_peak->dy, psf_fft_peak->dz);

  pFTGetPSFdx(psf_fft_fit->psf_fft_data, psf_fft_fit->dx);
  pFTGetPSFdy(psf_fft_fit->psf_fft_data, psf_fft_fit->dy);
  pFTGetPSFdz(psf_fft_fit->psf_fft_data, psf_fft_fit->dz);
  
  /* 
   * Calculate jacobian and hessian. 
   */
  
  /* There are just for convenience. */
  psf = peak->psf;
  dx = psf_fft_fit->dx;
  dy = psf_fft_fit->dy;
  dz = psf_fft_fit->dz;

  height = peak->params[HEIGHT];
  i = peak->yi * fit_data->image_size_x + peak->xi;
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    l = fit_data->roi_y_index[j]*fit_data->fit_size_x + fit_data->roi_x_index[j];
    
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];

    /* Calculate derivatives. */
    jt[0] = rqei*psf[j];
    jt[1] = rqei*height*(dx[l]);
    jt[2] = rqei*height*(dy[l]);
    jt[3] = rqei*height*(dz[l]);
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
 * ftFitCalcPeakShape()
 *
 * Calculate peak shape.
 *
 * fit_data - pointer to a fitData structure.
 */
void ftFitCalcPeakShape(fitData *fit_data)
{
  int i,j;
  
  peakData *peak;
  psfFFTPeak *psf_fft_peak;
  psfFFTFit *psf_fft_fit;

  peak = fit_data->working_peak;
  psf_fft_peak = (psfFFTPeak *)peak->peak_model;
  psf_fft_fit = (psfFFTFit *)fit_data->fit_model;
  
  /* 
   * Calculate PSF shape using the psf_fft library.
   */  
  psf_fft_peak->dx = peak->params[XCENTER] - (double)peak->xi;
  psf_fft_peak->dy = peak->params[YCENTER] - (double)peak->yi; 
  psf_fft_peak->dz = peak->params[ZCENTER];

  /* Translate PSF by dx, dy, dz. */
  pFTTranslate(psf_fft_fit->psf_fft_data, psf_fft_peak->dx, psf_fft_peak->dy, psf_fft_peak->dz);

  /* Get PSF values, save with the peak. */
  pFTGetPSF(psf_fft_fit->psf_fft_data, psf_fft_fit->psf);
  for(i=0;i<fit_data->roi_n_index;i++){
    j = fit_data->roi_y_index[i]*fit_data->fit_size_x + fit_data->roi_x_index[i];
    peak->psf[i] = psf_fft_fit->psf[j];
  }
}


/*
 * ftFitCleanup()
 *
 * Frees the fitData structure.
 *
 * fit_data - pointer to a fitData structure.
 */
void ftFitCleanup(fitData *fit_data)
{
  int i;
  psfFFTPeak *psf_fft_peak;
  psfFFTFit *psf_fft_fit;

  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->max_nfit;i++){
      psf_fft_peak = fit_data->fit[i].peak_model;
      free(psf_fft_peak);
    }
  }
  
  free(fit_data->working_peak->peak_model);
    
  psf_fft_fit = (psfFFTFit *)fit_data->fit_model;
  free(psf_fft_fit->dx);
  free(psf_fft_fit->dy);
  free(psf_fft_fit->dz);
  free(psf_fft_fit->psf);
  pFTCleanup(psf_fft_fit->psf_fft_data);
  free(psf_fft_fit);

  mFitCleanup(fit_data);
}


/*
 * ftFitCopyPeak()
 *
 * Copies the contents of peak structure into another peak structure.
 *
 * original - pointer to a peakData structure.
 * copy - pointer to a peakData structure.
 */
void ftFitCopyPeak(fitData *fit_data, peakData *original, peakData *copy)
{
  psfFFTPeak *psf_fft_copy, *psf_fft_original;

  psf_fft_copy = (psfFFTPeak *)copy->peak_model;
  psf_fft_original = (psfFFTPeak *)original->peak_model;

  /* Allocate storage, if necessary. */
  if(copy->psf == NULL){
    copy->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
  }
  
  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(fit_data, original, copy);

  /* Copy the parts that are specific to psf FFT. */
  psf_fft_copy->dx = psf_fft_original->dx;
  psf_fft_copy->dy = psf_fft_original->dy;
  psf_fft_copy->dz = psf_fft_original->dz;
}


/*
 * ftFitFreePeaks()
 *
 * Frees a peakData array.
 *
 * peaks - Pointer to an array of peakData.
 * n_peaks - The size of the array.
 */
void ftFitFreePeaks(peakData *peaks, int n_peaks)
{
  int i;

  for(i=0;i<n_peaks;i++){
    free(peaks[i].peak_model);
  }
}


/*
 * ftFitInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * psf_fft_data - Pointer to a psfFFT structure.
 * rqe - Pixel relative quantum efficiency.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* ftFitInitialize(psfFFT *psf_fft_data, double *rqe, double *scmos_calibration, double tol, int im_size_x, int im_size_y)
{
  int xy_size;
  fitData* fit_data;
  psfFFTFit *psf_fft_fit;

  fit_data = mFitInitialize(rqe, scmos_calibration, tol, im_size_x, im_size_y);
  fit_data->jac_size = 5;

  /*
   * Initialize fit model.
   */
  fit_data->fit_model = (psfFFTFit *)malloc(sizeof(psfFFTFit));
  psf_fft_fit = (psfFFTFit *)fit_data->fit_model;
  psf_fft_fit->psf_fft_data = psf_fft_data;
  psf_fft_fit->psf_x = pFTGetXSize(psf_fft_data);
  psf_fft_fit->psf_y = pFTGetYSize(psf_fft_data);
  psf_fft_fit->psf_z = pFTGetZSize(psf_fft_data);

  /*
   * FIXME: Only tested with even size PSFs. The PSF measurement tools
   * all return even size PSFs.
   */
  fit_data->xoff = 0.5*((double)psf_fft_fit->psf_x) - 1.0;
  fit_data->yoff = 0.5*((double)psf_fft_fit->psf_y) - 1.0;

  /*
   * Enforce that the PSF is square in X/Y.
   */
  if(psf_fft_fit->psf_x != psf_fft_fit->psf_y){
    printf("PSF must be square in X/Y!");
    exit(EXIT_FAILURE);
  }
  mFitInitializeROIIndexing(fit_data, psf_fft_fit->psf_x);
  
  xy_size = psf_fft_fit->psf_x * psf_fft_fit->psf_y;
  
  psf_fft_fit->dx = (double *)malloc(sizeof(double)*xy_size);
  psf_fft_fit->dy = (double *)malloc(sizeof(double)*xy_size);
  psf_fft_fit->dz = (double *)malloc(sizeof(double)*xy_size);
  psf_fft_fit->psf = (double *)malloc(sizeof(double)*xy_size);

  /* Allocate storage for the working peak. */
  fit_data->working_peak->peak_model = (psfFFTPeak *)malloc(sizeof(psfFFTPeak));

  /* Set function pointers. */
  fit_data->fn_alloc_peaks = &ftFitAllocPeaks;
  fit_data->fn_calc_JH = &ftFitCalcJH3D;
  fit_data->fn_calc_peak_shape = &ftFitCalcPeakShape;
  fit_data->fn_check = &mFitCheck;
  fit_data->fn_copy_peak = &ftFitCopyPeak;
  fit_data->fn_update = &ftFitUpdate3D;
  
  return fit_data;
}


/*
 * ftFitNewPeaks
 *
 * fit_data - Pointer to a fitData structure.
 * peak_params - Input values for the peak parameters.
 * n_peaks - The number of peaks.
 */
void ftFitNewPeaks(fitData *fit_data, double *peak_params, char *p_type, int n_peaks)
{
  int i,j;
  int start,stop;
  peakData *peak;

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

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      
      /* Allocate space for saving the PSF. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
      }

      /* Calculate (integer) peak locations. */
      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);

      /* Arbitrary initial values for BACKGROUND, HEIGHT. */
      peak->params[BACKGROUND] = 1.0;
      peak->params[HEIGHT] = 1.0;

      /* Copy into working peak. */
      ftFitCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	ftFitCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }

      /* Calculate peak shape (of working peak). */
      ftFitCalcPeakShape(fit_data);

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
	  printf("Warning peak %d has negative estimated height! %.2f\n", (i-start), fit_data->working_peak->params[HEIGHT]);
	  fit_data->working_peak->status = ERROR;
	  ftFitCopyPeak(fit_data, fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      ftFitCopyPeak(fit_data, fit_data->working_peak, peak);      
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

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      peak->params[BACKGROUND] = peak_params[j+3];
      peak->params[HEIGHT] = peak_params[j+4];

      /* Allocate space for saving the PSF. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
      }
      
      /* Calculate (integer) peak locations. */
      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);

      /* Copy into working peak. */
      ftFitCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	ftFitCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      ftFitCalcPeakShape(fit_data);
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      ftFitCopyPeak(fit_data, fit_data->working_peak, peak);
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
    ftFitCopyPeak(fit_data, peak, fit_data->working_peak);
    fit_data->fn_error_fn(fit_data);
    ftFitCopyPeak(fit_data, fit_data->working_peak, peak);
  }

  fit_data->nfit = stop;
}


/*
 * ftFitUpdate3D()
 *
 * Update for a PSF FFT fitting.
 *
 * fit_data - pointer to a fitData structure.
 * delta - the deltas for different parameters.
 */
void ftFitUpdate3D(fitData *fit_data, double *delta)
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

  ftFitZRangeCheck(fit_data);
}

/*
 * ftFitZRangeCheck()
 *
 * Keep Z in a fixed range. The FFT has periodic boundary conditions
 * so we need to stay in the range z_size / 2 (the PSF is the middle
 * slice of the FFT).
 */
void ftFitZRangeCheck(fitData *fit_data)
{
  double max_z;
  peakData *peak;
  psfFFTFit *psf_fft_fit;

  peak = fit_data->working_peak;
  psf_fft_fit = (psfFFTFit *)fit_data->fit_model;
  max_z = ((double)psf_fft_fit->psf_z)*0.5 - 1.0e-12;
  
  if(peak->params[ZCENTER] < -max_z){
    peak->params[ZCENTER] = -max_z;
  }
  if(peak->params[ZCENTER] > max_z){
    peak->params[ZCENTER] = max_z;
  }  
}
