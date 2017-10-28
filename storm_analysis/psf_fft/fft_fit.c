/*
 * Fit multiple, possible overlapping, FFT PSFs to 
 * image data.
 *
 * Hazen 10/17
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fft_fit.h"

/*
 * ftFitAddPeak()
 *
 * Calculate peak shape and add the working peak to the 
 * foreground and background data arrays.
 *
 * fit_data - pointer to a fitData structure.
 */
void ftFitAddPeak(fitData *fit_data)
{
  int j,k,l,m,n;
  double bg,height;
  double *psf;
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
  pFTGetPSF(psf_fft_fit->psf_fft_data, psf_fft_peak->psf);
  
  /* 
   * Add peak to the foreground and background arrays. 
   */
  psf = psf_fft_peak->psf;
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  height = peak->params[HEIGHT];
  for (j=0;j<peak->size_y;j++){
    for (k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;

      /*
       * With this indexing we are taking the transpose of the PSF. This
       * convention is also followed in calcJH3D() and subtractPeak().
       */
      n = j * peak->size_y + k;
      fit_data->f_data[m] += height*psf[n];
      fit_data->bg_counts[m] += 1;
      fit_data->bg_data[m] += bg + fit_data->scmos_term[m];
    }
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
  int i,j,k,l,m,n,o;
  double height,fi,t1,t2,xi;
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
  psf = psf_fft_peak->psf;
  dx = psf_fft_fit->dx;
  dy = psf_fft_fit->dy;
  dz = psf_fft_fit->dz;

  height = peak->params[HEIGHT];
  i = peak->yi * fit_data->image_size_x + peak->xi;
  for(j=0;j<peak->size_y;j++){
    for(k=0;k<peak->size_x;k++){
      l = i + j * fit_data->image_size_x + k;
      o = j * peak->size_y + k;
      
      fi = fit_data->f_data[l] + fit_data->bg_data[l] / ((double)fit_data->bg_counts[l]);
      xi = fit_data->x_data[l];

      /* Calculate derivatives. */
      jt[0] = psf[o];
      jt[1] = height*(dx[o]);
      jt[2] = height*(dy[o]);
      jt[3] = height*(dz[o]);
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
    for(i=0;i<fit_data->nfit;i++){
      psf_fft_peak = (psfFFTPeak *)(fit_data->fit[i].peak_model);
      free(psf_fft_peak->psf);
      free(psf_fft_peak);
    }
    free(fit_data->fit);
  }
  
  psf_fft_peak = (psfFFTPeak *)(fit_data->working_peak->peak_model);
  free(psf_fft_peak->psf);
  free(psf_fft_peak);
    
  psf_fft_fit = (psfFFTFit *)fit_data->fit_model;
  free(psf_fft_fit->dx);
  free(psf_fft_fit->dy);
  free(psf_fft_fit->dz);
  pFTCleanup(psf_fft_fit->psf_fft_data);

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
void ftFitCopyPeak(peakData *original, peakData *copy)
{
  int i;
  psfFFTPeak *psf_fft_copy, *psf_fft_original;

  psf_fft_copy = (psfFFTPeak *)copy->peak_model;
  psf_fft_original = (psfFFTPeak *)original->peak_model;

  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(original, copy);

  /* Copy the parts that are specific to psf FFT. */
  psf_fft_copy->dx = psf_fft_original->dx;
  psf_fft_copy->dy = psf_fft_original->dy;
  psf_fft_copy->dz = psf_fft_original->dz;

  for(i=0;i<(copy->size_x*copy->size_y);i++){
    psf_fft_copy->psf[i] = psf_fft_original->psf[i];
  }
}


/*
 * ftFitInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * psf_fft_data - Pointer to a psfFFT structure.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * clamp - The starting clamp values for each peak.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* ftFitInitialize(psfFFT *psf_fft_data, double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y)
{
  int xy_size;
  fitData* fit_data;
  psfFFTFit *psf_fft_fit;

  fit_data = mFitInitialize(scmos_calibration, clamp, tol, im_size_x, im_size_y);
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

  fit_data->xoff = 0.5*((double)psf_fft_fit->psf_x);
  fit_data->yoff = 0.5*((double)psf_fft_fit->psf_y);

  xy_size = psf_fft_fit->psf_x * psf_fft_fit->psf_y;
  
  psf_fft_fit->dx = (double *)malloc(sizeof(double)*xy_size);
  psf_fft_fit->dy = (double *)malloc(sizeof(double)*xy_size);
  psf_fft_fit->dz = (double *)malloc(sizeof(double)*xy_size);

  /* Allocate storage for the working peak. */
  fit_data->working_peak->peak_model = (psfFFTPeak *)malloc(sizeof(psfFFTPeak));
  ((psfFFTPeak *)fit_data->working_peak->peak_model)->psf = (double *)malloc(sizeof(double)*xy_size);

  /* Set function pointers. */
  fit_data->fn_add_peak = &ftFitAddPeak;
  fit_data->fn_calc_JH = &ftFitCalcJH3D;
  fit_data->fn_check = &mFitCheck;
  fit_data->fn_copy_peak = &ftFitCopyPeak;
  fit_data->fn_subtract_peak = &ftFitSubtractPeak;  
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
void ftFitNewPeaks(fitData *fit_data, double *peak_params, int n_peaks)
{
  int i;
  peakData *peak;
  psfFFTPeak *psf_fft_peak;

  /*
   * Free old peaks, if necessary.
   */
  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->nfit;i++){
      psf_fft_peak = ((psfFFTPeak *)fit_data->fit[i].peak_model);
      free(psf_fft_peak->psf);
    }
    free(fit_data->fit);
  }

  /*
   * Generic initializations.
   */
  mFitNewPeaks(fit_data, peak_params, n_peaks);
  
  /*
   * psf FFT specific initializations.
   */
  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    peak->peak_model = (psfFFTPeak *)malloc(sizeof(psfFFTPeak));
    psf_fft_peak = (psfFFTPeak *)peak->peak_model;

    /* Initial location. */
    peak->params[HEIGHT]     = peak_params[i*NPEAKPAR+HEIGHT];
    peak->params[XCENTER]    = peak_params[i*NPEAKPAR+XCENTER] - fit_data->xoff;
    peak->params[YCENTER]    = peak_params[i*NPEAKPAR+YCENTER] - fit_data->yoff;
    peak->params[BACKGROUND] = peak_params[i*NPEAKPAR+BACKGROUND];
    peak->params[ZCENTER]    = peak_params[i*NPEAKPAR+ZCENTER] - fit_data->zoff;

    /* 
     * These are not used, but need to be initialized so that they do not 
     * come out as confusing garbage. 
     */
    peak->params[XWIDTH] = 0.5;
    peak->params[YWIDTH] = 0.5;
    
    /*
     * Note: Even though these are the same for every peak (as the spline
     *       does not change size during fitting), they are duplicated
     *       for each peak for the benefit of the mFit functions.
     */
    peak->size_x = ((psfFFTFit *)fit_data->fit_model)->psf_x;
    peak->size_y = ((psfFFTFit *)fit_data->fit_model)->psf_y;

    /* Allocate space for saving the PSF. */
    psf_fft_peak->psf = (double *)malloc(sizeof(double)*peak->size_x*peak->size_y);

    /* Calculate (integer) peak locations. */
    peak->xi = (int)round(peak->params[XCENTER]);
    peak->yi = (int)round(peak->params[YCENTER]);

    /*
     * Add the peak to the fit. 
     */
    ftFitCopyPeak(peak, fit_data->working_peak);
    ftFitAddPeak(fit_data);
    ftFitCopyPeak(fit_data->working_peak, peak);    
  }

  /*
   * Initial error calculation. 
   */
  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    ftFitCopyPeak(peak, fit_data->working_peak);
    mFitCalcErr(fit_data);
    ftFitCopyPeak(fit_data->working_peak, peak);
  }
}


/*
 * ftFitSubtractPeak()
 *
 * Subtract the working peak out of the current fit, basically 
 * this just undoes addPeak().
 *
 * fit_data - pointer to a fitData structure.
 */
void ftFitSubtractPeak(fitData *fit_data)
{
  int j,k,l,m,n;
  double bg,height;
  double *psf;
  peakData *peak;
  psfFFTPeak *psf_fft_peak;

  peak = fit_data->working_peak;
  psf_fft_peak = (psfFFTPeak *)peak->peak_model;

  psf = psf_fft_peak->psf;
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  height = peak->params[HEIGHT];
  for (j=0;j<peak->size_y;j++){
    for (k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      n = j*peak->size_y + k;
      fit_data->f_data[m] -= height*psf[n];
      fit_data->bg_counts[m] -= 1;
      fit_data->bg_data[m] -= (bg + fit_data->scmos_term[m]);
    }
  }
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
  if(fabs(peak->params[XCENTER] - (double)peak->xi) > HYSTERESIS){
    peak->xi = (int)round(peak->params[XCENTER]);
  }
  if(fabs(peak->params[YCENTER] - (double)peak->yi) > HYSTERESIS){
    peak->yi = (int)round(peak->params[YCENTER]);
  }

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
