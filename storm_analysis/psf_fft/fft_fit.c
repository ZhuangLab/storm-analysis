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

  peak->added++;
  
  if(TESTING){
    if(peak->added != 1){
      printf("Peak count error detected in ftFitAddPeak()! %d\n", peak->added);
      exit(EXIT_FAILURE);
    }
  }
  
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
      n = j * peak->size_x + k;
      fit_data->f_data[m] += height*psf[n];
      fit_data->bg_counts[m] += 1;
      fit_data->bg_data[m] += bg + fit_data->scmos_term[m];
    }
  }
}


/*
 * ftFitAllocPeaks()
 *
 * Allocate storage for psfFFTPeaks. Note that this does not allocate
 * space for the psf element, which is done in ftFitNewPeaks().
 */
struct peakData *ftFitAllocPeaks(int n_peaks)
{
  int i;
  peakData *new_peaks;

  new_peaks = (peakData *)malloc(sizeof(peakData)*n_peaks);  
  for(i=0;i<n_peaks;i++){
    new_peaks[i].peak_model = (psfFFTPeak *)malloc(sizeof(psfFFTPeak));
    ((psfFFTPeak *)new_peaks[i].peak_model)->psf = NULL;
  }
  return new_peaks;
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
      o = j * peak->size_x + k;
      
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
 * ftFitCalcPeakShape()
 *
 * Calculate peak shape.
 *
 * fit_data - pointer to a fitData structure.
 */
void ftFitCalcPeakShape(fitData *fit_data)
{
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
  psfFFTPeak *psf_fft_peak;

  for(i=0;i<n_peaks;i++){
    psf_fft_peak = (psfFFTPeak *)(peaks[i].peak_model);
    free(psf_fft_peak->psf);
    free(psf_fft_peak);
  }
  free(peaks);
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
  fit_data->fn_alloc_peaks = &ftFitAllocPeaks;
  fit_data->fn_calc_JH = &ftFitCalcJH3D;
  fit_data->fn_calc_peak_shape = &ftFitCalcPeakShape;
  fit_data->fn_check = &mFitCheck;
  fit_data->fn_copy_peak = &ftFitCopyPeak;
  fit_data->fn_peak_sum = &ftFitPeakSum;
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
void ftFitNewPeaks(fitData *fit_data, double *peak_params, char *p_type, int n_peaks)
{
  int i,j,k,l,m,n,xc,yc;
  int start,stop;
  double sp,sx,t1;
  peakData *peak;
  psfFFTPeak *psf_fft_peak;

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
      psf_fft_peak = (psfFFTPeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      
      /*
       * Note: Even though these are the same for every peak (as the PSF FFT
       *       does not change size during fitting), they are duplicated for 
       *       each peak for the benefit of the mFit functions.
       */
      peak->size_x = ((psfFFTFit *)fit_data->fit_model)->psf_x;
      peak->size_y = ((psfFFTFit *)fit_data->fit_model)->psf_y;

      /* Allocate space for saving the PSF. */
      if(psf_fft_peak->psf == NULL){
	psf_fft_peak->psf = (double *)malloc(sizeof(double)*peak->size_x*peak->size_y);
      }
      
      /* Calculate (integer) peak locations. */
      peak->xi = (int)round(peak->params[XCENTER]);
      peak->yi = (int)round(peak->params[YCENTER]);

      /* Estimate background. */
      xc = (int)round(peak_params[j]);
      yc = (int)round(peak_params[j+1]);      
      peak->params[BACKGROUND] = fit_data->bg_estimate[yc * fit_data->image_size_x + xc];      

      /* Arbitrary initial value for HEIGHT. */
      peak->params[HEIGHT] = 1.0;      

      /* Copy into working peak. */
      ftFitCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	ftFitCopyPeak(fit_data->working_peak, peak);
	continue;
      }

      if(!strcmp(p_type, "finder")){

	/* Calculate peak shape (of working peak). */
	ftFitCalcPeakShape(fit_data);

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
	psf_fft_peak = (psfFFTPeak *)fit_data->working_peak->peak_model;
	k = peak->yi * fit_data->image_size_x + peak->xi;
	sp = 0.0;  /* This is fi*fi*fi. */
	sx = 0.0;  /* This is fi*fi*xi. */
	for(l=0;l<peak->size_y;l++){
	  for(m=0;m<peak->size_x;m++){
	    n = l * fit_data->image_size_x + m + k;
	    t1 = psf_fft_peak->psf[l*peak->size_x+m];
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
	  printf("Warning peak %d has negative estimated height! %.2f\n", (i-start), fit_data->working_peak->params[HEIGHT]);
	  fit_data->working_peak->status = ERROR;
	  ftFitCopyPeak(fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      ftFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      ftFitCopyPeak(fit_data->working_peak, peak);      
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
      psf_fft_peak = (psfFFTPeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      peak->params[BACKGROUND] = peak_params[j+3];
      peak->params[HEIGHT] = peak_params[j+4];

      /*
       * Note: Even though these are the same for every peak (as the pupil
       *       function does not change size during fitting), they are 
       *       duplicated for each peak for the benefit of the mFit functions.
       */
      peak->size_x = ((psfFFTFit *)fit_data->fit_model)->psf_x;
      peak->size_y = ((psfFFTFit *)fit_data->fit_model)->psf_y;

      /* Allocate space for saving the PSF. */
      if(psf_fft_peak->psf == NULL){
	psf_fft_peak->psf = (double *)malloc(sizeof(double)*peak->size_x*peak->size_y);
      }
      
      /* Calculate (integer) peak locations. */
      peak->xi = (int)round(peak->params[XCENTER]);
      peak->yi = (int)round(peak->params[YCENTER]);

      /* Copy into working peak. */
      ftFitCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	ftFitCopyPeak(fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      ftFitCalcPeakShape(fit_data);
      ftFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      ftFitCopyPeak(fit_data->working_peak, peak);
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
    ftFitCopyPeak(peak, fit_data->working_peak);
    mFitCalcErr(fit_data);
    ftFitCopyPeak(fit_data->working_peak, peak);
  }

  fit_data->nfit = stop;

  /* Reset the clamp values on all the peaks. */
  if(USECLAMP){
    mFitResetClampValues(fit_data);
  }  
}


/*
 * ftFitPeakSum()
 *
 * Return the integral of the PSF.
 */
double ftFitPeakSum(peakData *peak)
{
  int i;
  double sum;
  psfFFTPeak *psf_fft_peak;
    
  psf_fft_peak = (psfFFTPeak *)peak->peak_model;

  sum = 0.0;
  for(i=0;i<(peak->size_x * peak->size_y);i++){
    sum += psf_fft_peak->psf[i];
  }
  sum = sum*peak->params[HEIGHT];

  return sum;
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

  peak->added--;

  if(TESTING){
    if(peak->added != 0){
      printf("Peak count error detected in ftFitSubtractPeak()! %d\n", peak->added);
      exit(EXIT_FAILURE);
    }
  }
  
  psf = psf_fft_peak->psf;
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  height = peak->params[HEIGHT];
  for (j=0;j<peak->size_y;j++){
    for (k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      n = j*peak->size_x + k;
      fit_data->f_data[m] -= height*psf[n];
      fit_data->bg_counts[m] -= 1;
      fit_data->bg_data[m] -= (bg + fit_data->scmos_term[m]);
    }
  }
}


/*
 * ftFitUpdate()
 *
 * Updates peak location with hysteresis.
 *
 * peak - pointer to a peakData structure.
 */
void ftFitUpdate(peakData *peak)
{
  if(fabs(peak->params[XCENTER] - (double)peak->xi) > HYSTERESIS){
    peak->xi = (int)round(peak->params[XCENTER]);
  }
  if(fabs(peak->params[YCENTER] - (double)peak->yi) > HYSTERESIS){
    peak->yi = (int)round(peak->params[YCENTER]);
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
  ftFitUpdate(peak);

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
