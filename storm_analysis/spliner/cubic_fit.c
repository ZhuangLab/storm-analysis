/*
 * Fit multiple, possible overlapping, cubic splines simultaneously to
 * image data.
 *
 * Hazen 09/18
 */

/* Include */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cubic_fit.h"


/*
 * cfAllocPeaks()
 *
 * Allocate storage for cfPeaks. Note that this does not allocate
 * space for the peak_values element, which is done in cfNewPeaks().
 *
 * FIXME: Not sure this is the best approach as both cfNewPeaks() and
 *        cfCopyPeaks() need to check whether peak_values is initialized.
 *        The advantage is that this function does not need to know
 *        the peak size.
 */
void cfAllocPeaks(peakData *new_peaks, int n_peaks)
{
  int i;

  for(i=0;i<n_peaks;i++){
    new_peaks[i].peak_model = (splinePeak *)malloc(sizeof(splinePeak));
  }
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
  int i,j,k,l,m,n,o, x_start, y_start;
  double height,fi,rqei,t1,t2,xi;
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
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    l = fit_data->roi_y_index[j];
    m = fit_data->roi_x_index[j];
      
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];

    jt[0] = rqei*peak->psf[j];
    jt[1] = -rqei*height*dxfAt2D(spline_fit->spline_data,l+y_start,m+x_start);
    jt[2] = -rqei*height*dyfAt2D(spline_fit->spline_data,l+y_start,m+x_start);
    jt[3] = rqei;
    
    /* Calculate jacobian. */
    t1 = 2.0*(1.0 - xi/fi);
    for(n=0;n<4;n++){
      jacobian[n] += t1*jt[n];
    }
    
    /* Calculate hessian. */
    t2 = 2.0*xi/(fi*fi);
    for(n=0;n<4;n++){
      for(o=n;o<4;o++){
	hessian[n*4+o] += t2*jt[n]*jt[o];
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
  int i,j,k,l,m,n,o,x_start,y_start,zi;
  double height,fi,rqei,t1,t2,xi;
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
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    l = fit_data->roi_y_index[j];
    m = fit_data->roi_x_index[j];

    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];

    jt[0] = rqei*peak->psf[j];
    jt[1] = -rqei*height*dxfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[2] = -rqei*height*dyfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[3] = rqei*height*dzfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[4] = rqei;
    
    /* Calculate jacobian. */
    t1 = 2.0*(1.0 - xi/fi);
    for(n=0;n<5;n++){
      jacobian[n] += t1*jt[n];
    }
	  
    /* Calculate hessian. */
    t2 = 2.0*xi/(fi*fi);
    for(n=0;n<5;n++){
      for(o=n;o<5;o++){
	hessian[n*5+o] += t2*jt[n]*jt[o];
      }
    }
  }
}


/* 
 * cfCalcJH3DALS()
 *
 * Calculate Jacobian and Hessian for the 3D spline (Anscombe least-squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void cfCalcJH3DALS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m,n,o,x_start,y_start,zi;
  double height,fi,t1,t2;
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
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    l = fit_data->roi_y_index[j];
    m = fit_data->roi_x_index[j];

    fi = fit_data->t_fi[k];

    /*
     * This is the LM algorithm according to wikipedia.
     *
     * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
     *
     * fi = 2.0 * (rqe_i * fit_fn_i(theta) + variance_i + 3.0/8.0) ^ 1/2
     * dfi/dtheta = (rqe_i * fit_fn_i(theta) + variance_i + 3.0/8.0) ^ -1/2 * rqe_i * dfit_fn_i(theta)/dtheta
     */
    t1 = fit_data->rqe[k]*2.0/fi;

    jt[0] = t1*peak->psf[j];
    jt[1] = -t1*height*dxfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[2] = -t1*height*dyfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[3] = t1*height*dzfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[4] = t1;

    /* Calculate jacobian. */
    t2 = (fi - fit_data->as_xi[k]);
    for(n=0;n<5;n++){
      jacobian[n] += t2*jt[n];
    }
    
    /* Calculate hessian. */
    for(n=0;n<5;n++){
      for(o=n;o<5;o++){
	hessian[n*5+o] += jt[n]*jt[o];
      }
    }
  }
}


/* 
 * cfCalcJH3DLS()
 *
 * Calculate Jacobian and Hessian for the 3D spline (Least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void cfCalcJH3DLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m,n,o,x_start,y_start,zi;
  double height,fi,rqei,t1;
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
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    l = fit_data->roi_y_index[j];
    m = fit_data->roi_x_index[j];

    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];

    /*
     * This is the LM algorithm according to wikipedia.
     *
     * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
     *
     * fi = rqe_i * fit_fn_i(theta) + variance_i
     * dfi/dtheta = rqe_i * dfit_fn_i(theta)/dtheta
     */
    jt[0] = rqei*peak->psf[j];
    jt[1] = -rqei*height*dxfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[2] = -rqei*height*dyfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[3] = rqei*height*dzfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[4] = rqei;

    /* Calculate jacobian. */
    t1 = (fi - fit_data->x_data[k]);
    for(n=0;n<5;n++){
      jacobian[n] += t1*jt[n];
    }
    
    /* Calculate hessian. */
    for(n=0;n<5;n++){
      for(o=n;o<5;o++){
	hessian[n*5+o] += jt[n]*jt[o];
      }
    }
  }
}


/* 
 * cfCalcJH3DLS()
 *
 * Calculate Jacobian and Hessian for the 3D spline (Fit weighted least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void cfCalcJH3DFWLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m,n,o,x_start,y_start,zi;
  double height,fi,rqei,t1,t2,xi;
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
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    l = fit_data->roi_y_index[j];
    m = fit_data->roi_x_index[j];

    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];

    /*
     * This is the LM algorithm according to wikipedia.
     *
     * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
     *
     * fi = rqe_i * fit_fn_i(theta) + variance_i
     * dfi/dtheta = rqe_i * dfit_fn_i(theta)/dtheta
     */
    jt[0] = rqei*peak->psf[j];
    jt[1] = -rqei*height*dxfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[2] = -rqei*height*dyfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[3] = rqei*height*dzfAt3D(spline_fit->spline_data,zi,l+y_start,m+x_start);
    jt[4] = rqei;

    /* Calculate jacobian. */
    t1 = (fi - xi)/fi;
    for(n=0;n<5;n++){
      jacobian[n] += (2.0*t1 - t1*t1)*jt[n];
    }
    
    /* Calculate hessian. */
    t2 = 1.0/fi;
    for(n=0;n<5;n++){
      for(o=n;o<5;o++){
	hessian[n*5+o] += 2.0*t2*(xi*t2 - t1 + t1*t1)*jt[n]*jt[o];
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
  int i,j,k,x_start,y_start;
  double xd,yd,zd;
  peakData *peak;
  splinePeak *spline_peak;
  splineFit *spline_fit;

  peak = fit_data->working_peak;
  spline_peak = (splinePeak *)peak->peak_model;
  spline_fit = (splineFit *)fit_data->fit_model;

  /*
   * To handle pixel bias issues, the difference between the peak
   * center and the integer center can vary over the range [-0.5, 1.5].
   * This complicates things here as the spline dx and dy need to
   * be in the range [0.0, 1.0]. So we use the spline as if it was
   * one pixel smaller than it actually is and adjust which pixel
   * we start on.
   */
  xd = 1.5 - (peak->params[XCENTER] - peak->xi);
  yd = 1.5 - (peak->params[YCENTER] - peak->yi);
  zd = (peak->params[ZCENTER] - spline_peak->zi);

  x_start = 0;
  if(xd >= 1.0){
    xd -= 1.0;
    x_start = 1;
  }

  y_start = 0;
  if(yd >= 1.0){
    yd -= 1.0;
    y_start = 1;
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
    for(i=0;i<fit_data->roi_n_index;i++){
      j = fit_data->roi_y_index[i];
      k = fit_data->roi_x_index[i];
      peak->psf[i] = fAt3D(spline_fit->spline_data,spline_peak->zi,j+y_start,k+x_start);
    }
  }
  else{
    computeDelta2D(spline_fit->spline_data, yd, xd);
    for(i=0;i<fit_data->roi_n_index;i++){
      j = fit_data->roi_y_index[i];
      k = fit_data->roi_x_index[i];    
      peak->psf[i] = fAt2D(spline_fit->spline_data,j+y_start,k+x_start);
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
    for(i=0;i<fit_data->max_nfit;i++){
      spline_peak = (splinePeak *)(fit_data->fit[i].peak_model);
      free(spline_peak);
    }
  }
  
  spline_peak = (splinePeak *)(fit_data->working_peak->peak_model);
  free(spline_peak);
    
  spline_fit = (splineFit *)fit_data->fit_model;
  splineCleanup(spline_fit->spline_data);
  free(spline_fit);

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
void cfCopyPeak(fitData *fit_data, peakData *original, peakData *copy)
{
  splinePeak *spline_copy, *spline_original;

  spline_copy = (splinePeak *)copy->peak_model;
  spline_original = (splinePeak *)original->peak_model;

  /* Allocate storage, if necessary. */
  if(copy->psf == NULL){
    copy->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
  }
  
  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(fit_data, original, copy);

  /* Copy the parts that are specific to Spliner. */
  spline_copy->x_start = spline_original->x_start;
  spline_copy->y_start = spline_original->y_start;
  spline_copy->zi = spline_original->zi;

  spline_copy->x_delta = spline_original->x_delta;
  spline_copy->y_delta = spline_original->y_delta;
  spline_copy->z_delta = spline_original->z_delta;
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

  for(i=0;i<n_peaks;i++){
    free(peaks[i].peak_model);
  }
}


/*
 * cfInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * spline_data - Pointer to a splineData structure.
 * rqe - Pixel relative quantum efficiency.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* cfInitialize(splineData *spline_data, double *rqe, double *scmos_calibration, double tol, int im_size_x, int im_size_y)
{
  int sx,sy;
  fitData* fit_data;

  fit_data = mFitInitialize(rqe, scmos_calibration, tol, im_size_x, im_size_y);

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
   *  1. Spline sizes are n + 1. 
   *  2. This assumes that the splines are square.
   */
  sx = spline_data->xsize - 1;
  sy = spline_data->ysize - 1;

  /*
   * Since sx, sy are always even not sure if the odd case is 
   * correct.
   */
  if((sx%2)==1){
    fit_data->xoff = (double)(sx/2) - 0.5;
    fit_data->yoff = (double)(sy/2) - 0.5;
  }
  else{
    fit_data->xoff = (double)(sx/2) - 1.0;
    fit_data->yoff = (double)(sy/2) - 1.0;
  }
  
  /*
   * Enforce that the spline is square in X/Y.
   */
  if(sy != sx){
    printf("PSF must be square in X/Y!\n");
    exit(EXIT_FAILURE);
  }
  mFitInitializeROIIndexing(fit_data, sx);
  
  /* Allocate storage for the working peak. */
  fit_data->working_peak->peak_model = (splinePeak *)malloc(sizeof(splinePeak));
  fit_data->working_peak->psf = (double *)malloc(sizeof(double)*sx*sy);

  /* Set function pointers. */
  fit_data->fn_alloc_peaks = &cfAllocPeaks;
  fit_data->fn_calc_peak_shape = &cfCalcPeakShape;
  fit_data->fn_check = &mFitCheck;
  fit_data->fn_copy_peak = &cfCopyPeak;
  fit_data->fn_free_peaks = &cfFreePeaks;
  
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
  fit_data->fn_error_fn = &mFitCalcErr;
  fit_data->fn_update = &cfUpdate3D;
}


/*
 * cfInitialize3DALS()
 *
 * Initializes 3D spline fitting (Anscombe least squares).
 *
 * fit_data - pointer to a fitData structure.
 */
void cfInitialize3DALS(fitData *fit_data)
{
  fit_data->jac_size = 5;
  
  fit_data->fn_calc_JH = &cfCalcJH3DALS;
  fit_data->fn_error_fn = &mFitCalcErrALS;
  fit_data->fn_update = &cfUpdate3D;
}


/*
 * cfInitialize3DLS()
 *
 * Initializes 3D spline fitting (Least squares).
 *
 * fit_data - pointer to a fitData structure.
 */
void cfInitialize3DLS(fitData *fit_data)
{
  fit_data->jac_size = 5;

  fit_data->fn_calc_JH = &cfCalcJH3DLS;
  fit_data->fn_error_fn = &mFitCalcErrLS;
  fit_data->fn_update = &cfUpdate3D;
}


/*
 * cfInitialize3DFWLS()
 *
 * Initializes 3D spline fitting (Fit weighted least squares).
 *
 * fit_data - pointer to a fitData structure.
 */
void cfInitialize3DFWLS(fitData *fit_data)
{
  fit_data->jac_size = 5;

  fit_data->fn_calc_JH = &cfCalcJH3DFWLS;
  fit_data->fn_error_fn = &mFitCalcErrFWLS;
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
  int i,j;
  int start,stop;
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
      
      /* Allocate space for saving the peak shape. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
      }
      
      /* 
       * Calculate (integer) peak locations. 
       *
       * Note: Spliner expects the integer position to be floor() of the floating 
       *       point value, so don't change this to round().
       */
      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);
      spline_peak->zi = (int)peak->params[ZCENTER];
      
      /* Arbitrary initial values for BACKGROUND, HEIGHT. */
      peak->params[BACKGROUND] = 1.0;
      peak->params[HEIGHT] = 1.0;

      /* Copy into working peak. */
      cfCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	cfCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }
      
      /* Calculate peak shape (of working peak). */
      cfCalcPeakShape(fit_data);

      /* Estimate best starting background. */
      mFitEstimatePeakBackground(fit_data);

      if(!strcmp(p_type, "finder")){

	/* Estimate best starting height. */
	mFitEstimatePeakHeight(fit_data);

	/* 
	 * This is for the benefit of multi-plane fitting where some peaks might
	 * not be present in a particular plane, but we don't want to throw them
	 * out here as they presumably are present in other planes.
	 */
	if(fit_data->working_peak->params[HEIGHT] < fit_data->minimum_height){
	  fit_data->working_peak->params[HEIGHT] = fit_data->minimum_height;
	}

	/* Check that the initial height is positive, error it out if not. */
	if(fit_data->working_peak->params[HEIGHT] <= 0.0){
	  printf("%.2f\n", fit_data->working_peak->params[HEIGHT]);
	  printf("Warning peak %d has negative estimated height!\n", (i-start));
	  fit_data->working_peak->status = ERROR;
	  cfCopyPeak(fit_data, fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      cfCopyPeak(fit_data, fit_data->working_peak, peak);      
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
      spline_peak = (splinePeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      peak->params[BACKGROUND] = peak_params[j+3];
      peak->params[HEIGHT] = peak_params[j+4];
      
      /* Allocate space for saving the peak shape. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
      }
      
      /* Calculate (integer) peak locations. */
      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);
      spline_peak->zi = (int)peak->params[ZCENTER];

      /* Copy into working peak. */
      cfCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	cfCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      cfCalcPeakShape(fit_data);
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      cfCopyPeak(fit_data, fit_data->working_peak, peak);
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
    cfCopyPeak(fit_data, peak, fit_data->working_peak);
    fit_data->fn_error_fn(fit_data);
    cfCopyPeak(fit_data, fit_data->working_peak, peak);

    if(VERBOSE){
      printf("newPeak error %d %.6e\n", i, fit_data->working_peak->error);
    }
  }

  fit_data->nfit = stop;
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

  mFitUpdate(peak);
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

  mFitUpdate(peak);

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
