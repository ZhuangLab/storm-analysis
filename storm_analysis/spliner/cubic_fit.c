/*
 * Fit multiple, possible overlapping, cubic splines simultaneously to
 * image data.
 *
 * Hazen 01/14
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cubic_spline.h"
#include "../sa_library/multi_fit.h"


/* Structures */
typedef struct
{
  int zi;                     /* Location of the spline in z. */

  int x_start;                /* Spline offset in x (0 or 1). */
  int y_start;                /* Spline offset in y (0 or 1). */
  
  double x_delta;             /* Peak x delta (0.0 - 1.0). */
  double y_delta;             /* Peak y delta (0.0 - 1.0). */
  double z_delta;             /* Peak z delta (0.0 - 1.0). */
  
  double *peak_values;        /* The peak shape. */
} splinePeak;

  
typedef struct
{
  int fit_type;               /* 2D or 3D spline fit, (S2D or S3D). */

  int spline_size_x;          /* The size of the spline in x (in pixels). */
  int spline_size_y;          /* The size of the spline in y (in pixels). */
  int spline_size_z;          /* The size of the spline in z. */
  
  splineData *spline_data;    /* Spline data structure. */
} splineFit;


/* Functions */
void addPeak(fitData *, peakData *);
void cleanup(fitData *);
void fitDataUpdate(fitData *, peakData *, double *);
fitData* initialize(splineData *, double *, double *, double, int, int);
void iterateSpline(fitData *);
void newPeaks(fitData *, double *, int);
void subtractPeak(fitData *, peakData *);
void updateSpline2D(fitData *, peakData *);
void updateSpline3D(fitData *, peakData *);


/* LAPACK Functions */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);


/*
 * addPeak()
 *
 * Calculate peak shape and add the peak to the 
 * foreground and background data arrays.
 *
 * fit_data - pointer to a fitData structure.
 * peak - pointer to a peakData structure.
 */
void addPeak(fitData *fit_data, peakData *peak)
{
  int j,k,l,m,psx,psy,x_start,y_start;
  double bg,height,xd,yd,zd;
  splinePeak *spline_peak;
  splineFit *spline_fit;
  
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
  height = peak->params[HEIGHT];
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

  
  /* Add peak to the foreground and background arrays. */
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
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
 * cleanup()
 *
 * Frees the fitData structure.
 *
 * fit_data - pointer to a fitData structure.
 */
void cleanup(fitData *fit_data)
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
  spline_fit = (splineFit *)fit_data->fit_model;
  splineCleanup(spline_fit->spline_data);

  free(fit_data->bg_counts);
  free(fit_data->bg_data);
  free(fit_data->f_data);
  free(fit_data->scmos_term);
  free(fit_data->x_data);
  free(fit_data);
}


/*
 * fitDataUpdate()
 *
 * Updates fit data given deltas.
 *
 * Also checks for out-of-bounds parameters.
 *
 * fit_data - pointer to a fitData structure.
 * peak - pointer to the peakData structure to update.
 * delta - the deltas for different parameters.
 */
void fitDataUpdate(fitData *fit_data, peakData *peak, double *delta)
{
  int xi,yi;
  double maxz;
  splineFit *spline_fit;

  spline_fit = (splineFit *)fit_data->fit_model;
  
  /* Update the peak parameters. */
  mFitUpdateParams(peak, delta);

  /* Update peak (integer) location with hysteresis. */
  if(fabs(peak->params[XCENTER] - (double)peak->xi - 0.5) > HYSTERESIS){
    peak->xi = (int)peak->params[XCENTER];
  }
  if(fabs(peak->params[YCENTER] - (double)peak->yi - 0.5) > HYSTERESIS){
    peak->yi = (int)peak->params[YCENTER];
  }
  
  /*
   * Check that the peak hasn't moved to close to the 
   * edge of the image. Flag the peak as bad if it has.
   */
  xi = peak->xi;
  yi = peak->yi;
  if((xi < 0)||(xi >= (fit_data->image_size_x - peak->size_x))||(yi < 0)||(yi >= (fit_data->image_size_y - peak->size_y))){
    peak->status = ERROR;
    fit_data->n_margin++;
    if(TESTING){
      printf("object outside margins, %d, %d, %d\n", peak->index, xi, yi);
    }
  }
  
  /* 
   * Check for negative height. 
   */
  if(peak->params[HEIGHT] < 0.0){
    peak->status = ERROR;
    fit_data->n_neg_height++;
    if(TESTING){
      printf("negative height, %d, %.3f\n", peak->index, peak->params[HEIGHT]);
    }
  }

  /* 
   * Update peak (integer) z location and also check for z out of range. 
   */
  if(spline_fit->fit_type == S3D){
    if(peak->params[ZCENTER] < 1.0e-12){
      peak->params[ZCENTER] = 1.0e-12;
    }
    maxz = ((double)spline_fit->spline_size_z) - 1.0e-12;
    if(peak->params[ZCENTER] > maxz){
      peak->params[ZCENTER] = maxz;
    }
    ((splinePeak *)peak->peak_model)->zi = (int)(peak->params[ZCENTER]);
  }
}


/*
 * initialize()
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
fitData* initialize(splineData *spline_data, double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y)
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

  return fit_data;
}


/*
 * iterateSpline()
 *
 * Performs a single cycle of fit improvement.
 */
void iterateSpline(fitData *fit_data)
{
  int i;
  peakData *peak;

  if(((splineFit *)fit_data->fit_model)->fit_type == S3D){
    for(i=0;i<fit_data->nfit;i++){
      peak = &fit_data->fit[i];
      updateSpline3D(fit_data, peak);
    }
  }
  else{
    for(i=0;i<fit_data->nfit;i++){
      peak = &fit_data->fit[i];
      updateSpline2D(fit_data, peak);
    }
  }

  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    mFitCalcErr(fit_data, peak);
  }
}


/*
 * newPeaks
 *
 * fit_data - Pointer to a fitData structure.
 * peak_params - Input values for the peak parameters.
 * n_peaks - The number of peaks.
 */
void newPeaks(fitData *fit_data, double *peak_params, int n_peaks)
{
  int i,j;
  peakData *peak;
  splinePeak *spline_peak;

  mFitNewPeaks(fit_data);

  /*
   * Free old peaks, if necessary.
   */
  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->nfit;i++){
      spline_peak = ((splinePeak *)fit_data->fit[i].peak_model);
      free(spline_peak->peak_values);
      free(spline_peak);
    }
    free(fit_data->fit);
  }

  /*
   * Initialize peaks (localizations).
   */
  fit_data->nfit = n_peaks;
  fit_data->fit = (peakData *)malloc(sizeof(peakData)*n_peaks);
  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    peak->peak_model = (splinePeak *)malloc(sizeof(splinePeak));
    spline_peak = (splinePeak *)peak->peak_model;

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

    /* Initial location. */
    peak->params[HEIGHT]     = peak_params[i*NPEAKPAR+HEIGHT];
    peak->params[XCENTER]    = peak_params[i*NPEAKPAR+XCENTER] - fit_data->xoff;
    peak->params[YCENTER]    = peak_params[i*NPEAKPAR+YCENTER] - fit_data->yoff;
    peak->params[BACKGROUND] = peak_params[i*NPEAKPAR+BACKGROUND];
    peak->params[ZCENTER]    = peak_params[i*NPEAKPAR+ZCENTER] - fit_data->zoff;

    /* These are not used, but need to be initialized so that they do not come out as confusing garbage. */
    peak->params[XWIDTH] = 0.5;
    peak->params[YWIDTH] = 0.5;
    
    /* Initial clamp values. */
    for(j=0;j<NFITTING;j++){
      peak->clamp[j] = fit_data->clamp_start[j];
      peak->sign[j] = 0;
    }

    /* Spliner specific initializations. */

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
    peak->xi = (int)peak->params[XCENTER];
    peak->yi = (int)peak->params[YCENTER];
    spline_peak->zi = (int)peak->params[ZCENTER];

    /* Calculate peak and add it into the fit. */
    addPeak(fit_data, peak);
  }

  /* Initial error calculation. */
  for(i=0;i<fit_data->nfit;i++){
    mFitCalcErr(fit_data, &fit_data->fit[i]);
  }
}


/*
 * subtractPeak()
 *
 * Subtract the peak out of the current fit, basically 
 * this just undoes addPeak().
 *
 * fit_data - pointer to a fitData structure.
 * peak - pointer to a peakData structure.
 */
void subtractPeak(fitData *fit_data, peakData *peak)
{
  int j,k,l,m;
  double bg,height;
  splinePeak *spline_peak;
  
  spline_peak = (splinePeak *)peak->peak_model;
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  height = peak->params[HEIGHT];
  for (j=0;j<peak->size_x;j++){
    for (k=0;k<peak->size_y;k++){
      m = j * fit_data->image_size_x + k + l;
      fit_data->f_data[m] -= height*spline_peak->peak_values[j*peak->size_x+k];
      fit_data->bg_counts[m] -= 1;
      fit_data->bg_data[m] -= (bg + fit_data->scmos_term[m]);
    }
  }
}


/*
 * updateSpline2D()
 *
 * Update current fit for a 2D spline.
 *
 * fit_data - pointer to a fitData structure.
 * peak - pointer to a peakData structure.
 */
void updateSpline2D(fitData *fit_data, peakData *peak)
{
  // These are for Lapack
  int o = 4, nrhs = 1, lda = 4, ldb = 4, info;

  int i,j,k,l,m,n;
  int x_start, y_start;
  double height,fi,t1,t2,xi;
  double delta[NPEAKPAR];
  double jt[4];
  double jacobian[4];
  double hessian[16];  
  splinePeak *spline_peak;
  splineFit *spline_fit;
  
  if(peak->status==RUNNING){

    /*
     * Initializations.
     */
    spline_peak = (splinePeak *)peak->peak_model;
    spline_fit = (splineFit *)fit_data->fit_model;
    x_start = spline_peak->x_start;
    y_start = spline_peak->y_start;
    
    for(i=0;i<NPEAKPAR;i++){
      delta[i] = 0.0;
    }
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

    /* Subtract the old peak out of the foreground and background arrays. */
    subtractPeak(fit_data, peak);
      
    /* Use Lapack to solve AX=B to calculate update vector. */
    dposv_( "Lower", &o, &nrhs, hessian, &lda, jacobian, &ldb, &info );

    if(info!=0){
      peak->status = ERROR;
      fit_data->n_dposv++;
      if(TESTING){
	printf("fitting error! %d %d %d\n", peak->index, info, ERROR);
      }
    }
    else{
      /* Update params. */
      delta[HEIGHT]     = jacobian[0];
      delta[XCENTER]    = jacobian[1];
      delta[YCENTER]    = jacobian[2];
      delta[BACKGROUND] = jacobian[3];
      fitDataUpdate(fit_data, peak, delta);

      /* Add the new peak to the foreground and background arrays. */
      if (peak->status != ERROR){
	addPeak(fit_data, peak);
      }
    }
  }
}


/*
 * updateSpline3D()
 *
 * Update current fit for a 3D spline.
 *
 * fit_data - pointer to a fitData structure.
 * peak - pointer to a peakData structure.
 */
void updateSpline3D(fitData *fit_data, peakData *peak)
{
  // These are for Lapack
  int o = 5, nrhs = 1, lda = 5, ldb = 5, info;

  int i,j,k,l,m,n,zi;
  int x_start, y_start;
  double height,fi,t1,t2,xi;
  double delta[NPEAKPAR];
  double jt[5];
  double jacobian[5];
  double hessian[25];  
  splinePeak *spline_peak;
  splineFit *spline_fit;
  
  if(peak->status==RUNNING){

    /*
     * Initializations.
     */
    spline_peak = (splinePeak *)peak->peak_model;
    spline_fit = (splineFit *)fit_data->fit_model;
      
    x_start = spline_peak->x_start;
    y_start = spline_peak->y_start;
    zi = spline_peak->zi;
    
    for(i=0;i<NPEAKPAR;i++){
      delta[i] = 0.0;
    }
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

    /* Subtract the old peak out of the foreground and background arrays. */
    subtractPeak(fit_data, peak);
      
    /* Use Lapack to solve AX=B to calculate update vector. */
    dposv_( "Lower", &o, &nrhs, hessian, &lda, jacobian, &ldb, &info );

    if(info!=0){
      peak->status = ERROR;
      fit_data->n_dposv++;
      if(TESTING){
	printf("fitting error! %d %d %d\n", peak->index, info, ERROR);
      }
    }
    else{
      /* Update params. */
      delta[HEIGHT]     = jacobian[0];
      delta[XCENTER]    = jacobian[1];
      delta[YCENTER]    = jacobian[2];
      delta[ZCENTER]    = jacobian[3];
      delta[BACKGROUND] = jacobian[4]; 
      fitDataUpdate(fit_data, peak, delta);

      /* Add the new peak to the foreground and background arrays. */
      if (peak->status != ERROR){
	addPeak(fit_data, peak);
      }
    }
  }
}
