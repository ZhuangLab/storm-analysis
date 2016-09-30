/*
 * Routine(s) for attempting to (MLE) fit multiple gaussians to an image.
 * The approach follows Laurence and Chromy, Nature Methods, 2010.
 *
 * 01/11
 *
 *
 * Generalized for 2D & 3D fitting.
 *
 * 07/11
 *
 *
 * Speed things up by only updating foreground and background data
 * for those peaks that are still running (i.e. moving around).
 *
 * 10/11
 *
 *
 * Remove (useless) OFFSET term.
 *
 * 03/13
 *
 *
 * Add hystersis to minimize a bad interaction between the parameter
 * clamp and moving / changing the size of the AOI.
 *
 * 07/13
 *
 *
 * Add scmos_term to enable analysis of sCMOS data.
 *
 * 10/13
 *
 *
 * Remove static variables so that it is thread safe.
 *
 * 09/16
 *
 * Hazen
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall multi_fit.c
 *  gcc -shared -Wl,-soname,multi_fit.so.1 -o multi_fit.so.1.0.1 multi_fit.o -lc -llapack
 *
 * Windows:
 *  gcc -c multi_fit.c
 *  gcc -shared -o multi_fit.dll multi_fit.o -llapack -Lc:\Users\Hazen\lib
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "multi_fit.h"

/* Define */
#define TESTING 0
#define VERBOSE 0

#define HYSTERESIS 0.6 /* In order to move the AOI or change it's size,
			  the new value must differ from the old value
			  by at least this much (<= 0.5 is no hysteresis). */

#define MARGIN 10      /* Margin around the edge of the image. This
			  also sets the maximum number of pixels that will
			  be fit.
			  
			  FIXME: This should be adjustable. */


/* 
 * Structures
 *
 * I read somewhere (and long ago) that it is better to organize 
 * structures with fixed types first and pointers at the end, so 
 * that is the convention I try to follow.
 */

typedef struct
{
  int offset;
  int status;
  int wx;
  int wy;
  int xc;
  int yc;
  double error;
  double error_old;
  double wx_term;
  double wy_term;

  int sign[NFITTING];

  double clamp[NFITTING];
  double params[NFITTING];  /* [height x-center x-width y-center y-width background] */  
  double xt[2*MARGIN+1];
  double ext[2*MARGIN+1];
  double yt[2*MARGIN+1];
  double eyt[2*MARGIN+1];
} peakData;

typedef struct
{
  int nfit;                 /* number of peaks to fit. */
  int image_size_x;         /* size in x (fast axis). */
  int image_size_y;         /* size in y (slow axis). */
  int zfit;                 /* fit with wx, wy as fixed functions of z. */

  /* These are for diagnostics. */
  int n_dposv;              /* number lost to an error trying to solve Ax = b. */
  int n_margin;             /* number lost because they were too close to the edge of the image. */
  int n_neg_fi;             /* number lost to a negative fi. */
  int n_neg_height;         /* number lost to negative height. */
  int n_neg_width;          /* number lost to negative width. */

  double tolerance;         /* fit tolerance. */
  double min_z;             /* minimum z value. */
  double max_z;             /* maximum z value. */

  int *bg_counts;           /* number of peaks covering a particular pixel. */
  
  double *bg_data;          /* background data. */
  double *f_data;           /* fit (foreground) data. */
  double *scmos_term;       /* sCMOS calibration term for each pixel (var/gain^2). */
  double *x_data;           /* image data. */

  double clamp_start[7];    /* starting values for the peak clamp values. */
  double wx_z_params[5];    /* x width versus z parameters. */
  double wy_z_params[5];    /* y width versus z parameters. */

  peakData *fit;
} fitData;


/* Function Declarations */
void addPeak(fitData *, peakData *);
void calcErr(fitData *);
void calcFit(fitData *);
int calcWidth(double, int);
void calcWidthsFromZ(fitData *, peakData *);
void cleanup(fitData *);
void fitDataUpdate(fitData *, peakData *, double *);
//double getError(fitData *);
void getResidual(fitData *, double *);
void getResults(fitData *, double *);
int getUnconverged(fitData *);
fitData* initialize(double *, double *, double, int, int);
void initializeZParameters(fitData *, double *, double *, double, double);
void iterate2DFixed(fitData *);
void iterate2D(fitData *);
void iterate3D(fitData *);
void iterateZ(fitData *);
void newImage(fitData *, double *);
void newPeaks(fitData *, double *, int);
void subtractPeak(fitData *, peakData *);
void update2DFixed(fitData *);
void update2D(fitData *);
void update3D(fitData *);
void updateZ(fitData *);

/* LAPACK Functions */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);



/* Functions */


/*
 * addPeak()
 *
 * Add a peak to the foreground and background data arrays.
 *
 * fit_data - pointer to a fitData structure.
 * peak - pointer to a peakData structure.
 */
void addPeak(fitData *fit_data, peakData *peak)
{
  int j,k,l,m,n,wx,wy,xc,yc;
  double bg,mag,tmp,xt,yt;

  xc = peak->xc;
  yc = peak->yc;
  peak->offset = yc * fit_data->image_size_x + xc;

  wx = peak->wx;
  wy = peak->wy;

  for(j=(xc-wx);j<=(xc+wx);j++){
    xt = (double)j - peak->params[XCENTER];
    n = j-xc+wx;
    peak->xt[n] = xt;
    peak->ext[n] = exp(-xt*xt*peak->params[XWIDTH]);
  }
  for(j=(yc-wy);j<=(yc+wy);j++){
    yt = (double)j - peak->params[YCENTER];
    n = j-yc+wy;
    peak->yt[n] = yt;
    peak->eyt[n] = exp(-yt*yt*peak->params[YWIDTH]);
  }

  /* gaussian function */
  l = peak->offset;
  bg = peak->params[BACKGROUND];
  mag = peak->params[HEIGHT];
  for(j=-wy;j<=wy;j++){
    tmp = peak->eyt[j+wy];
    for(k=-wx;k<=wx;k++){
      m = j * fit_data->image_size_x + k + l;
      fit_data->f_data[m] += mag * tmp * peak->ext[k+wx];
      fit_data->bg_counts[m] += 1;
      fit_data->bg_data[m] += bg + fit_data->scmos_term[m];
    }
  }

}


/*
 * calcErr()
 *
 * Calculate error in fit for each gaussian.
 *
 * fit_data - pointer to a fitData structure.
 */
void calcErr(fitData *fit_data)
{
  int i,j,k,l,m,wx,wy;
  double err,fi,xi;
  peakData *peak;

  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    if(peak->status == RUNNING){
      l = peak->offset;
      wx = peak->wx;
      wy = peak->wy;
      err = 0.0;
      for(j=-wy;j<=wy;j++){
	for(k=-wx;k<=wx;k++){
	  m = (j * fit_data->image_size_x) + k + l;
	  fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
	  if (fi <= 0.0){
	    if(VERBOSE){
	      printf(" Negative f detected! %.3f %.3f %.3f %.3f %d\n", peak->params[BACKGROUND], fi, fit_data->f_data[m], fit_data->bg_data[m], fit_data->bg_counts[m]);
	    }
	    peak->status = ERROR;
	    fit_data->n_neg_fi++;
	    j = wy + 1;
	    k = wx + 1;
	  }
	  xi = fit_data->x_data[m];
	  if (xi <= 0.0){
	    if(VERBOSE){
	      printf(" Negative x detected! %.3f\n", xi);
	    }
	  }
	  err += 2*(fi-xi)-2*xi*log(fi/xi);
	}
      }
      peak->error_old = peak->error;
      peak->error = err;
      if (VERBOSE){
	printf("%d %f %f %f\n", i, peak->error_old, peak->error, fit_data->tolerance);
      }
      if(((fabs(err - peak->error_old)/err) < fit_data->tolerance) && (peak->status != ERROR)){
	peak->status = CONVERGED;
      }
    }
  }
}


/*
 * calcFit()
 *
 * Calculate fit from gaussian parameters. This assumes
 * that all the peak parameters are reasonable.
 */
void calcFit(fitData *fit_data)
{
  int i;
  peakData *peak;

  // zero matrices.
  for(i=0;i<(fit_data->image_size_x * fit_data->image_size_y);i++){
    fit_data->f_data[i] = 1.0;
    fit_data->bg_counts[i] = 0;
    fit_data->bg_data[i] = 0;
  }

  // update fit matrix with values from fits.
  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    if (peak->status != ERROR){
      addPeak(fit_data, peak);
    }
  }
}


/*
 * calcWidth(peak_width)
 *
 * Given a peak_width, returns the appropriate 
 * bounding box to use for fitting.
 */
int calcWidth(double peak_width, int old_w)
{
  int new_w;
  double tmp;

  if(peak_width < 0.0){
    if(TESTING){
      printf(" Got negative peak width! %.3f", peak_width);
    }
    return 1;
  }
  else{
    new_w = old_w;
    tmp = 4.0*sqrt(1.0/(2.0*peak_width));
    if(fabs(tmp - (double)old_w - 0.5) > HYSTERESIS){
      new_w = (int)tmp;
    }
    if(new_w > MARGIN){
      new_w = MARGIN;
    }
    return new_w;
  }
}


/*
 * calcWidthsFromZ(cur)
 *
 * Updates wx, wy given z.
 *
 * fit_data - The fitD
 * peak - Peak data structure to update.
 */
void calcWidthsFromZ(fitData *fit_data, peakData *peak)
{
  double z0,z1,z2,z3,tmp;

  // wx
  z0 = (peak->params[ZCENTER] - fit_data->wx_z_params[1]) / fit_data->wx_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  z3 = z2*z0;
  tmp = 1.0 + z1 + fit_data->wx_z_params[3] * z2 + fit_data->wx_z_params[4] * z3;
  peak->wx_term = tmp*tmp;
  peak->params[XWIDTH] = 2.0/(fit_data->wx_z_params[0] * tmp);

  // wy
  z0 = (peak->params[ZCENTER] - fit_data->wy_z_params[1]) / fit_data->wy_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  z3 = z2*z0;
  tmp = 1.0 + z1 + fit_data->wy_z_params[3] * z2 + fit_data->wy_z_params[4] * z3;
  peak->wy_term = tmp*tmp;
  peak->params[YWIDTH] = 2.0/(fit_data->wy_z_params[0] * tmp);
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
  free(fit_data->fit);
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
  int i,xc,yc;

  // update
  for(i=0;i<NFITTING;i++){
    if(VERBOSE){
      printf("%.3e %.3f | ", delta[i], peak->clamp[i]);
    }

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

    // update values based on delta & clamp.
    if (delta[i] != 0.0){
      peak->params[i] -= delta[i]/(1.0 + fabs(delta[i])/peak->clamp[i]);
    }
  }
  if(VERBOSE){
    printf("\n");
  }

  // Update peak (integer) center (w/ hysteresis).
  if(fabs(peak->params[XCENTER] - (double)peak->xc - 0.5) > HYSTERESIS){
    peak->xc = (int)peak->params[XCENTER];
  }
  if(fabs(peak->params[YCENTER] - (double)peak->yc - 0.5) > HYSTERESIS){
    peak->yc = (int)peak->params[YCENTER];
  }

  // Check that the peak hasn't moved to close to the 
  // edge of the image. Flag the peak as bad if it has.
  xc = peak->xc;
  yc = peak->yc;
  if((xc <= MARGIN)||(xc >= (fit_data->image_size_x - MARGIN-1))||(yc <= MARGIN)||(yc >= (fit_data->image_size_y - MARGIN - 1))){
    peak->status = ERROR;
    fit_data->n_margin++;
    if(TESTING){
      printf("object outside margins, %.3f, %.3f\n", peak->params[XCENTER], peak->params[YCENTER]);
    }
  }
  
  // check for negative background or height
  if(peak->params[HEIGHT]<0.0){
    peak->status = ERROR;
    fit_data->n_neg_height++;
    if(TESTING){
      printf("negative height, %.3f, %.3f (%.3f, %.3f)\n", peak->params[BACKGROUND], peak->params[HEIGHT], peak->params[XCENTER], peak->params[YCENTER]);
    }
  }

  // check for negative widths
  if((peak->params[XWIDTH]<0.0)||(peak->params[YWIDTH]<0.0)){
    peak->status = ERROR;
    fit_data->n_neg_width++;
    if(TESTING){
      printf("negative widths, %.3f, %.3f (%.3f, %.3f)\n", peak->params[XWIDTH], peak->params[YWIDTH], peak->params[XCENTER], peak->params[YCENTER]);
    }
  }

  // Option 1: Peak errors out if z is out of range.
  /*
  if((cur->params[ZCENTER]<MINZ)||(cur->params[ZCENTER]>MAXZ)){
    cur->status = ERROR;
    if(TESTING){
      printf("z value out of range, %.3f (%.3f, %.3f)\n", cur->params[ZCENTER], cur->params[XCENTER], cur->params[YCENTER]);
    }
  }
  */

  // Option 2: Clamp z value range.
  if (fit_data->zfit){
    if(peak->params[ZCENTER] < fit_data->min_z){
      peak->params[ZCENTER] = fit_data->min_z;
    }

    if(peak->params[ZCENTER] > fit_data->max_z){
      peak->params[ZCENTER] = fit_data->max_z;
    }
  }

}


/*
 * getError()
 *
 * Return the current error in the fit.
 */
/* FIXME: not used, remove?
double getError()
{
  int i;
  double err;
  fitData *cur;

  err = 0.0;
  for(i=0;i<nfit;i++){
    cur = &fit[i];
    err += cur->error;
  }

  return err;
}
*/

/*
 * getResidual(residual).
 *
 * Returns image - fit.
 *
 * fit_data - Pointer to a fitData structure.
 * residual - Pre-allocated space to store the residual values.
 *            This should be square & the same size as the image.
 */
void getResidual(fitData *fit_data, double *residual)
{
  int i;

  calcFit(fit_data);
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
void getResults(fitData *fit_data, double *peak_params)
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
    peak_params[i*NPEAKPAR+XCENTER]    = peak->params[XCENTER];
    peak_params[i*NPEAKPAR+YCENTER]    = peak->params[YCENTER];
    peak_params[i*NPEAKPAR+BACKGROUND] = peak->params[BACKGROUND];
    peak_params[i*NPEAKPAR+ZCENTER]    = peak->params[ZCENTER];

    peak_params[i*NPEAKPAR+STATUS] = (double)peak->status;
    peak_params[i*NPEAKPAR+IERROR] = peak->error;
  }
}


/*
 * getUnconverged()
 *
 * Return the number of fits that have not yet converged.
 *
 * fit_data - Pointer to a fitData structure.
 */
int getUnconverged(fitData *fit_data)
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
 * initialize(image, params, tol, im_size, n)
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
fitData* initialize(double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y)
{
  int i;
  fitData* fit_data;

  /* Initialize fitData structure. */
  fit_data = (fitData*)malloc(sizeof(fitData));
  fit_data->image_size_x = im_size_x;
  fit_data->image_size_y = im_size_y;
  fit_data->tolerance = tol;
  fit_data->zfit = 0;
  
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

  /* 
   * Default z fitting range. This is only relevant for Z mode peak 
   * fitting. These values are set to the correct values by calling
   * initializeZParameters().
   */
  fit_data->min_z = -1.0e-6;
  fit_data->max_z = 1.0e+6;

  return fit_data;
}

/*
 * initializeZParameters(wx_vs_z, wy_vs_z)
 *
 * Initializes fitting for z with wx, wy dependence on z.
 *
 * fit_data - pointer to a fitData structure.
 * wx_vs_z - [wo, c, d, A, B] wx vs z coefficients.
 * wy_vs_z - [wo, c, d, A, B] wy vs z coefficients.
 * z_min - minimum allowed z value.
 * z_max - maximum allowed z value.
 */
void initializeZParameters(fitData* fit_data, double *wx_vs_z, double *wy_vs_z, double z_min, double z_max)
{
  int i;

  fit_data->zfit = 1;
  for(i=0;i<5;i++){
    fit_data->wx_z_params[i] = wx_vs_z[i];
    fit_data->wy_z_params[i] = wy_vs_z[i];
  }
  fit_data->wx_z_params[0] = fit_data->wx_z_params[0] * fit_data->wx_z_params[0];
  fit_data->wy_z_params[0] = fit_data->wy_z_params[0] * fit_data->wy_z_params[0];

  fit_data->min_z = z_min;
  fit_data->max_z = z_max;
}


/*
 * iterate2DFixed()
 *
 * Performs a single cycle of fit improvement with fixed x, y width.
 *
 */
void iterate2DFixed(fitData *fit_data)
{
  update2DFixed(fit_data);
  calcErr(fit_data);
}


/*
 * iterate2D()
 *
 * Performs a single cycle of fit improvement with equal x, y width.
 *
 */
void iterate2D(fitData *fit_data)
{
  update2D(fit_data);
  calcErr(fit_data);
}


/*
 * iterate3D()
 *
 * Performs a single cycle of fit improvement.
 */
void iterate3D(fitData *fit_data)
{
  update3D(fit_data);
  calcErr(fit_data);
}


/*
 * iterateZ()
 *
 * Performs a single cycle of fit improvement with x, y width
 * determined by the z parameter.
 */
void iterateZ(fitData *fit_data)
{
  updateZ(fit_data);
  calcErr(fit_data);
}


/*
 * newImage
 *
 * fit_data - Pointer to a fitData structure.
 * new_image - Pointer to the image data of size image_size_x by image_size_y.
 */
void newImage(fitData *fit_data, double *new_image)
{
  int i;

  for(i=0;i<(fit_data->image_size_x*fit_data->image_size_y);i++){
    fit_data->x_data[i] = new_image[i];
  }  
}


/*
 * newPeaks
 *
 * fit_data - Pointer to a fitData structure.
 * peak_data - Input values for the peak parameters.
 * n_peaks - The number of peaks.
 */
void newPeaks(fitData *fit_data, double *peak_data, int n_peaks)
{
  int i,j;
  peakData *peak;

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

  /*
   * Initialize peaks (localizations).
   */
  fit_data->nfit = n_peaks;
  fit_data->fit = (peakData *)malloc(sizeof(peakData)*n_peaks);
  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    peak->status = (int)(peak_data[i*NPEAKPAR+STATUS]);
    if(peak->status==RUNNING){
      peak->error = 0.0;
      peak->error_old = 0.0;
    }
    else {
      peak->error = peak_data[i*NPEAKPAR+IERROR];
      peak->error_old = peak->error;
    }
    
    peak->params[HEIGHT]     = peak_data[i*NPEAKPAR+HEIGHT];
    peak->params[XCENTER]    = peak_data[i*NPEAKPAR+XCENTER];
    peak->params[YCENTER]    = peak_data[i*NPEAKPAR+YCENTER];
    peak->params[BACKGROUND] = peak_data[i*NPEAKPAR+BACKGROUND];
    peak->params[ZCENTER]    = peak_data[i*NPEAKPAR+ZCENTER];

    if(fit_data->zfit){
      calcWidthsFromZ(fit_data, peak);
    }
    else{
      peak->params[XWIDTH] = 1.0/(2.0*peak_data[i*NPEAKPAR+XWIDTH]*peak_data[i*NPEAKPAR+XWIDTH]);
      peak->params[YWIDTH] = 1.0/(2.0*peak_data[i*NPEAKPAR+YWIDTH]*peak_data[i*NPEAKPAR+YWIDTH]);
    }

    peak->xc = (int)peak->params[XCENTER];
    peak->yc = (int)peak->params[YCENTER];
    peak->wx = calcWidth(peak->params[XWIDTH],-10.0);
    peak->wy = calcWidth(peak->params[YWIDTH],-10.0);

    for(j=0;j<NFITTING;j++){
      peak->sign[j] = 0;
    }
  }

  calcFit(fit_data);
  calcErr(fit_data);
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
  int j,k,l,m,wx,wy;
  double bg,mag,tmp;

  wx = peak->wx;
  wy = peak->wy;

  /* gaussian function */
  l = peak->offset;
  bg = peak->params[BACKGROUND];
  mag = peak->params[HEIGHT];
  for(j=-wy;j<=wy;j++){
    tmp = peak->eyt[j+wy];
    for(k=-wx;k<=wx;k++){
      m = j * fit_data->image_size_x + k + l;
      fit_data->f_data[m] -= mag * tmp * peak->ext[k+wx];
      fit_data->bg_counts[m] -= 1;
      fit_data->bg_data[m] -= (bg + fit_data->scmos_term[m]);
    }
  } 
}


/*
 * update2DFixed()
 *
 * Update current fits given fixed x & y width.
 *
 * This procedure is also responsible for flagging peaks
 * that might be bad & that should be removed from fitting.
 *
 * fit_data - pointer to a fitData structure.
 */
void update2DFixed(fitData *fit_data)
{
  // Lapack
  int n = 4, nrhs = 1, lda = 4, ldb = 4, info;

  // Local
  int i,j,k,l,m,wx,wy;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double delta[NPEAKPAR];
  double jt[4];
  double jacobian[4];
  double hessian[16];
  peakData *peak;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    if(peak->status==RUNNING){
      for(j=0;j<4;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<16;j++){
	hessian[j] = 0.0;
      }
      l = peak->offset;
      wx = peak->wx;
      wy = peak->wy;
      a1 = peak->params[HEIGHT];
      width = peak->params[XWIDTH];
      for(j=-wy;j<=wy;j++){
	yt = peak->yt[j+wy];
	eyt = peak->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j * fit_data->image_size_x + k + l;
	  fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
	  xi = fit_data->x_data[m];
	  xt = peak->xt[k+wx];
	  ext = peak->ext[k+wx];
	  e_t = ext*eyt;

	  jt[0] = e_t;
	  jt[1] = 2.0*a1*width*xt*e_t;
	  jt[2] = 2.0*a1*width*yt*e_t;
	  jt[3] = 1.0;
	  
	  // calculate jacobian
	  t1 = 2.0*(1.0 - xi/fi);
	  jacobian[0] += t1*jt[0];
	  jacobian[1] += t1*jt[1];
	  jacobian[2] += t1*jt[2];
	  jacobian[3] += t1*jt[3];
	  
	  // calculate hessian
	  t2 = 2.0*xi/(fi*fi);

	  if (0){
	    // hessian with second derivative terms.
	    hessian[0] += t2*jt[0]*jt[0];
	    hessian[1] += t2*jt[0]*jt[1]+t1*2.0*xt*width*e_t;
	    hessian[2] += t2*jt[0]*jt[2]+t1*2.0*yt*width*e_t;
	    hessian[3] += t2*jt[0]*jt[3];
	    
	    // hessian[4]
	    hessian[5] += t2*jt[1]*jt[1]+t1*(-2.0*a1*width*e_t+4.0*a1*width*width*xt*xt*e_t);
	    hessian[6] += t2*jt[1]*jt[2]+t1*4.0*a1*xt*yt*width*width*e_t;
	    hessian[7] += t2*jt[1]*jt[3];
	    
	    // hessian[8]
	    // hessian[9]
	    hessian[10] += t2*jt[2]*jt[2]+t1*(-2.0*a1*width*e_t+4.0*a1*width*width*yt*yt*e_t);
	    hessian[11] += t2*jt[2]*jt[3];
	    
	    // hessian[12]
	    // hessian[13]
	    // hessian[14]
	    hessian[15]  += t2*jt[3]*jt[3];
	  }
	  else{
	    // calculate hessian without second derivative terms.
	    hessian[0] += t2*jt[0]*jt[0];
	    hessian[1] += t2*jt[0]*jt[1];
	    hessian[2] += t2*jt[0]*jt[2];
	    hessian[3] += t2*jt[0]*jt[3];
	    
	    // hessian[4]
	    hessian[5] += t2*jt[1]*jt[1];
	    hessian[6] += t2*jt[1]*jt[2];
	    hessian[7] += t2*jt[1]*jt[3];
	    
	    // hessian[8]
	    // hessian[9]
	    hessian[10] += t2*jt[2]*jt[2];
	    hessian[11] += t2*jt[2]*jt[3];
	    
	    // hessian[12]
	    // hessian[13]
	    // hessian[14]
	    hessian[15] += t2*jt[3]*jt[3];
	  }
	}
      }

      // subtract the old peak out of the foreground and background arrays.
      subtractPeak(fit_data, peak);
      
      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      if(info!=0){
	peak->status = ERROR;
	fit_data->n_dposv++;
	if(TESTING){
	  printf("fitting error! %d %d %d\n", i, info, ERROR);
	}
      }
      else{
	// update params
	delta[HEIGHT]     = jacobian[0];
	delta[XCENTER]    = jacobian[1];
	delta[YCENTER]    = jacobian[2];
	delta[BACKGROUND] = jacobian[3];

	if(VERBOSE){
	  printf("%d\n", i);
	}
	fitDataUpdate(fit_data, peak, delta);

	// add the new peak to the foreground and background arrays.
	if (peak->status != ERROR){
	  addPeak(fit_data, peak);
	}
      }
    }
  }
  if(VERBOSE){
    printf("\n");
  }
}


/*
 * update2D()
 *
 * Update current fits given equal width in x and y.
 *
 * This procedure is also responsible for flagging peaks
 * that might be bad & that should be removed from fitting.
 *
 * fit_data - pointer to a fitData structure.
 */
void update2D(fitData *fit_data)
{
  // Lapack
  int n = 5, nrhs = 1, lda = 5, ldb = 5, info;

  // Local
  int i,j,k,l,m,wx,wy;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double delta[NPEAKPAR];
  double jt[5];
  double jacobian[5];
  double hessian[25];
  peakData *peak;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    if(peak->status==RUNNING){
      for(j=0;j<5;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<25;j++){
	hessian[j] = 0.0;
      }
      l = peak->offset;
      wx = peak->wx;
      wy = peak->wy;
      a1 = peak->params[HEIGHT];
      width = peak->params[XWIDTH];
      for(j=-wy;j<=wy;j++){
	yt = peak->yt[j+wy];
	eyt = peak->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j * fit_data->image_size_x + k + l;
	  fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
	  xi = fit_data->x_data[m];
	  xt = peak->xt[k+wx];
	  ext = peak->ext[k+wx];
	  e_t = ext*eyt;

	  jt[0] = e_t;
	  jt[1] = 2.0*a1*width*xt*e_t;
	  jt[2] = 2.0*a1*width*yt*e_t;
	  jt[3] = -a1*xt*xt*e_t-a1*yt*yt*e_t;
	  jt[4] = 1.0;
	  
	  // calculate jacobian
	  t1 = 2.0*(1.0 - xi/fi);
	  jacobian[0] += t1*jt[0];
	  jacobian[1] += t1*jt[1];
	  jacobian[2] += t1*jt[2];
	  jacobian[3] += t1*jt[3];
	  jacobian[4] += t1*jt[4];
	  
	  // calculate hessian
	  t2 = 2.0*xi/(fi*fi);

	  // calculate hessian without second derivative terms.
	  hessian[0] += t2*jt[0]*jt[0];
	  hessian[1] += t2*jt[0]*jt[1];
	  hessian[2] += t2*jt[0]*jt[2];
	  hessian[3] += t2*jt[0]*jt[3];
	  hessian[4] += t2*jt[0]*jt[4];
	  
	  // hessian[5]
	  hessian[6] += t2*jt[1]*jt[1];
	  hessian[7] += t2*jt[1]*jt[2];
	  hessian[8] += t2*jt[1]*jt[3];
	  hessian[9] += t2*jt[1]*jt[4];
	    
	  // hessian[10]
	  // hessian[11]
	  hessian[12] += t2*jt[2]*jt[2];
	  hessian[13] += t2*jt[2]*jt[3];
	  hessian[14] += t2*jt[2]*jt[4];
	  
	  // hessian[15]
	  // hessian[16]
	  // hessian[17]
	  hessian[18] += t2*jt[3]*jt[3];
	  hessian[19] += t2*jt[3]*jt[4];

	  // hessian[20]
	  // hessian[21]
	  // hessian[22]
	  // hessian[23]
	  hessian[24] += t2*jt[4]*jt[4];
	}
      }
      

      // subtract the old peak out of the foreground and background arrays.
      subtractPeak(fit_data, peak);

      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      if(info!=0){
	peak->status = ERROR;
	fit_data->n_dposv++;
	if(TESTING){
	  printf("fitting error! %d %d %d\n", i, info, ERROR);
	  printf("  %f %f %f %f %f\n", delta[HEIGHT], delta[XCENTER], delta[YCENTER], delta[XWIDTH], delta[BACKGROUND]);
	}
      }
      else{
	// update params
	delta[HEIGHT]     = jacobian[0];
	delta[XCENTER]    = jacobian[1];
	delta[YCENTER]    = jacobian[2];
	delta[XWIDTH]     = jacobian[3];
	delta[YWIDTH]     = jacobian[3];
	delta[BACKGROUND] = jacobian[4];

	fitDataUpdate(fit_data, peak, delta);

	// add the new peak to the foreground and background arrays.
	// recalculate peak fit area as the peak width may have changed.
	if (peak->status != ERROR){
	  peak->wx = calcWidth(peak->params[XWIDTH], peak->wx);
	  peak->wy = peak->wx;
	  addPeak(fit_data, peak);
	}
      }
    }
  }
}


/*
 * update3D()
 *
 * Update current fits allowing all parameters to change.
 *
 * This procedure is also responsible for flagging peaks
 * that might be bad & that should be removed from fitting.
 *
 * fit_data - pointer to a fitData structure.
 */
void update3D(fitData *fit_data)
{
  // Lapack
  int n = 6, nrhs = 1, lda = 6, ldb = 6, info;

  // Local
  int i,j,k,l,m,wx,wy;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,a3,a5;
  double delta[NPEAKPAR];
  double jt[6];
  double jacobian[6];
  double hessian[36];
  peakData *peak;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    if(peak->status==RUNNING){
      for(j=0;j<6;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<36;j++){
	hessian[j] = 0.0;
      }
      l = peak->offset;
      wx = peak->wx;
      wy = peak->wy;
      a1 = peak->params[HEIGHT];
      a3 = peak->params[XWIDTH];
      a5 = peak->params[YWIDTH];
      for(j=-wy;j<=wy;j++){
	yt = peak->yt[j+wy];
	eyt = peak->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j * fit_data->image_size_x + k + l;
	  fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
	  xi = fit_data->x_data[m];
	  xt = peak->xt[k+wx];
	  ext = peak->ext[k+wx];
	  e_t = ext*eyt;
	  
	  jt[0] = e_t;
	  jt[1] = 2.0*a1*a3*xt*e_t;
	  jt[2] = -a1*xt*xt*e_t;
	  jt[3] = 2.0*a1*a5*yt*e_t;
	  jt[4] = -a1*yt*yt*e_t;
	  jt[5] = 1.0;
	    	  
	  // calculate jacobian
	  t1 = 2.0*(1.0 - xi/fi);
	  jacobian[0] += t1*jt[0];
	  jacobian[1] += t1*jt[1];
	  jacobian[2] += t1*jt[2];
	  jacobian[3] += t1*jt[3];
	  jacobian[4] += t1*jt[4];
	  jacobian[5] += t1*jt[5];

	  // calculate hessian
	  t2 = 2.0*xi/(fi*fi);

	  if (0){
	    // FIXME: not complete
	    // hessian with second derivative terms.
	    hessian[0] += t2*jt[0]*jt[0];
	    hessian[1] += t2*jt[0]*jt[1]+t1*2.0*xt*a3*ext*eyt;
	    hessian[2] += t2*jt[0]*jt[2];
	    hessian[3] += t2*jt[0]*jt[3]+t1*2.0*yt*a5*ext*eyt;
	    hessian[4] += t2*jt[0]*jt[4];
	    hessian[5] += t2*jt[0]*jt[5];
	    
	    // hessian[6]
	    hessian[7]  += t2*jt[1]*jt[1]+t1*(-2.0*a1*a3*ext*eyt+4.0*a1*a3*a3*xt*xt*ext*eyt);
	    hessian[8]  += t2*jt[1]*jt[2]+t1*4.0*a1*xt*yt*a3*a3*ext*eyt;
	    hessian[9]  += t2*jt[1]*jt[3];
	    hessian[10] += t2*jt[1]*jt[4];
	    hessian[11] += t2*jt[1]*jt[5];
	    
	    // hessian[12]
	    // hessian[13]
	    hessian[14] += t2*jt[2]*jt[2]+t1*(-2.0*a1*a3*ext*eyt+4.0*a1*a3*a3*yt*yt*ext*eyt);
	    hessian[15] += t2*jt[2]*jt[3];
	    hessian[16] += t2*jt[2]*jt[4];
	    hessian[17] += t2*jt[2]*jt[5];
	    
	    // hessian[18]
	    // hessian[19]
	    // hessian[20]
	    hessian[21] += t2*jt[3]*jt[3];
	    hessian[22] += t2*jt[3]*jt[4];
	    hessian[23] += t2*jt[3]*jt[5];
	    
	    // hessian[24]
	    // hessian[25]
	    // hessian[26]
	    // hessian[27]
	    hessian[28] += t2*jt[4]*jt[4];
	    hessian[29] += t2*jt[4]*jt[5];

	    // hessian[30]
	    // hessian[31]
	    // hessian[32]
	    // hessian[33]
	    // hessian[34]
	    hessian[35] += t2*jt[5]*jt[5];
	  }
	  else {
	    // hessian without second derivative terms.
	    hessian[0] += t2*jt[0]*jt[0];
	    hessian[1] += t2*jt[0]*jt[1];
	    hessian[2] += t2*jt[0]*jt[2];
	    hessian[3] += t2*jt[0]*jt[3];
	    hessian[4] += t2*jt[0]*jt[4];
	    hessian[5] += t2*jt[0]*jt[5];
	    
	    // hessian[6]
	    hessian[7]  += t2*jt[1]*jt[1];
	    hessian[8]  += t2*jt[1]*jt[2];
	    hessian[9]  += t2*jt[1]*jt[3];
	    hessian[10] += t2*jt[1]*jt[4];
	    hessian[11] += t2*jt[1]*jt[5];
	    
	    // hessian[12]
	    // hessian[13]
	    hessian[14] += t2*jt[2]*jt[2];
	    hessian[15] += t2*jt[2]*jt[3];
	    hessian[16] += t2*jt[2]*jt[4];
	    hessian[17] += t2*jt[2]*jt[5];
	    
	    // hessian[18]
	    // hessian[19]
	    // hessian[20]
	    hessian[21] += t2*jt[3]*jt[3];
	    hessian[22] += t2*jt[3]*jt[4];
	    hessian[23] += t2*jt[3]*jt[5];
	    
	    // hessian[24]
	    // hessian[25]
	    // hessian[26]
	    // hessian[27]
	    hessian[28] += t2*jt[4]*jt[4];
	    hessian[29] += t2*jt[4]*jt[5];

	    // hessian[30]
	    // hessian[31]
	    // hessian[32]
	    // hessian[33]
	    // hessian[34]
	    hessian[35] += t2*jt[5]*jt[5];

	    // Ignore off-diagonal terms.
	    // This approach converges incredibly slowly.
	    /*
	    hessian[0]  += t2*jt[0]*jt[0];
	    hessian[7]  += t2*jt[1]*jt[1];
	    hessian[14] += t2*jt[2]*jt[2];
	    hessian[21] += t2*jt[3]*jt[3];
	    hessian[28] += t2*jt[4]*jt[4];
	    hessian[35] += t2*jt[5]*jt[5];
	    */
	  }
	}

      }
      
      // subtract the old peak out of the foreground and background arrays.
      subtractPeak(fit_data, peak);

      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      if(info!=0){
	peak->status = ERROR;
	fit_data->n_dposv++;
	if(TESTING){
	  printf("fitting error! %d %d %d\n", i, info, ERROR);
	}
      }

      else{
	// update params
	delta[HEIGHT]     = jacobian[0];
	delta[XCENTER]    = jacobian[1];
	delta[XWIDTH]     = jacobian[2];
	delta[YCENTER]    = jacobian[3];
	delta[YWIDTH]     = jacobian[4];
	delta[BACKGROUND] = jacobian[5];

	fitDataUpdate(fit_data, peak, delta);

	// add the new peak to the foreground and background arrays.
	if (peak->status != ERROR){
	  peak->wx = calcWidth(peak->params[XWIDTH], peak->wx);
	  peak->wy = calcWidth(peak->params[YWIDTH], peak->wy);
	  addPeak(fit_data, peak);
	}
      }
    }
  }
}


/*
 * updateZ()
 *
 * Update current fits given x, y width determined by z parameter.
 *
 * This procedure is also responsible for flagging peaks
 * that might be bad & that should be removed from fitting.
 *
 * fit_data - pointer to a fitData structure.
 */
void updateZ(fitData *fit_data)
{
  // Lapack
  int n = 5, nrhs = 1, lda = 5, ldb = 5, info;

  // Local
  int i,j,k,l,m,wx,wy;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,a3,a5;
  double z0,z1,z2,zt,gx,gy;
  double delta[NPEAKPAR];
  double jt[5];
  double jacobian[5];
  double hessian[25];
  peakData *peak;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<fit_data->nfit;i++){
    peak = &fit_data->fit[i];
    if(peak->status==RUNNING){
      for(j=0;j<5;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<25;j++){
	hessian[j] = 0.0;
      }
      l = peak->offset;
      wx = peak->wx;
      wy = peak->wy;
      a1 = peak->params[HEIGHT];
      a3 = peak->params[XWIDTH];
      a5 = peak->params[YWIDTH];

      // dwx/dz calcs
      z0 = (peak->params[ZCENTER] - fit_data->wx_z_params[1]) / fit_data->wx_z_params[2];
      z1 = z0*z0;
      z2 = z1*z0;
      zt = 2.0*z0 + 3.0*fit_data->wx_z_params[3]*z1 + 4.0*fit_data->wx_z_params[4]*z2;
      gx = -2.0*zt/(fit_data->wx_z_params[0]*peak->wx_term);

      // dwy/dz calcs
      z0 = (peak->params[ZCENTER] - fit_data->wy_z_params[1]) / fit_data->wy_z_params[2];
      z1 = z0*z0;
      z2 = z1*z0;
      zt = 2.0*z0 + 3.0*fit_data->wy_z_params[3]*z1 + 4.0*fit_data->wy_z_params[4]*z2;
      gy = -2.0*zt/(fit_data->wy_z_params[0]*peak->wy_term);
      for(j=-wy;j<=wy;j++){
	yt = peak->yt[j+wy];
	eyt = peak->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j*fit_data->image_size_x + k + l;
	  fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
	  xi = fit_data->x_data[m];
	  xt = peak->xt[k+wx];
	  ext = peak->ext[k+wx];
	  e_t = ext*eyt;

	  // first derivatives
	  jt[0] = e_t;
	  jt[1] = 2.0*a1*a3*xt*e_t;
	  jt[2] = 2.0*a1*a5*yt*e_t;
	  jt[3] = -a1*xt*xt*gx*e_t-a1*yt*yt*gy*e_t;
	  jt[4] = 1.0;
	  
	  // calculate jacobian
	  t1 = 2.0*(1.0 - xi/fi);
	  jacobian[0] += t1*jt[0];
	  jacobian[1] += t1*jt[1];
	  jacobian[2] += t1*jt[2];
	  jacobian[3] += t1*jt[3];
	  jacobian[4] += t1*jt[4];
	  
	  // calculate hessian
	  t2 = 2.0*xi/(fi*fi);

	  // calculate hessian without second derivative terms.
	  hessian[0] += t2*jt[0]*jt[0];
	  hessian[1] += t2*jt[0]*jt[1];
	  hessian[2] += t2*jt[0]*jt[2];
	  hessian[3] += t2*jt[0]*jt[3];
	  hessian[4] += t2*jt[0]*jt[4];
	  
	  // hessian[5]
	  hessian[6] += t2*jt[1]*jt[1];
	  hessian[7] += t2*jt[1]*jt[2];
	  hessian[8] += t2*jt[1]*jt[3];
	  hessian[9] += t2*jt[1]*jt[4];
	    
	  // hessian[10]
	  // hessian[11]
	  hessian[12] += t2*jt[2]*jt[2];
	  hessian[13] += t2*jt[2]*jt[3];
	  hessian[14] += t2*jt[2]*jt[4];
	  
	  // hessian[15]
	  // hessian[16]
	  // hessian[17]
	  hessian[18] += t2*jt[3]*jt[3];
	  hessian[19] += t2*jt[3]*jt[4];

	  // hessian[20]
	  // hessian[21]
	  // hessian[22]
	  // hessian[23]
	  hessian[24] += t2*jt[4]*jt[4];
	}
      }

      // subtract the old peak out of the foreground and background arrays.
      subtractPeak(fit_data, peak);
      
      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      if(info!=0){
	peak->status = ERROR;
	fit_data->n_dposv++;
	if(TESTING){
	  printf("fitting error! %d %d %d\n", i, info, ERROR);
	}
      }
      else{
	// update params
	delta[HEIGHT]     = jacobian[0];
	delta[XCENTER]    = jacobian[1];
	delta[YCENTER]    = jacobian[2];
	delta[ZCENTER]    = jacobian[3];
	delta[BACKGROUND] = jacobian[4];

	fitDataUpdate(fit_data, peak, delta);

	// add the new peak to the foreground and background arrays.
	if (peak->status != ERROR){
	  // calculate new x,y width, update fit area.
	  calcWidthsFromZ(fit_data, peak);
	  peak->wx = calcWidth(peak->params[XWIDTH], peak->wx);
	  peak->wy = calcWidth(peak->params[YWIDTH], peak->wy);
	  addPeak(fit_data, peak);
	}
      }
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
