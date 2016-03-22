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
#define MARGIN 10      /* Margin around the edge of the image. */

#define MINZ -0.5
#define MAXZ 0.5

/* Structures */
typedef struct
{
  int sign[NPEAKPAR];
  int status;
  int offset;
  int wx;
  int wy;
  int xc;
  int yc;
  double clamp[NPEAKPAR];
  double error;
  double error_old;
  double params[NPEAKPAR];  /* [height x-center x-width y-center y-width background] */
  double wx_term;
  double wy_term;
  double xt[2*MARGIN+1];
  double ext[2*MARGIN+1];
  double yt[2*MARGIN+1];
  double eyt[2*MARGIN+1];
} fitData;


/* Function Declarations */
void addPeak(fitData *);
void calcErr();
void calcFit();
int calcWidth(double, int);
void calcWidthsFromZ(fitData *);
void cleanup();
void fitDataUpdate(fitData *, double *);
double getError();
void getResidual(double *);
void getResults(double *);
int getUnconverged();
void initialize(double *, double *, double *, double, int, int, int, int);
void initializeZParameters(double *, double *, double, double);
void iterate2DFixed();
void iterate2D();
void iterate3D();
void iterateZ();
void subtractPeak(fitData *);
void update2DFixed();
void update2D();
void update3D();
void updateZ();

/* LAPACK Functions */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);


/* Global Variables */
static int nfit;                 /* number of peaks to fit. */
static int image_size_x;         /* size in x (fast axis). */
static int image_size_y;         /* size in y (slow axis). */
static int *bg_counts;           /* number of peaks covering a particular pixel. */
static double tolerance;         /* fit tolerance. */
static double min_z = MINZ;      /* minimum z value. */
static double max_z = MAXZ;      /* maximum z value. */
static double *bg_data;          /* background data. */
static double *f_data;           /* fit (foreground) data. */
static double *x_data;           /* image data. */
static double *scmos_term;       /* sCMOS calibration term for each pixel (var/gain^2). */
static double wx_z_params[5];    /* x width versus z parameters. */
static double wy_z_params[5];    /* y width versus z parameters. */

static fitData *fit;


/* Functions */


/*
 * addPeak()
 *
 * Add a peak to the foreground and background data arrays.
 */
void addPeak(fitData *cur)
{
  int j,k,l,m,n,wx,wy,xc,yc;
  double bg,mag,tmp,xt,yt;

  xc = cur->xc;
  yc = cur->yc;
  cur->offset = yc*image_size_x+xc;

  wx = cur->wx;
  wy = cur->wy;

  for(j=(xc-wx);j<=(xc+wx);j++){
    xt = (double)j - cur->params[XCENTER];
    n = j-xc+wx;
    cur->xt[n] = xt;
    cur->ext[n] = exp(-xt*xt*cur->params[XWIDTH]);
  }
  for(j=(yc-wy);j<=(yc+wy);j++){
    yt = (double)j - cur->params[YCENTER];
    n = j-yc+wy;
    cur->yt[n] = yt;
    cur->eyt[n] = exp(-yt*yt*cur->params[YWIDTH]);
  }

  /* gaussian function */
  l = cur->offset;
  bg = cur->params[BACKGROUND];
  mag = cur->params[HEIGHT];
  for(j=-wy;j<=wy;j++){
    tmp = cur->eyt[j+wy];
    for(k=-wx;k<=wx;k++){
      m = j*image_size_x+k+l;
      f_data[m] += mag*tmp*cur->ext[k+wx];
      bg_counts[m] += 1;
      bg_data[m] += bg + scmos_term[m];
    }
  }

}


/*
 * calcErr()
 *
 * Calculate error in fit for each gaussian.
 */
void calcErr()
{
  int i,j,k,l,m,wx,wy;
  double err,fi,xi;

  for(i=0;i<nfit;i++){
    if(fit[i].status == RUNNING){
      l = fit[i].offset;
      wx = fit[i].wx;
      wy = fit[i].wy;
      err = 0.0;
      for(j=-wy;j<=wy;j++){
	for(k=-wx;k<=wx;k++){
	  m = (j*image_size_x)+k+l;
	  fi = f_data[m]+bg_data[m]/((double)bg_counts[m]);
	  if (fi <= 0.0){
	    if(VERBOSE){
	      printf(" Negative f detected! %.3f %.3f %.3f %.3f %d\n", fit[i].params[BACKGROUND], fi, f_data[m], bg_data[m], bg_counts[m]);
	    }
	    fit[i].status = ERROR;
	    j = wy + 1;
	    k = wx + 1;
	  }
	  xi = x_data[m];
	  if (xi <= 0.0){
	    if(VERBOSE){
	      printf(" Negative x detected! %.3f\n", x_data[m]);
	    }
	  }
	  err += 2*(fi-xi)-2*xi*log(fi/xi);
	}
      }
      fit[i].error_old = fit[i].error;
      fit[i].error = err;
      if (VERBOSE){
	printf("%d %f %f %f\n", i, fit[i].error_old, fit[i].error, tolerance);
      }
      if(((fabs(err - fit[i].error_old)/err) < tolerance)&&(fit[i].status!=ERROR)){
	fit[i].status = CONVERGED;
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
void calcFit()
{
  int i;
  fitData *cur;

  // zero matrices.
  for(i=0;i<(image_size_x*image_size_y);i++){
    f_data[i] = 1.0;
    bg_counts[i] = 0;
    bg_data[i] = 0;
  }

  // update fit matrix with values from fits.
  for(i=0;i<nfit;i++){
    cur = &fit[i];
    if (cur->status != ERROR){
      addPeak(cur);
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
 * cur - fit data structure to update.
 */
void calcWidthsFromZ(fitData *cur)
{
  double z0,z1,z2,z3,tmp;

  // wx
  z0 = (cur->params[ZCENTER]-wx_z_params[1])/wx_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  z3 = z2*z0;
  tmp = 1.0 + z1 + wx_z_params[3]*z2 + wx_z_params[4]*z3;
  cur->wx_term = tmp*tmp;
  cur->params[XWIDTH] = 2.0/(wx_z_params[0]*tmp);

  // wy
  z0 = (cur->params[ZCENTER]-wy_z_params[1])/wy_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  z3 = z2*z0;
  tmp = 1.0 + z1 + wy_z_params[3]*z2 + wy_z_params[4]*z3;
  cur->wy_term = tmp*tmp;
  cur->params[YWIDTH] = 2.0/(wy_z_params[0]*tmp);

  // printf("%.4f %.4f %.4f\n",cur->params[ZCENTER],cur->params[XWIDTH],cur->params[YWIDTH]);
}


/*
 * cleanup()
 *
 * Frees storage associated gaussian fitting.
 */
void cleanup()
{
  free(fit);
  free(x_data);
  free(f_data);
  free(scmos_term);
  free(bg_data);
  free(bg_counts);
}


/*
 * fitDataUpdate()
 *
 * Updates fit data given deltas.
 *
 * Also checks for out-of-bounds parameters.
 */
void fitDataUpdate(fitData *cur, double *delta)
{
  int i,xc,yc;

  // update
  for(i=0;i<NPEAKPAR;i++){
    if(VERBOSE){
      printf("%.3e %.3f | ", delta[i], cur->clamp[i]);
    }

    // update sign & clamp if the solution appears to be oscillating.
    if (cur->sign[i] != 0){
      if ((cur->sign[i] == 1) && (delta[i] < 0.0)){
	cur->clamp[i] *= 0.5;
      }
      else if ((cur->sign[i] == -1) && (delta[i] > 0.0)){
	cur->clamp[i] *= 0.5;
      }
    }
    if (delta[i] > 0.0){
      cur->sign[i] = 1;
    }
    else {
      cur->sign[i] = -1;
    }

    // update values based on delta & clamp.
    if (delta[i] != 0.0){
      cur->params[i] -= delta[i]/(1.0 + fabs(delta[i])/cur->clamp[i]);
    }
  }
  if(VERBOSE){
    printf("\n");
  }

  // Update peak (integer) center (w/ hysteresis).
  if(fabs(cur->params[XCENTER] - (double)cur->xc - 0.5) > HYSTERESIS){
    cur->xc = (int)cur->params[XCENTER];
  }
  if(fabs(cur->params[YCENTER] - (double)cur->yc - 0.5) > HYSTERESIS){
    cur->yc = (int)cur->params[YCENTER];
  }

  // Check that the peak hasn't moved to close to the 
  // edge of the image. Flag the peak as bad if it has.
  xc = cur->xc;
  yc = cur->yc;
  if((xc<=MARGIN)||(xc>=(image_size_x-MARGIN-1))||(yc<=MARGIN)||(yc>=(image_size_y-MARGIN-1))){
    cur->status = ERROR;
    if(TESTING){
      printf("object outside margins, %.3f, %.3f\n", cur->params[XCENTER], cur->params[YCENTER]);
    }
  }
  
  // check for negative background or height
  //if((cur->params[BACKGROUND]<0.0)||(cur->params[HEIGHT]<0.0)){
  if(cur->params[HEIGHT]<0.0){
    cur->status = ERROR;
    if(TESTING){
      printf("negative height, %.3f, %.3f (%.3f, %.3f)\n", cur->params[BACKGROUND], cur->params[HEIGHT], cur->params[XCENTER], cur->params[YCENTER]);
    }
  }

  // check for negative widths
  if((cur->params[XWIDTH]<0.0)||(cur->params[YWIDTH]<0.0)){
    cur->status = ERROR;
    if(TESTING){
      printf("negative widths, %.3f, %.3f (%.3f, %.3f)\n", cur->params[XWIDTH], cur->params[YWIDTH], cur->params[XCENTER], cur->params[YCENTER]);
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
  if(cur->params[ZCENTER]<min_z){
    cur->params[ZCENTER] = min_z;
  }

  if(cur->params[ZCENTER]>max_z){
    cur->params[ZCENTER] = max_z;
  }

}


/*
 * getError()
 *
 * Return the current error in the fit.
 */
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

  // printf("%f\n",err);
  return err;
}


/*
 * getResidual(residual).
 *
 * Returns image - fit.
 *
 * residual - Pre-allocated space to store the residual values.
 *            This should be square & the same size as the image.
 */
void getResidual(double *residual)
{
  int i;

  calcFit();
  for(i=0;i<(image_size_x*image_size_y);i++){
    residual[i] = x_data[i] - f_data[i];
    //residual[i] = x_data[i] - (f_data[i] + scmos_term[i]);
  }
}


/*
 * getResults(params)
 *
 * Return the current fitting results.
 *
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
void getResults(double *peak_params)
{
  int i;
  fitData *cur;

  for(i=0;i<nfit;i++){
    cur = &fit[i];

    if(cur->status != ERROR){
      peak_params[i*NRESULTSPAR+XWIDTH] = sqrt(1.0/(2.0*cur->params[XWIDTH]));
      peak_params[i*NRESULTSPAR+YWIDTH] = sqrt(1.0/(2.0*cur->params[YWIDTH]));
    }
    else{
      peak_params[i*NRESULTSPAR+XWIDTH] = 1.0;
      peak_params[i*NRESULTSPAR+YWIDTH] = 1.0;
    }
    peak_params[i*NRESULTSPAR+HEIGHT]     = cur->params[HEIGHT];
    peak_params[i*NRESULTSPAR+XCENTER]    = cur->params[XCENTER];
    peak_params[i*NRESULTSPAR+YCENTER]    = cur->params[YCENTER];
    peak_params[i*NRESULTSPAR+BACKGROUND] = cur->params[BACKGROUND];
    peak_params[i*NRESULTSPAR+ZCENTER]    = cur->params[ZCENTER];

    peak_params[i*NRESULTSPAR+STATUS] = (double)cur->status;
    peak_params[i*NRESULTSPAR+IERROR] = cur->error;
  }
}


/*
 * getUnconverged()
 *
 * Return the number of fits that have not yet converged.
 */
int getUnconverged()
{
  int i,count;

  count = 0;
  for(i=0;i<nfit;i++){
    if(fit[i].status==RUNNING){
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
 * image - pointer to the image data.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * params - pointer to the initial values for the parameters for each point.
 *           1. height
 *           2. x-center
 *           3. x-sigma
 *           4. y-center
 *           5. y-sigma
 *           6. background
 *           7. z-center
 *           8. status
 *           9. error
 *               .. repeat ..
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 * n - number of parameters.
 * zfit - fitting wx, wy based on z?
 */
void initialize(double *image, double *scmos_calibration, double *params, double tol, int im_size_x, int im_size_y, int n, int zfit)
{
  int i,j;

  nfit = n;
  image_size_x = im_size_x;
  image_size_y = im_size_y;
  tolerance = tol;

  x_data = (double *)malloc(sizeof(double)*image_size_x*image_size_y);
  f_data = (double *)malloc(sizeof(double)*image_size_x*image_size_y);
  bg_data = (double *)malloc(sizeof(double)*image_size_x*image_size_y);
  bg_counts = (int *)malloc(sizeof(int)*image_size_x*image_size_y);
  scmos_term = (double *)malloc(sizeof(double)*image_size_x*image_size_y);

  for(i=0;i<(image_size_x*image_size_y);i++){
    x_data[i] = image[i];
    scmos_term[i] = scmos_calibration[i];
  }

  fit = (fitData *)malloc(sizeof(fitData)*nfit);
  for(i=0;i<nfit;i++){
    fit[i].status = (int)(params[i*NRESULTSPAR+STATUS]);
    if(fit[i].status==RUNNING){
      fit[i].error = 0.0;
      fit[i].error_old = 0.0;
    }
    else {
      fit[i].error = params[i*NRESULTSPAR+IERROR];
      fit[i].error_old = fit[i].error;
    }
    
    fit[i].params[HEIGHT]     = params[i*NRESULTSPAR+HEIGHT];
    fit[i].params[XCENTER]    = params[i*NRESULTSPAR+XCENTER];
    fit[i].params[YCENTER]    = params[i*NRESULTSPAR+YCENTER];
    fit[i].params[BACKGROUND] = params[i*NRESULTSPAR+BACKGROUND];
    fit[i].params[ZCENTER]    = params[i*NRESULTSPAR+ZCENTER];

    if(zfit){
      calcWidthsFromZ(&fit[i]);
    }
    else{
      fit[i].params[XWIDTH] = 1.0/(2.0*params[i*NRESULTSPAR+XWIDTH]*params[i*NRESULTSPAR+XWIDTH]);
      fit[i].params[YWIDTH] = 1.0/(2.0*params[i*NRESULTSPAR+YWIDTH]*params[i*NRESULTSPAR+YWIDTH]);
    }

    // printf("%d %.3f %.3f %.2f %.2f\n", i, fit[i].params[XCENTER], fit[i].params[YCENTER], fit[i].params[HEIGHT], fit[i].params[BACKGROUND]);
    fit[i].xc = (int)fit[i].params[XCENTER];
    fit[i].yc = (int)fit[i].params[YCENTER];
    fit[i].wx = calcWidth(fit[i].params[XWIDTH],-10.0);
    fit[i].wy = calcWidth(fit[i].params[YWIDTH],-10.0);

    //
    // FIXME: It would probably better to pass these in.
    //        Currently "tuned" for "standard" STORM images.
    //
    fit[i].clamp[HEIGHT] = 1000.0;
    fit[i].clamp[XCENTER] = 1.0;
    fit[i].clamp[XWIDTH] = 0.3;
    fit[i].clamp[YCENTER] = 1.0;
    fit[i].clamp[YWIDTH] = 0.3;
    fit[i].clamp[BACKGROUND] = 100.0;
    fit[i].clamp[ZCENTER] = 0.1;

    for(j=0;j<NPEAKPAR;j++){
      fit[i].sign[j] = 0;
    }
  }

  calcFit();
  calcErr();
}


/*
 * initializeZParameters(wx_vs_z, wy_vs_z)
 *
 * Initializes fitting for z with wx, wy dependence on z.
 *
 * wx_vs_z - [wo, c, d, A, B] wx vs z coefficients.
 * wy_vs_z - [wo, c, d, A, B] wy vs z coefficients.
 * z_min - minimum allowed z value.
 * z_max - maximum allowed z value.
 */
void initializeZParameters(double *wx_vs_z, double *wy_vs_z, double z_min, double z_max)
{
  int i;

  for(i=0;i<5;i++){
    wx_z_params[i] = wx_vs_z[i];
    wy_z_params[i] = wy_vs_z[i];
    //printf("iZ: %d %f %f\n", i, wx_vs_z[i], wy_vs_z[i]);
  }
  wx_z_params[0] = wx_z_params[0]*wx_z_params[0];
  wy_z_params[0] = wy_z_params[0]*wy_z_params[0];

  min_z = z_min;
  max_z = z_max;
}


/*
 * iterate2DFixed()
 *
 * Performs a single cycle of fit improvement with fixed x, y width.
 *
 */
void iterate2DFixed()
{
  update2DFixed();
  calcErr();
}


/*
 * iterate2D()
 *
 * Performs a single cycle of fit improvement with equal x, y width.
 *
 */
void iterate2D()
{
  update2D();
  calcErr();
}


/*
 * iterate3D()
 *
 * Performs a single cycle of fit improvement.
 */
void iterate3D()
{
  update3D();
  calcErr();
}


/*
 * iterateZ()
 *
 * Performs a single cycle of fit improvement with x, y width
 * determined by the z parameter.
 */
void iterateZ()
{
  updateZ();
  calcErr();
}


/*
 * subtractPeak()
 *
 * Subtract the peak out of the current fit, basically 
 * this just undoes addPeak().
 */
void subtractPeak(fitData *cur)
{
  int j,k,l,m,wx,wy;
  double bg,mag,tmp;

  wx = cur->wx;
  wy = cur->wy;

  /* gaussian function */
  l = cur->offset;
  bg = cur->params[BACKGROUND];
  mag = cur->params[HEIGHT];
  for(j=-wy;j<=wy;j++){
    tmp = cur->eyt[j+wy];
    for(k=-wx;k<=wx;k++){
      m = j*image_size_x+k+l;
      f_data[m] -= mag*tmp*cur->ext[k+wx];
      bg_counts[m] -= 1;
      bg_data[m] -= (bg + scmos_term[m]);
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
 */
void update2DFixed()
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
  fitData *cur;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<nfit;i++){
    cur = &fit[i];
    if(cur->status==RUNNING){
      for(j=0;j<4;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<16;j++){
	hessian[j] = 0.0;
      }
      l = cur->offset;
      wx = cur->wx;
      wy = cur->wy;
      a1 = cur->params[HEIGHT];
      width = cur->params[XWIDTH];
      for(j=-wy;j<=wy;j++){
	yt = cur->yt[j+wy];
	eyt = cur->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j*image_size_x+k+l;
	  fi = f_data[m]+bg_data[m]/((double)bg_counts[m]);
	  xi = x_data[m];
	  xt = cur->xt[k+wx];
	  ext = cur->ext[k+wx];
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
      subtractPeak(cur);
      
      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      if(info!=0){
	cur->status = ERROR;
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
	fitDataUpdate(cur, delta);

	// add the new peak to the foreground and background arrays.
	if (cur->status != ERROR){
	  addPeak(cur);
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
 */
void update2D()
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
  fitData *cur;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<nfit;i++){
    cur = &fit[i];
    if(cur->status==RUNNING){
      for(j=0;j<5;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<25;j++){
	hessian[j] = 0.0;
      }
      l = cur->offset;
      wx = cur->wx;
      wy = cur->wy;
      a1 = cur->params[HEIGHT];
      width = cur->params[XWIDTH];
      for(j=-wy;j<=wy;j++){
	yt = cur->yt[j+wy];
	eyt = cur->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j*image_size_x+k+l;
	  fi = f_data[m]+bg_data[m]/((double)bg_counts[m]);
	  xi = x_data[m];
	  xt = cur->xt[k+wx];
	  ext = cur->ext[k+wx];
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
      subtractPeak(cur);

      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      if(info!=0){
	cur->status = ERROR;
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

	fitDataUpdate(cur, delta);

	// add the new peak to the foreground and background arrays.
	// recalculate peak fit area as the peak width may have changed.
	if (cur->status != ERROR){
	  cur->wx = calcWidth(cur->params[XWIDTH],cur->wx);
	  cur->wy = cur->wx;
	  addPeak(cur);
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
 */
void update3D()
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
  fitData *cur;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<nfit;i++){
    cur = &fit[i];
    if(cur->status==RUNNING){
      for(j=0;j<6;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<36;j++){
	hessian[j] = 0.0;
      }
      l = cur->offset;
      wx = cur->wx;
      wy = cur->wy;
      a1 = cur->params[HEIGHT];
      a3 = cur->params[XWIDTH];
      a5 = cur->params[YWIDTH];
      for(j=-wy;j<=wy;j++){
	yt = cur->yt[j+wy];
	eyt = cur->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j*image_size_x+k+l;
	  fi = f_data[m]+bg_data[m]/((double)bg_counts[m]);
	  xi = x_data[m];
	  xt = cur->xt[k+wx];
	  ext = cur->ext[k+wx];
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
      
      /*
      printf("hessian:\n");
      for(j=0;j<6;j++){
	for(k=0;k<6;k++){
	  printf("%.4f ", hessian[j*6+k]);
	}
	printf("\n");
      }
      printf("\n");
      */
      // subtract the old peak out of the foreground and background arrays.
      subtractPeak(cur);

      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      if(info!=0){
	cur->status = ERROR;
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

	fitDataUpdate(cur, delta);

	// add the new peak to the foreground and background arrays.
	if (cur->status != ERROR){
	  cur->wx = calcWidth(cur->params[XWIDTH],cur->wx);
	  cur->wy = calcWidth(cur->params[YWIDTH],cur->wy);
	  addPeak(cur);
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
 */
void updateZ()
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
  fitData *cur;

  for(i=0;i<NPEAKPAR;i++){
    delta[i] = 0.0;
  }

  for(i=0;i<nfit;i++){
    cur = &fit[i];
    if(cur->status==RUNNING){
      for(j=0;j<5;j++){
	jacobian[j] = 0.0;
      }
      for(j=0;j<25;j++){
	hessian[j] = 0.0;
      }
      l = cur->offset;
      wx = cur->wx;
      wy = cur->wy;
      a1 = cur->params[HEIGHT];
      a3 = cur->params[XWIDTH];
      a5 = cur->params[YWIDTH];

      // dwx/dz calcs
      z0 = (cur->params[ZCENTER]-wx_z_params[1])/wx_z_params[2];
      z1 = z0*z0;
      z2 = z1*z0;
      zt = 2.0*z0+3.0*wx_z_params[3]*z1+4.0*wx_z_params[4]*z2;
      gx = -2.0*zt/(wx_z_params[0]*cur->wx_term);

      // dwy/dz calcs
      z0 = (cur->params[ZCENTER]-wy_z_params[1])/wy_z_params[2];
      z1 = z0*z0;
      z2 = z1*z0;
      zt = 2.0*z0+3.0*wy_z_params[3]*z1+4.0*wy_z_params[4]*z2;
      gy = -2.0*zt/(wy_z_params[0]*cur->wy_term);
      for(j=-wy;j<=wy;j++){
	yt = cur->yt[j+wy];
	eyt = cur->eyt[j+wy];
	for(k=-wx;k<=wx;k++){
	  m = j*image_size_x+k+l;
	  fi = f_data[m]+bg_data[m]/((double)bg_counts[m]);
	  xi = x_data[m];
	  xt = cur->xt[k+wx];
	  ext = cur->ext[k+wx];
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
      subtractPeak(cur);
      
      // Use Lapack to solve AX=B to calculate update vector
      dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );

      /*
      for(j=0;j<5;j++){
	for(k=0;k<5;k++){
	  printf("%.4f ", hessian[j*5+k]);
	}
	printf("\n");
      }
      */

      if(info!=0){
	cur->status = ERROR;
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

	fitDataUpdate(cur, delta);

	// add the new peak to the foreground and background arrays.
	if (cur->status != ERROR){
	  // calculate new x,y width, update fit area.
	  calcWidthsFromZ(cur);
	  cur->wx = calcWidth(cur->params[XWIDTH],cur->wx);
	  cur->wy = calcWidth(cur->params[YWIDTH],cur->wy);
	  addPeak(cur);
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
