/*
 * Fit multiple, possible overlapping, cubic splines simultaneously to
 * image data.
 *
 * Hazen 01/14
 *
 *
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall cubic_fit.c
 *  gcc -shared -Wl,-soname,cubic_fit.so.1 -o cubic_fit.so.1.0.1 cubic_fit.o -lc -llapack
 *
 * Windows:
 *  gcc -c cubic_fit.c
 *  gcc -shared -o cubic_fit.dll cubic_fit.o cubic_spline.o multi_fit.o -llapack -Lc:\Users\Hazen\lib
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cubic_spline.h"
#include "multi_fit_core.h"

/* Define */
#define DEBUG 0
#define TESTING 1

#define F2D 0
#define F3D 1

#define CF_HEIGHT 0
#define CF_XCENTER 1
#define CF_YCENTER 2
#define CF_BACKGROUND 3
#define CF_ZCENTER 4

/* Functions */
void freePeaks(void);
void getResults(double *);
int getUnconverged(void);
int getZOff(void);
void iterateSpline(void);
void newPeaks(double *, int, int);
void newPeaks2D(double *, int);
void newPeaks3D(double *, int);
void updateFitValues2D(fitData *);
void updateFitValues3D(fitData *);
void updatePeakParameters(fitData *);

/* Global variables */
static int n_fit_data;
static double xoff;
static double yoff;
static double zoff;
static fitData *new_fit_data;
static fitData *old_fit_data;

/*
 * freePeaks()
 *
 * Free fit_data and old_fit_data.
 */
void freePeaks(void)
{
  int i;

  if (new_fit_data != NULL){
    for(i=0;i<n_fit_data;i++){
      freeFitData(&(new_fit_data[i]));
      freeFitData(&(old_fit_data[i]));
    }
    free(new_fit_data);
    free(old_fit_data);
  }
}

/*
 * getResults()
 *
 * Return the current fitting results.
 *
 * peaks - pre-allocated space for storing the peak fitting parameters.
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
 */
void getResults(double *peaks)
{
  int i;
  double dx,dy,dz;
  fitData *a_peak;

  for(i=0;i<n_fit_data;i++){
    a_peak = &(new_fit_data[i]);

    if (DEBUG){
      printf(" gR: %.2f %.2f\n", a_peak->params[CF_XCENTER] + xoff, a_peak->params[CF_YCENTER] + yoff);
    }

    peaks[i*NRESULTSPAR+HEIGHT] = a_peak->params[CF_HEIGHT];
    peaks[i*NRESULTSPAR+XCENTER] = a_peak->params[CF_XCENTER] + xoff;
    peaks[i*NRESULTSPAR+YCENTER] = a_peak->params[CF_YCENTER] + yoff;
    if (a_peak->type == F3D){
      peaks[i*NRESULTSPAR+ZCENTER] = a_peak->params[CF_ZCENTER] + zoff;
    }
    peaks[i*NRESULTSPAR+BACKGROUND] = a_peak->params[CF_BACKGROUND];
    peaks[i*NRESULTSPAR+STATUS] = (double)a_peak->status;
    peaks[i*NRESULTSPAR+IERROR] = a_peak->error;

    /* Dummy values for peak width and height. */
    peaks[i*NRESULTSPAR+XWIDTH] = 1.0;
    peaks[i*NRESULTSPAR+YWIDTH] = 1.0;
  }
}

/*
 * getUnconverged()
 *
 * Return the number of fits that have not yet converged.
 */
int getUnconverged(void)
{
  int i,count;

  count = 0;
  for(i=0;i<n_fit_data;i++){
    if(new_fit_data[i].status==RUNNING){
      count++;
    }
  }

  return count;
}

/*
 * getZOff()
 *
 * Return what is being used for the z offset.
 */
int getZOff(void){
  return zoff;
}

/*
 * iterateSpline()
 *
 * Perform one cycle of update of all the peaks.
 */
void iterateSpline(void)
{
  int i;
  fitData *new_peak, *old_peak;
  
  for(i=0;i<n_fit_data;i++){
    new_peak = &(new_fit_data[i]);
    old_peak = &(old_fit_data[i]);
    if (new_peak->status == RUNNING){

      /* Add the peak & calculate the error. */
      new_peak->error = calcError(new_peak, new_peak->params[CF_BACKGROUND]);
      //printf("error: %f %f %f\n", new_peak->error, old_peak->error, new_peak->lambda);

      if((fabs(new_peak->error - old_peak->error)/new_peak->error) < tolerance){
	new_peak->status = CONVERGED;
      }
      else{
	/* 
	 * Check if the fit is getting worse. If it is then reset to the 
	 * previous parameters & increase the lambda parameter. This will 
	 * cause the search to move downhill with higher probability but 
	 * less speed. 
	 */
	if(new_peak->error > (1.5*old_peak->error)){
	  subtractPeak(new_peak);
	  copyFitData(old_peak, new_peak);
	  new_peak->lambda = 10.0 * new_peak->lambda;
	  addPeak(new_peak);
	}
	/*
	 * If the fit is getting better and lambda is greater than 1.0
	 * then we slowly decrease lambda.
	 */
	else if(new_peak->error < old_peak->error){
	  if(new_peak->lambda > 1.0){
	    new_peak->lambda = 0.8*new_peak->lambda;
	    if(new_peak->lambda < 1.0){
	      new_peak->lambda = 1.0;
	    }
	  }
	}

	/* Calculate updated fit values. */
	updateFit(new_peak, new_peak->params[CF_BACKGROUND]);

	/* 
	 * Increase lambda in the event of Cholesky failure. 
	 */
	if(new_peak->status == CHOLERROR){
	  new_peak->lambda = 10.0 * new_peak->lambda;
	  new_peak->status = RUNNING;
	}
	else{
	  subtractPeak(new_peak);
	  copyFitData(new_peak, old_peak);
	  updatePeakParameters(new_peak);
	  if(new_peak->type == F2D){
	    updateFitValues2D(new_peak);
	  }
	  else{
	    updateFitValues3D(new_peak);
	  }
	  if(new_peak->status != BADPEAK){
	    addPeak(new_peak);
	  }
	}
      }
    }
  }
}

/*
 * newPeaks()
 *
 * Create initial fit data structures for spline fitting from the peak array. The 
 * format of the initial values is the same as what was used for 3D-DAOSTORM.
 *
 * peaks - pointer to the initial values for the parameters for each peak.
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
 * n_peaks - The number of parameters (peaks).
 */
void newPeaks(double *peaks, int n_peaks, int fit_type)
{
  int i,j,n_fit_params,sx,sy,sz;
  double temp;

  if(fit_type == F2D){
    n_fit_params = 4;
  }
  else{
    n_fit_params = 5;
  }

  /* Delete old fit data structures, if they exist. */
  freePeaks();
  
  /* Reset the fit array. */
  resetFit();

  /* Allocate storage for fit data structure. */
  n_fit_data = n_peaks;
  new_fit_data = (fitData *)malloc(sizeof(fitData)*n_peaks);
  old_fit_data = (fitData *)malloc(sizeof(fitData)*n_peaks);
  sx = getXSize()/2 - 1;
  sy = getYSize()/2 - 1;
  sz = getZSize();

  /* Note the assumption here that the splines are square. */
  if((sx%2)==1){
    xoff = (double)(sx/2) - 1.0;
    yoff = (double)(sy/2) - 1.0;
  }
  else{
    xoff = (double)(sx/2) - 1.5;
    yoff = (double)(sy/2) - 1.5;
  }

  zoff = -(double)(sz/2);

  /* Initialize fit data structure. */
  for(i=0;i<n_fit_data;i++){
    new_fit_data[i].index = i;
    new_fit_data[i].n_params = n_fit_params;
    new_fit_data[i].size_x = sx;
    new_fit_data[i].size_y = sy;
    new_fit_data[i].size_z = sz;
    new_fit_data[i].status = (int)(peaks[i*NRESULTSPAR+STATUS]);
    new_fit_data[i].type = fit_type;

    new_fit_data[i].lambda = 1.0;

    mallocFitData(&(new_fit_data[i]), n_fit_params, sx*sy*(n_fit_params+1));

    for(j=0;j<n_fit_params;j++){
      new_fit_data[i].sign[j] = 0;
    }

    new_fit_data[i].clamp[CF_HEIGHT] = 100.0;
    new_fit_data[i].clamp[CF_XCENTER] = 1.0;
    new_fit_data[i].clamp[CF_YCENTER] = 1.0;
    new_fit_data[i].clamp[CF_BACKGROUND] = 10.0;
    if (fit_type == F3D){
      new_fit_data[i].clamp[CF_ZCENTER] = 2.0;
    }

    if (DEBUG){
      printf(" peaks: %.2f %.2f\n", peaks[i*NRESULTSPAR+XCENTER], peaks[i*NRESULTSPAR+YCENTER]);
    }
    new_fit_data[i].params[CF_HEIGHT] = peaks[i*NRESULTSPAR+HEIGHT];
    new_fit_data[i].params[CF_XCENTER] = peaks[i*NRESULTSPAR+XCENTER] - xoff;
    new_fit_data[i].params[CF_YCENTER] = peaks[i*NRESULTSPAR+YCENTER] - yoff;
    new_fit_data[i].params[CF_BACKGROUND] = peaks[i*NRESULTSPAR+BACKGROUND];
    if (fit_type == F3D){
      new_fit_data[i].params[CF_ZCENTER] = peaks[i*NRESULTSPAR+ZCENTER] - zoff;
    }

    new_fit_data[i].xi = (int)(new_fit_data[i].params[CF_XCENTER]);
    new_fit_data[i].yi = (int)(new_fit_data[i].params[CF_YCENTER]);
    if (fit_type == F3D){
      new_fit_data[i].zi = (int)(new_fit_data[i].params[CF_ZCENTER]);
    }

    if (fit_type == F2D){
      updateFitValues2D(&(new_fit_data[i]));
    }
    else{
      updateFitValues3D(&(new_fit_data[i]));
    }

    /*
     * The derivative with respect to background term is always 1.0
     */
    for(j=0;j<(sx*sy);j++){
      new_fit_data[i].values[(CF_BACKGROUND+1)*sx*sy+j] = 1.0;
    }

    addPeak(&(new_fit_data[i]));

    mallocFitData(&(old_fit_data[i]), n_fit_params, sx*sy*(n_fit_params+1));
    copyFitData(&(new_fit_data[i]), &(old_fit_data[i]));

    if(new_fit_data[i].status == RUNNING){
      if (TESTING){
	temp = calcError(&(new_fit_data[i]), new_fit_data[i].params[CF_BACKGROUND]);
	new_fit_data[i].error = temp;
      }
      else{
	new_fit_data[i].error = 1.0e+12;
      }
      old_fit_data[i].error = 1.0e+13;
    }
    else{
      new_fit_data[i].error = peaks[i*NRESULTSPAR+IERROR];
      old_fit_data[i].error = new_fit_data[i].error;
    }
  }
}

/*
 * newPeaks2D()
 *
 * A thin wrapper of newPeaks().
 *
 * peaks - pointer to the initial values for the parameters for each peak.
 * n_peaks - The number of parameters (peaks).
 */
void newPeaks2D(double *peaks, int n_peaks)
{
  newPeaks(peaks, n_peaks, F2D);
}

/*
 * newPeaks3D()
 *
 * A thin wrapper of newPeaks().
 *
 * peaks - pointer to the initial values for the parameters for each peak.
 * n_peaks - The number of parameters (peaks).
 */
void newPeaks3D(double *peaks, int n_peaks)
{
  newPeaks(peaks, n_peaks, F3D);
}

/*
 * updateFitValues2D()
 *
 * Use the cubic_spline library to update the fit values for a peak
 * based on the current parameters.
 *
 * a_peak - The peak to update.
 */
void updateFitValues2D(fitData *a_peak)
{
  int i,j,l,psx,psy,size,xstart,ystart;
  double height,temp,xc,yc;

  psx = a_peak->size_x;
  psy = a_peak->size_y;
  size = psx*psy;

  /* Calculate new fit values. */
  height = a_peak->params[CF_HEIGHT];
  xc = 2.0*(2.0 - (a_peak->params[CF_XCENTER] - a_peak->xi));
  yc = 2.0*(2.0 - (a_peak->params[CF_YCENTER] - a_peak->yi));

  xstart = 0;
  while(xc>1.0){
    xstart += 1;
    xc -= 1.0;
  }

  ystart = 0;
  while(yc>1.0){
    ystart += 1;
    yc -= 1.0;
  }

  computeDelta2D(yc, xc);

  for(i=0;i<psy;i++){
    for(j=0;j<psx;j++){
      l = i*psx+j;
      temp = fAt2D(2*i+ystart,2*j+xstart);
      a_peak->values[l] = height*temp;
      a_peak->values[(CF_HEIGHT+1)*size+l] = temp;
      a_peak->values[(CF_XCENTER+1)*size+l] = -1.0*height*dxfAt2D(2*i+ystart,2*j+xstart);
      a_peak->values[(CF_YCENTER+1)*size+l] = -1.0*height*dyfAt2D(2*i+ystart,2*j+xstart);
    }
  }
}

/*
 * updateFitValues3D()
 *
 * Use the cubic_spline library to update the fit values for a peak
 * based on the current parameters.
 *
 * a_peak - The peak to update.
 */
void updateFitValues3D(fitData *a_peak)
{
  int i,j,l,psx,psy,size,xstart,ystart,zi;
  double background,height,temp,xc,yc,zc;

  psx = a_peak->size_x;
  psy = a_peak->size_y;
  size = psx*psy;
  zi = a_peak->zi;

  /* Calculate new fit values. */
  height = a_peak->params[CF_HEIGHT];
  xc = 2.0*(2.0 - (a_peak->params[CF_XCENTER] - a_peak->xi));
  yc = 2.0*(2.0 - (a_peak->params[CF_YCENTER] - a_peak->yi));
  zc = (a_peak->params[CF_ZCENTER] - zi);

  xstart = 0;
  while(xc>1.0){
    xstart += 1;
    xc -= 1.0;
  }

  ystart = 0;
  while(yc>1.0){
    ystart += 1;
    yc -= 1.0;
  }

  computeDelta3D(zc, yc, xc);

  for(i=0;i<psy;i++){
    for(j=0;j<psx;j++){
      l = i*psx+j;
      temp = fAt3D(zi,2*i+ystart,2*j+xstart);
      a_peak->values[l] = height*temp;
      a_peak->values[(CF_HEIGHT+1)*size+l] = temp;
      a_peak->values[(CF_XCENTER+1)*size+l] = -1.0*height*dxfAt3D(zi,2*i+ystart,2*j+xstart);
      a_peak->values[(CF_YCENTER+1)*size+l] = -1.0*height*dyfAt3D(zi,2*i+ystart,2*j+xstart);
      a_peak->values[(CF_ZCENTER+1)*size+l] = 1.0*height*dzfAt3D(zi,2*i+ystart,2*j+xstart);
    }
  }
}

/*
 * updatePeakParameters()
 *
 * Update the peak parameters based on the update vector.
 *
 * a_peak - The peak to update.
 */
void updatePeakParameters(fitData *a_peak)
{
  int i,xc,yc;
  double dx,dy,update,maxz;

  /* Update current parameters. */
  if (DEBUG){
    printf("update: (%d)", a_peak->index);
  }
  for(i=0;i<a_peak->n_params;i++){

    /* Update sign & clamp if the solution appears to be oscillating. */
    if (a_peak->sign[i] != 0){
      if ((a_peak->sign[i] == 1) && (a_peak->delta[i] < 0.0)){
	a_peak->clamp[i] *= 0.5;
      }
      else if ((a_peak->sign[i] == -1) && (a_peak->delta[i] > 0.0)){
	a_peak->clamp[i] *= 0.5;
      }
    }
    if (a_peak->delta[i] > 0.0){
      a_peak->sign[i] = 1;
    }
    else {
      a_peak->sign[i] = -1;
    }

    /* Update values based on delta & clamp. */
    update = a_peak->delta[i]/(1.0 + fabs(a_peak->delta[i])/a_peak->clamp[i]);
    a_peak->params[i] -= update;
    if (DEBUG){
      printf(" %.3f", a_peak->params[i]);
      if (a_peak->sign[i] > 0){
	printf("+");
      }
      else{
	printf("-");
      }
    }
  }
  if (DEBUG){
    printf("\n");
  }

  /* 
   * Update peak position with hysteresis.
   */
  dx = a_peak->params[CF_XCENTER] - (double)a_peak->xi;
  if(dx < 0.25){
    a_peak->xi -= 1;
  }
  if(dx > 1.75){
    a_peak->xi += 1;
  }
  
  dy = a_peak->params[CF_YCENTER] - (double)a_peak->yi;
  if(dy < 0.25){
    a_peak->yi -= 1;
  }
  if(dy > 1.75){
    a_peak->yi += 1;
  }

  /*
  a_peak->xi = (int)(a_peak->params[CF_XCENTER]);
  a_peak->yi = (int)(a_peak->params[CF_YCENTER]);
  */

  /*
   * Check that the peak hasn't moved to close to the 
   * edge of the image. Flag the peak as bad if it has.
   */
  xc = a_peak->xi;
  yc = a_peak->yi;
  if((xc < 0)||(xc >= (image_size_x-a_peak->size_x))||(yc < 0)||(yc >= (image_size_y-a_peak->size_y))){
    a_peak->status = BADPEAK;
    if(TESTING){
      printf("object outside margins, %d, %d\n", xc, yc);
    }
  }
  
  /* 
   * Check for negative height. 
   */
  if(a_peak->params[CF_HEIGHT] < 0.0){
    a_peak->status = BADPEAK;
    if(TESTING){
      printf("negative height, %.3f\n", a_peak->params[CF_HEIGHT]);
    }
  }

  /* 
   * Check for z out of range. 
   */
  if(a_peak->type == F3D){
    if(a_peak->params[CF_ZCENTER] < 0.0){
      a_peak->params[CF_ZCENTER] = 0.0;
    }
    maxz = ((double)a_peak->size_z) - 1.0e-12;
    if(a_peak->params[CF_ZCENTER] > maxz){
      a_peak->params[CF_ZCENTER] = maxz;
    }
    a_peak->zi = (int)(a_peak->params[CF_ZCENTER]);
  }
}
