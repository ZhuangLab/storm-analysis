/*
 * Fit multiple, possibly overlapping, functions simultaneously to 
 * image data.
 *
 * The MLE fitter is based on the approach described in
 * Laurence and Chromy, Nature Methods, 2010.
 *
 * This is an overhaul and re-factoring of multi_fit.c to (primarily)
 * make it more general. At some point in the future this may replace
 * multi_fit.c. As with multi_fit.c this library is not thread safe.
 *
 * Expected program flow:
 *  1. Call initializeMultiFit() with sCMOS calibration data (or an 
 *     array of zeros) if this is not relevant, and the image size.
 *
 *  2. Call newImage() with the image to analyze.
 *
 *  3. Initialize the fitData structures & subtract initial peaks 
 *     from the image.
 *
 *  4. ..Fitting.. addPeak(), calcError(), updateFit(), subtractPeak().
 *
 *  5. Get updated peak parameters.
 *
 *  6. Call getResidual() to get the residual image.
 *
 *  7. Free fitData structures with fitDataFree(). The image and
 *     residual arrays should not be freed as they are recycled.
 *
 *  8. Repeat 2-7 to analyze additional frames (or to re-analyze the
 *     current frame with new peak parameters).
 *
 *  9. Call multiFitCleanup().
 * 
 * Hazen 12/13
 *
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall multi_fit_core.c
 *  gcc -shared -Wl,-soname,multi_fit_core.so.1 -o multi_fit_core.so.1.0.1 multi_fit_core.o -lc -llapack
 *
 * Windows:
 *  gcc -c multi_fit_core.c
 *  gcc -shared -o multi_fit_core.dll multi_fit_core.o -llapack -Lc:\Users\Hazen\lib
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "multi_fit_core.h"

/* Define */
#define TESTING 0

/* (External) LAPACK functions. */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);

/* (External) Global variables. */
int image_size_x;    /* Image size in x (fast axis). */
int image_size_y;    /* Image size in y (slow axis). */
int n_fit_data;      /* Number of fitData structures. */
double tolerance;    /* Fit tolerance for convergence. */
fitData *fit_data;   /* fitData structures. */

/* (Internal) Global variables. */
static double *image;       /* Image. */
static double *fit;         /* Fit array. */
static double *scmos_data;  /* sCMOS calibration data for each pixel (var/gain^2). */


/*
 * addPeak()
 *
 * Add the peak back to fit array.
 *
 * a_peak - The peak to add to the fit array.
 */
void addPeak(fitData *a_peak)
{
  int i,j,k,psx;
  //double bg_term;

  psx = a_peak->size_x;
  i = a_peak->yi * image_size_x + a_peak->xi;
  for (j=0;j<a_peak->size_y;j++){
    for (k=0;k<psx;k++){
      fit[i+j*image_size_x+k] += a_peak->values[j*psx+k];
    }
  }
}

/*
 * calcError()
 *
 * Calculate the fit error of the peak. Technically this is
 * actually the total error in the pixels that are covered
 * by the peak. When peaks overlap substantially they will
 * have similar errors.
 *
 * a_peak - The peak to calculate the fit error of.
 * background - The background value of the peak.
 */
double calcError(fitData *a_peak, double background)
{
  int i,j,k,psx;
  double err,fi,xi;

  err = 0.0;
  psx = a_peak->size_x;
  i = a_peak->yi * image_size_x + a_peak->xi;
  for (j=0;j<a_peak->size_y;j++){
    for (k=0;k<psx;k++){
      fi = fit[i+j*image_size_x+k] + background;
      if (fi <= 0.0){
	if (TESTING){
	  printf(" Negative f detected! %.3f %d %d %d\n", fi, j, k, a_peak->index);
	}
	a_peak->status = BADPEAK;
	return -1.0;
      }
      xi = image[i+j*image_size_x+k];
      err += 2.0*((fi-xi)-xi*log(fi/xi));
    }
  }
  return err;
}

/*
 * copyFitData()
 *
 * Copy a fitData structure.
 *
 * original - The original structure.
 * copy - The copy.
 */
void copyFitData(fitData *original, fitData *copy)
{
  int i;

  copy->index = original->index;
  copy->n_params = original->n_params;
  copy->size_x = original->size_x;
  copy->size_y = original->size_y;
  copy->status = original->status;
  copy->xi = original->xi;
  copy->yi = original->yi;
  copy->error = original->error;
  copy->lambda = original->lambda;

  for(i=0;i<original->n_params;i++){
    copy->sign[i] = original->sign[i];
    copy->clamp[i] = original->clamp[i];
    copy->delta[i] = original->delta[i];
    copy->params[i] = original->params[i];
  }

  for(i=0;i<(original->size_x*original->size_y*(original->n_params+1));i++){
    copy->values[i] = original->values[i];
  }
}

/*
 * freeFitData()
 *
 * Free a fitData structure.
 *
 * a_peak - The fitData structure to free.
 */
void freeFitData(fitData *a_peak)
{
  free(a_peak->sign);
  free(a_peak->clamp);
  free(a_peak->delta);
  free(a_peak->params);
  free(a_peak->values);
}

/*
 * getResidual()
 *
 * Return the residual image in user supplied storage. The
 * residual includes the background. It is the image minus
 * the fit peaks.
 *
 * res - Storage for the residual array.
 */
void getResidual(double *res)
{
  int i;

  for (i=0;i<(image_size_y*image_size_x);i++){
    res[i] = image[i] - fit[i];
  }
}

/*
 * initializeMultiFit()
 *
 * Initialize storage and the scmos_data array.
 *
 * new_scmos_data - New sCMOS calibration data.
 * new_tolerance - New fit tolerance value.
 * new_image_size_x - New image size in x.
 * new_image_size_y - New image size in y.
 */
void initializeMultiFit(double *new_scmos_data, double new_tolerance, int new_image_size_x, int new_image_size_y)
{
  int i;

  tolerance = new_tolerance;
  image_size_x = new_image_size_x;
  image_size_y = new_image_size_y;

  image = (double *)malloc(sizeof(double)*image_size_x*image_size_y);
  fit = (double *)malloc(sizeof(double)*image_size_x*image_size_y);
  scmos_data = (double *)malloc(sizeof(double)*image_size_x*image_size_y);

  for (i=0;i<(image_size_x*image_size_y);i++){
    scmos_data[i] = new_scmos_data[i];
  }
}

/*
 * mallocFitData()
 *
 * Allocate pointers in a fitData structure.
 *
 * a_peak - The fitData structure.
 * n_params - The (maximum) number of parameters.
 * n_vals - The (maximum) number of fit_vals.
 */
void mallocFitData(fitData *a_peak, int n_params, int n_vals)
{
  a_peak->sign = (int *)malloc(sizeof(int)*n_params);
  a_peak->clamp = (double *)malloc(sizeof(double)*n_params);
  a_peak->delta = (double *)malloc(sizeof(double)*n_params);
  a_peak->params = (double *)malloc(sizeof(double)*n_params);

  a_peak->values = (double *)malloc(sizeof(double)*n_vals);
}

/*
 * multiFitCleanup()
 *
 * Free (locally) allocated storage.
 */
void multiFitCleanup(void)
{
  if (image != NULL){
    free(image);
    free(fit);
    free(scmos_data);
  }
}

/*
 * newImage()
 *
 * Initialize the image.
 *
 * new_image - New image.
 */
void newImage(double *new_image)
{
  int i;

  for (i=0;i<(image_size_x*image_size_y);i++){
    image[i] = new_image[i];
    fit[i] = 0.0;
  }
}

/*
 * resetFit()
 *
 * Zero the fit array so that you can fit the same image with
 * a new list of peaks.
 */
void resetFit()
{
  int i;

  for (i=0;i<(image_size_x*image_size_y);i++){
    fit[i] = 0.0;
  }
}

/*
 * subtractPeak()
 *
 * Subtract the peak from the image residual.
 *
 * a_peak - The peak to add back to the image residual.
 */
void subtractPeak(fitData *a_peak)
{
  int i,j,k,psx;

  psx = a_peak->size_x;
  i = a_peak->yi * image_size_x + a_peak->xi;
  for (j=0;j<a_peak->size_y;j++){
    for (k=0;k<psx;k++){
      fit[i+j*image_size_x+k] -= a_peak->values[j*psx+k];
    }
  }
}

/*
 * updateFit()
 *
 * Performs (part of) one iteration of peak fitting for the peak.
 *
 * 1. Calculate jacobian and hessian.
 *
 * 2. Calculate update term using the LAPACK dposv_ function.
 *
 *
 * Expected program flow:
 *  1. The calling function adds the peak to the residual image.
 * 
 *  2. The calling function calls calcError() to determine if the
 *     peak parameters need additional refinement.
 *
 *  3. This function calculates the update term for additional
 *     refinement of the peak parameters.
 *
 *  4. The calling function subtracts the original peak from the 
 *     residual image before calling this function again with 
 *     another peak.
 *
 *  5. The calling function updates the peak fitData based on the
 *     update term.
 *
 * a_peak - (input) Peak to update.
 * background - The background value of the peak.
 */
void updateFit(fitData *a_peak, double background)
{
  /* Lapack variables */
  int info, lda, ldb, nrhs;

  /* Local */
  int i,j,k,l,m,n,np,psx,size;
  double fi,xi,t1,t2;
  double df_dp[MAXPARAMS];
  double jacobian[MAXPARAMS];
  double hessian[MAXPARAMS*MAXPARAMS];

  np = a_peak->n_params;

  /* Zero jacobian and hessian. */
  for(i=0;i<np;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<(np*np);i++){
    hessian[i] = 0.0;
  }

  /* Calculate jacobian and hessian. */
  size = a_peak->size_y*a_peak->size_x;
  psx = a_peak->size_x;
  i = a_peak->yi * image_size_x + a_peak->xi;
  for(j=0;j<a_peak->size_y;j++){
    for(k=0;k<psx;k++){

      /* setup */
      l = j*psx+k;
      xi = image[i+j*image_size_x+k];
      fi = fit[i+j*image_size_x+k] + background;
      for(m=0;m<np;m++){
	df_dp[m] = a_peak->values[(m+1)*size+l];
      }

      /* jacobian */
      t1 = (1.0 - xi/fi);
      for(l=0;l<np;l++){
	jacobian[l] += t1*df_dp[l];
      }

      /* hessian */
      t2 = xi/(fi*fi);
      for(l=0;l<np;l++){
	for(m=l;m<np;m++){
	  hessian[l*np+m] += t2*df_dp[l]*df_dp[m];
	}
      }
    }
  }

  /* Multiply diagonal terms by lambda. */
  for(i=0;i<np;i++){
    hessian[i*np+i] = hessian[i*np+i] * a_peak->lambda;
  }

  /* Print jacobian and hessian. */
  if(0){
    printf("Jacobian:\n");
    for(i=0;i<np;i++){
      printf(" %.3f", jacobian[i]);
    }
    printf("\n\n");
    printf("Hessian:\n");
    for(i=0;i<np;i++){
      for(j=0;j<np;j++){
	printf(" %.3f", hessian[i*np+j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  /* Lapack setup */
  n = np;
  nrhs = 1;
  lda = np;
  ldb = np;

  /* Use Lapack to solve Ax=B to calculate update vector. */
  dposv_( "Lower", &n, &nrhs, hessian, &lda, jacobian, &ldb, &info );
  
  if(info!=0){
    a_peak->status = CHOLERROR;
    if(TESTING){
      printf("Cholesky error for peak %d\n", a_peak->index);
    }
  }

  /* Copy results into peak update vector. */
  for (i=0;i<a_peak->n_params;i++){
    a_peak->delta[i] = jacobian[i];
  }
}
