/*
 * Uses the homotopy C library for image analysis.
 *
 * FIXME? This assumes that L1FLT is of type double. 
 * If it is not then this will fail. Not sure if we
 * ever use floats though, so it is not obvious that
 * this is a problem.
 *
 * Hazen 04/13
 *
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "homotopy_common.h"
#include "homotopy_storm.h"
#include "homotopy_imagea_common.h"

/* Define */
#define WARNINGS 0
#define MAXITERS 800
#define NUMBERNONZERO 50

/* Structures */

/* Function Declarations */
void removePeak(double *, int *, int, int, double, double, double *, double *, double *);

/* Global Variables */
static int box_size;
static int box_size_hres;
static int image_x;
static int image_y;
static int overlap;
static int step;
static int xvec_size;
static int yvec_size;


/*
 * analyzeImage().
 *
 * Analyzes an entire image. It is assumed that the setImageParameters() has
 * already been called with an appropriate A matrix to use for
 * analysis of the image.
 *
 * hres - pre-allocated storage for the analysis results, should
 *        be initialized to zeros.
 * image - the image data.
 */
void analyzeImage(double *hres, double *image)
{
  int border,done,i,j,k,l,x_loc,x_start,y_loc,y_start;
  double alpha,lambda,tmp;

  /* initializations. */
  border = overlap*scale;
  done = 0;
  x_start = 0;
  y_start = 0;

  while(!done){
    
    /* copy image into yvec and compute l2 target. */
    alpha = 0.0;
    k = 0;
    for(i=0;i<box_size;i++){
      l = (y_start+i)*image_x + x_start;
      for(j=0;j<box_size;j++){
	tmp = image[l+j];
	alpha += fabs(tmp);
	yvec[k] = tmp;
	k++;
      }
    }

    /* calculate target alpha. */
    alpha = epsilon*sqrt(alpha);
    
    /* solve. */
    newYVector(yvec);
    lambda = solve(alpha, MAXITERS);
    if (WARNINGS && (lambda == 0.0)){
      printf(" Failed at %d %d %f %f\n", x_start, y_start, yvec[0], yvec[1]);
    }

    /* get solution. */
    for(i=0;i<xvec_size;i++){
      xvec[i] = 0.0;
    }
    getXVector(xvec);

    /* copy solution into hres vector. */
    x_loc = x_start*scale+border;
    y_loc = y_start*scale+border;
    for(i=0;i<(box_size_hres-2*border);i++){
      k = (y_loc+i)*hres_x;
      l = i*(box_size_hres-2*border);
      for(j=0;j<(box_size_hres-2*border);j++){
	hres[k+(x_loc+j)] = xvec[l+j];
      }
    }

    /* update x_start, y_start. */
    x_start += step;
    if(x_start==(image_x-2*overlap)){
      x_start = 0;
      y_start += step;
    }
    else if((x_start + box_size)>image_x){
      x_start = image_x - box_size;
    }

    if(y_start==(image_y-2*overlap)){
      done = 1;
    }
    else if((y_start + box_size)>image_y){
      y_start = image_y - box_size;
    }

  }
}

/*
 * setImageParameters()
 *
 * Initial setting of image parameters, since (for movies) 
 * they will all be the same..
 *
 * A - the PSF matrix, of size (M,N)
 *     M - box_size*box_size
 *     N - N ( >M )
 * n_cols_A - number of columns (i.e. N) in the A matrix.
 * has_background - 1 of A matrix includes a background term, zero otherwise.
 * pos_only - only add positive elements to the active set.
 * a_epsilon - l2 error scaling factor (see Bo's CS Nature Methods paper).
 * a_box_size - size of analysis sub-region.
 * a_image_x - image size in x.
 * a_image_y - image size in y.
 * a_overlap - amount of overlap between adjacent boxes.
 * a_scale - super-resolution overlapping.
 */
void setImageParameters(double *A, int n_cols_A, int has_background, int pos_only, double a_epsilon, int a_box_size, int a_image_x, int a_image_y, int a_overlap, int a_scale)
{
  /* Set some global variables. */
  epsilon = a_epsilon;

  box_size = a_box_size;
  image_x = a_image_x;
  image_y = a_image_y;
  overlap = a_overlap;
  scale = a_scale;

  box_size_hres = box_size*scale;

  hres_x = image_x*scale;
  hres_y = image_y*scale;
  
  step = box_size - 2*overlap;

  /* check that L1FLT is double and print a warning message if not. */
  if (sizeof(L1FLT)!=sizeof(double)){
    printf("!! L1FLT not equal to double !!\n");
  }

  /* Initialize homotopy library. */
  yvec_size = box_size*box_size;
  xvec_size = n_cols_A;
  if(has_background){
    initialize(A,yvec_size,xvec_size,1,pos_only,NUMBERNONZERO);
  }
  else{
    initialize(A,yvec_size,xvec_size,0,pos_only,NUMBERNONZERO);
  }
  
  /* Initialize storage. */
  xvec = (double *)malloc(sizeof(double)*xvec_size);
  yvec = (double *)malloc(sizeof(double)*yvec_size);
}

/*
 * See the accompanying license.txt file.
 *
 * Copyright (c) 2013 Zhuang Lab, Harvard University
 */
