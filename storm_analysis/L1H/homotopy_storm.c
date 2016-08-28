/*
 * C library for L1-homotopy continuation.
 *
 * This is based on:
 *   http://www.ceremade.dauphine.fr/~peyre/numerical-tour/tours/sparsity_2_sparsespikes_bp/#3
 *
 * And:
 *   Fast Solution of l1-norm Minimization Problems When the Solution May be Sparse
 *   David L. Donoho and Yaakov Tsaig, October 2006.
 *
 * This version has been cleaned up a bit and optimized for analysis of STORM images.
 * The assumptions here are:
 *   1. x is positive only.
 *   2. The last term of x is assumed to be the background term.
 *   3. The computation of G1 is shortcut based on values in C.
 *
 * Hazen 3/13
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "homotopy_common.h"
#include "homotopy_storm.h"

/* Define */
#define DEBUG 0
#define MAXGAMMA 1.0e+6
#define MINCFACTOR 0.5
#define PRECISION 1.0e-6
#define PROFILING 1
#define VERBOSE 0
#define VERYVERBOSE 0
#define WARNINGS 1

#define ONSET 0
#define OFFSET 1

/* Structures */

/* Function Declarations */
void computeC(void);
void computeD(void);
void computeG1(void);
void computeG3(void);
L1FLT l2Error(L1FLT *);
void update(void);
void updateOnIndices(void);

/* LAPACK Functions */
/* Cholesky solver */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);

/* Global Variables */
static int cholesky_failure;
static int max_c_index;
static int max_non_zero;
static int ncols;
static int nrows;
static int n_onset;
static int which_g1;
static int which_g3;

static int *offset;
static int *onset_indices;

static L1FLT gamma_g1;
static L1FLT gamma_g3;
static L1FLT gamma_s;
static L1FLT lambda;

static L1FLT *a_mat;
static L1FLT *a_mat_trans;
static L1FLT *a_y_vec;
static L1FLT *c_vec;
static L1FLT *d_vec;
static L1FLT *g_mat;
static L1FLT *x_vec;
static L1FLT *y_vec;

static L1FLT *work1;
static L1FLT *work2;

static double *double_work1;
static double *double_d_vec;

/*
 * cleanup()
 *
 * Frees all the allocated storage.
 */
void cleanup(void)
{
  freeCommon();

  free(offset);
  free(onset_indices);

  free(a_mat);
  free(a_mat_trans);
  free(a_y_vec);
  free(c_vec);
  free(d_vec);
  free(g_mat);
  free(x_vec);
  free(y_vec);

  free(work1);
  free(work2);

  free(double_work1);
  free(double_d_vec);
}

/*
 * computeC()
 *
 * Compute c vector, lambda and the location of the
 * maximum element in the c vector given current x vector.
 */
void computeC(void)
{
  int i,j,k;
  L1FLT tmp;

  if (PROFILING){
    startClock();
  }

  if (DEBUG){
    printf("(computeC) x_vec: ");
    for(i=0;i<n_onset;i++){
      printf("%d %f ", onset_indices[i], x_vec[onset_indices[i]]);
    }
    printf("\n");
  }

  if (1){
    for(i=0;i<ncols;i++){
      work1[i] = 0.0;
    }
    
    for(i=0;i<n_onset;i++){
      j = onset_indices[i];
      tmp = x_vec[j];
      j = j*ncols;
      for(k=0;k<ncols;k++){
	// using the fact that g_mat is symmetric (this goes across a row).
	work1[k] += g_mat[j+k] * tmp;
      }
    }
  }
  else{

    for(i=0;i<nrows;i++){
      work2[i] = 0.0;
    }

    for(i=0;i<n_onset;i++){
      j = onset_indices[i];
      tmp = x_vec[j];
      j = j*nrows;
      for(k=0;k<nrows;k++){
	work2[k] += a_mat_trans[j+k] * tmp;
      }
    }

    for(i=0;i<ncols;i++){
      work1[i] = 0.0;
      j = i*nrows;
      for(k=0;k<nrows;k++){
	work1[i] += a_mat_trans[j+k] * work2[k];
      }
    }
  }

  lambda = 0.0;
  max_c_index = 0;
  for(i=0;i<ncols;i++){
    c_vec[i] = a_y_vec[i] - work1[i];
    if(c_vec[i]>lambda){
      lambda = c_vec[i];
      max_c_index = i;
    }
  }

  if (DEBUG){
    printf("(computeC) onset size: %d\n", n_onset);
    for(i=0;i<ncols;i++){
      printf("%f ", c_vec[i]);
    }
    printf("\n");
  }

  if (PROFILING){
    stopClock(0);
  }
}

 /*
  * computeD()
  *
  * Compute the update direction vector.
  *
  * FIXME: 
  *  1. Work1 is symmetric, some calculations are redundant.
  *  2. Cholesky factorization rank-one updates.
  */
void computeD(void)
{
  int i,j;
  int n,nrhs,lda,ldb,info; /* for LAPACK */
  L1FLT tmp;

  if (PROFILING){
    startClock();
  }

  /* 
   * Compute A'[S] * A[S].
   * Since we already calculate A' * A, we just need to
   * pull the appropriate elements out of the G matrix.
   */
  for(i=0;i<n_onset;i++){
    for(j=0;j<n_onset;j++){
      double_work1[i*n_onset+j] = (double)g_mat[onset_indices[i]*ncols+onset_indices[j]];
    }
  }

  /* compute sign of c vector */
  for(i=0;i<n_onset;i++){
    tmp = c_vec[onset_indices[i]];
    if(tmp>0.0){
      double_d_vec[i] = 1.0;
    }
    else if(tmp<0.0){
      double_d_vec[i] = -1.0;
    }
    else{
      double_d_vec[i] = 0.0;
    }
  }

  if(n_onset==1){
    double_d_vec[0] = double_d_vec[0]/double_work1[0];
  }
  else{
    /* use LAPACK cholesky solver to determine d */
    nrhs = 1;
    n = lda = ldb = n_onset;
    dposv_( "Lower", &n, &nrhs, double_work1, &lda, double_d_vec, &ldb, &info );
    if((info!=0)&&(WARNINGS)){
      printf("  Cholesky solver failed with error: %d\n", info);
      cholesky_failure = 1;
    }
  }

  for(i=0;i<n_onset;i++){
    d_vec[i] = (L1FLT)double_d_vec[i];
  }

  if (DEBUG){
    printf("(computeD) c_vec (%d):\n", n_onset);
    for(i=0;i<ncols;i++){
      printf("%f ", c_vec[i]);
    }
    printf("\n");

    printf("(computeD) d_vec\n");
    for(i=0;i<n_onset;i++){
      printf("%f ", d_vec[i]);
    }
    printf("\n");

    printf("(computeD) onset_indices\n");
    for(i=0;i<n_onset;i++){
      printf("%d ", onset_indices[i]);
    }
    printf("\n");
  }

  if (PROFILING){
    stopClock(1);
  }
}

/*
 * computeG1()
 *
 * Compute gamma1 at current lambda value.
 */
void computeG1(void)
{
  int i,j,k;
  L1FLT w,tmp;

  if (PROFILING){
    startClock();
  }

  if (0){

    for(i=0;i<ncols;i++){
      work1[i] = 0.0;
    }
    
    for(i=0;i<n_onset;i++){
      j = ncols*onset_indices[i];
      tmp = d_vec[i];
      for(k=0;k<ncols;k++){
	work1[k] += g_mat[j+k] * tmp;
      }
    }

    gamma_g1 = MAXGAMMA;
    for(i=0;i<ncols;i++){
      if((offset[i]==OFFSET)&&(c_vec[i]>0.0)){
	tmp = (lambda - c_vec[i])/(1.0 - work1[i]);
	if((tmp>0.0)&&(tmp<gamma_g1)){
	  gamma_g1 = tmp;
	  which_g1 = i;
	}
      }
    }
  }
  else{

    /* compute v */

    for(i=0;i<nrows;i++){
      work1[i] = 0.0;
    }

    for(i=0;i<n_onset;i++){
      j = onset_indices[i]*nrows;
      for(k=0;k<nrows;k++){
	work1[k] += a_mat_trans[j+k]*d_vec[i];
      }
    }

    if (DEBUG){
      printf("(computeG1) work1: ");
      for(i=0;i<nrows;i++){
	printf("%f ", work1[i]);
      }
      printf("\n");
    }

    /* compute g1 */
    gamma_g1 = MAXGAMMA;
    for(i=0;i<ncols;i++){
      if((offset[i]==OFFSET)&&(c_vec[i]>0.0)){
	w = 0.0;
	j = i*nrows;
	for(k=0;k<nrows;k++){
	  w += a_mat_trans[j+k]*work1[k];
	}
	tmp = (lambda - c_vec[i])/(1.0 - w);
	if((tmp>0.0)&&(tmp<gamma_g1)){
	  gamma_g1 = tmp;
	  which_g1 = i;
	}
      }
    }
  }

  if (PROFILING){
    stopClock(2);
  }
}

/*
 * computeG3()
 *
 * Compute gamma3 vector.
 */
void computeG3(void)
{
  int i,j;
  L1FLT tmp;

  if (PROFILING){
    startClock();
  }

  gamma_g3 = MAXGAMMA;
  for(i=0;i<(n_onset-1);i++){
    j = onset_indices[i];
    tmp = -1.0*x_vec[j]/d_vec[i];
    if((tmp>0.0)&&(tmp<gamma_g3)){
      gamma_g3 = tmp;
      which_g3 = j;
    }
  }

  if (PROFILING){
    stopClock(3);
  }
}

/*
 * getXVector()
 *
 * Returns the current x vector in user supplied storage.
 *
 * xvec - user supplied storage for the x vector, assumed initialized to zero.
 */
void getXVector(L1FLT *xvec)
{
  int i,j;

  for(i=0;i<n_onset;i++){
    j = onset_indices[i];
    xvec[j] = x_vec[j];
  }
}

/*
 * initialize()
 *
 * Allocate space for computations and perform some preliminary calculations.
 *
 * A - the A matrix.
 * rows - the number of rows in the A matrix.
 * cols - the number of columns in the A matrix.
 *  (ignored) bgterm - the last column is a background term (i.e. not subject to compression, not removable).
 *  (ignored) pos_only - only add positive elements to the active set.
 * number_nonzero - the maximum number of non zero elements in the final solution.
 *
 * Since the problem is presumably under-determined rows <= cols.
 */
void initialize(L1FLT *A, int rows, int cols, int bgterm, int pos_only, int number_nonzero)
{
  int i,j,k,l;
  L1FLT sum;

  max_non_zero = number_nonzero;
  ncols = cols;
  nrows = rows;

  /* Initialize that which is common to all solvers.. */
  initCommon();

  /* Allocate memory. */
  offset = (int *)malloc(sizeof(int)*ncols);
  onset_indices = (int *)malloc(sizeof(int)*(max_non_zero+1));

  a_mat = (L1FLT *)malloc(sizeof(L1FLT)*nrows*ncols);
  a_mat_trans = (L1FLT *)malloc(sizeof(L1FLT)*nrows*ncols);
  a_y_vec = (L1FLT *)malloc(sizeof(L1FLT)*ncols);
  c_vec = (L1FLT *)malloc(sizeof(L1FLT)*ncols);
  d_vec = (L1FLT *)malloc(sizeof(L1FLT)*max_non_zero);
  g_mat = (L1FLT *)malloc(sizeof(L1FLT)*ncols*ncols);
  x_vec = (L1FLT *)malloc(sizeof(L1FLT)*ncols);
  y_vec = (L1FLT *)malloc(sizeof(L1FLT)*nrows);
  double_d_vec = (double *)malloc(sizeof(double)*max_non_zero);

  if((max_non_zero*max_non_zero)>ncols){
    work1 = (L1FLT *)malloc(sizeof(L1FLT)*max_non_zero*max_non_zero);
    work2 = (L1FLT *)malloc(sizeof(L1FLT)*max_non_zero*max_non_zero);
    double_work1 = (double *)malloc(sizeof(double)*max_non_zero*max_non_zero);
  }
  else{
    work1 = (L1FLT *)malloc(sizeof(L1FLT)*ncols);
    work2 = (L1FLT *)malloc(sizeof(L1FLT)*ncols);
    double_work1 = (double *)malloc(sizeof(double)*ncols);
  }

  /* Initializations and preliminary calculations. */

  for(i=0;i<nrows;i++){
    for(j=0;j<ncols;j++){
      a_mat[i*ncols+j] = A[i*ncols+j];
      a_mat_trans[j*nrows+i] = A[i*ncols+j];
    }
  }

  /* G matrix, A' * A */
  for(i=0;i<ncols;i++){
    for(j=0;j<ncols;j++){
      sum = 0.0;
      for(k=0;k<nrows;k++){
	l = k*ncols;
	sum += a_mat[l+i]*a_mat[l+j];
      }
      g_mat[i*ncols+j] = sum;
    }
  }

}

/*
 * l2Error()
 *
 * Calculate the l2Error given a x vector.
 *
 * xvec - a x_vector
 */
L1FLT l2Error(L1FLT *xvec)
{
  int i,j,k;
  L1FLT sum,tmp;

  if (PROFILING){
    startClock();
  }

  /* Calculate Ax */
  for(i=0;i<nrows;i++){
    work1[i] = 0.0;
  }

  for(i=0;i<n_onset;i++){
    j = onset_indices[i];
    tmp = xvec[j];
    j = j*nrows;
    for(k=0;k<nrows;k++){
      work1[k] += a_mat_trans[j+k]*tmp;
    }
  }

  /* Calculate error */
  sum = 0.0;
  for(i=0;i<nrows;i++){
    tmp = work1[i] - y_vec[i];
    sum += tmp*tmp;
  }

  if (PROFILING){
    stopClock(4);
  }

  return sum;
}

/*
 * newYVector()
 *
 * Setup problem with a new y vector. This does the following:
 *  1. Copy the y vector into local storage.
 *  2. Compute the product Ay and save this (this is also the initial C vector).
 *  3. Compute the maximum element of Ay and its location.
 *  4. Initialize offset vectors.
 *
 * yvec - the new y vector (assumed of size nrows).
 */
void newYVector(L1FLT *yvec)
{
  int i,j,k;
  L1FLT sum;

  if (PROFILING){
    startClock();
  }

  /* 1. Copy into y_vec. */
  for(i=0;i<nrows;i++){
    y_vec[i] = yvec[i];
  }

  /* 2/3. Calculate a_y_vec, c_vec, find maximum. */
  lambda = 0.0;
  max_c_index = 0;
  for(i=0;i<ncols;i++){
    sum = 0.0;
    j = i*nrows;
    for(k=0;k<nrows;k++){
      sum += a_mat_trans[j+k]*y_vec[k];
    }
    a_y_vec[i] = sum;
    c_vec[i] = sum;
    if(sum>lambda){
      lambda = sum;
      max_c_index = i;
    }
  }

  if(VERYVERBOSE){
    printf("first nrows elements of ay, c, y:\n");
    for(i=0;i<nrows;i++){
      printf(" %d %.20f %.20f %.20f\n", i, a_y_vec[i], c_vec[i], y_vec[i]);
    }
    printf("\n");
  }

  /* 4. Initialize offset vector. */
  for(i=0;i<ncols;i++){
    offset[i] = OFFSET;
  }
  offset[max_c_index] = ONSET;

  if (PROFILING){
    stopClock(5);
  }

}

/*
 * solve()
 *
 * Follow solution path until l2 error is equal to epsilon.
 *
 * epsilon - desired l2 error.
 * max_iters - maximum number of iterations before giving up.
 */
L1FLT solve(L1FLT epsilon, int max_iters)
{
  int i,is_bad,iters,j;
  L1FLT cur_error,low,high,mid;

  if (PROFILING){
    startClock();
  }

  cholesky_failure = 0;
  epsilon = epsilon*epsilon;
  is_bad = 0;
  n_onset = 0;

  for(i=0;i<ncols;i++){
    x_vec[i] = 0.0;
  }

  if(VERBOSE){
    printf(" l2 error target: %f\n\n", epsilon);
  }

  /* Check if we need to do anything. */
  cur_error = l2Error(x_vec);
  if(cur_error<epsilon){
    return lambda;
  }

  /* Update onset and offset using information calculated in newYVector(). */
  updateOnIndices();

  /* Follow solution until l2 error is less than epsilon. */
  iters = 0;
  while(((cur_error-epsilon)>PRECISION)&&(n_onset<max_non_zero)&&(iters<max_iters)&&(!cholesky_failure)){

    if(VERYVERBOSE){
      printf("iteration %d start (%d %.20f %.20f)\n",iters,max_c_index,lambda,l2Error(x_vec));
    }

    computeD();
    computeG1();
    computeG3();
    update();

    cur_error = l2Error(x_vec);

    /* Report current state. */
    if(VERYVERBOSE){
      printf(" iteration %d end (%d, %.20f, %.20f)\n",iters,n_onset,lambda,l2Error(x_vec));

      printf("  On set: <index> <on_index> <x> <d_vec>\n");
      for(i=0;i<n_onset;i++){
	printf("   %d %d %.14f %.14f\n",i,onset_indices[i],x_vec[onset_indices[i]],d_vec[i]);
      }
      printf("\n");
    }

    if((cur_error-epsilon)>PRECISION){
      updateOnIndices();
      computeC();
    }
    
    iters++;
  }

  if (PROFILING){
    updateIterations(iters);
  }

  updateFailureCounter(0);

  /* Print warning message in the event of too many non-zero elements. */
  if(n_onset>=max_non_zero){
    if(WARNINGS){
      printf("Max non zero value reached, solution is not optimal.\n");
      printf("  %d elements in the on set\n", n_onset);
      printf("  %d iterations performed\n", iters);
      printf("  current lambda is %f\n", lambda);
      printf("  current error is %f (%f target)\n\n", cur_error, epsilon);
    }
    updateFailureCounter(1);
    is_bad = 1;
  }

  /* Print warning message in the event that the 
     maximum number of iterations was reached. */
  if(iters==max_iters){
    if(WARNINGS){
      printf("Maximum number of iterations reached, solution is not optimal.\n");
      printf("  %d elements in the on set\n", n_onset);
      printf("  %d iterations performed\n", iters);
      printf("  current lambda is %f\n", lambda);
      printf("  current error is %f (%f target)\n\n", cur_error, epsilon);
    }
    updateFailureCounter(2);
    is_bad = 1;
  }

  /* Print cholesky failure warning */
  if(cholesky_failure){
    if(WARNINGS){
      printf(" Cholesky failure\n");
    }
    updateFailureCounter(3);
    is_bad = 1;
  }

  if(!is_bad){

    /* Mid-point bisection to find optimal gamma value. */
    low = 0.0;
    high = gamma_s;
    for(i=0;i<ncols;i++){
      work2[i] = 0.0;
    }
    while(((high-low)/gamma_s)>1.0e-3){
      mid = 0.5*(high+low);
      if(VERYVERBOSE){
	printf(" mpd: %f %f %f %f\n", low, mid, high, l2Error(work2));
      }
      for(i=0;i<n_onset;i++){
	j = onset_indices[i];
	work2[j] = x_vec[j] - d_vec[i]*mid;
      }
      if(l2Error(work2)>epsilon){
	high = mid;
      }
      else{
	low = mid;
      }
    }

    /* Update x vector and calculate corresponding lambda value. */
    for(i=0;i<n_onset;i++){
      j = onset_indices[i];
      x_vec[j] = work2[j];
    }
    computeC();
  }
  else{
    n_onset = 0;
    lambda = 0.0;
  }

  /* Print solution information. */
  if(VERBOSE){
    if(is_bad){
      printf(" Failed to find solution\n\n");
    }
    else{
      printf(" (Storm) Solved in %d iterations, %d non-zero, lambda = %f, l2_error = %f\n", iters, n_onset, lambda, l2Error(x_vec));
      printf(" On set: <index> <on_index> <x>\n");
      for(i=0;i<n_onset;i++){
	printf("  %d %d %f\n",i,onset_indices[i],x_vec[onset_indices[i]]);
      }
      printf("\n");
    }
  }

  if (PROFILING){
    stopClock(6);
  }

  return lambda;
}

/*
 * update()
 *
 * Finds the minimum value (> 0.0) of the g1, g2 and g3 vectors.
 * Updates the offset vector am_off accordingly.
 */
void update(void)
{
  int i,index,which_g;

  which_g = 1;
  index = which_g1;
  gamma_s = gamma_g1;

  if(gamma_g3<gamma_g1){
    which_g = 3;
    index = which_g3;
    gamma_s = gamma_g3;
  }

  /* update x vector */
  for(i=0;i<n_onset;i++){
    x_vec[onset_indices[i]] += gamma_s*d_vec[i];
  }

  /* update am_off */
  if(which_g==3){
    if(VERYVERBOSE){
      printf(" Removing: %d %.20f\n",index,gamma_s);
    }
    offset[index] = OFFSET;
  }
  else{
    if(VERYVERBOSE){
      printf(" Adding: %d %.20f\n",index,gamma_s);
    }
    offset[index] = ONSET;
  }
}

/*
 * updateOnIndices()
 *
 * Updates the on_indices array and number_on.
 */
void updateOnIndices(void)
{
  int i;

  n_onset = 0;
  for(i=0;i<ncols;i++){
    if(offset[i]==ONSET){
      onset_indices[n_onset] = i;
      n_onset++;
    }
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
