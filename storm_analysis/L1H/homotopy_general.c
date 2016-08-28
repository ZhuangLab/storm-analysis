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
 * This version is "general" in that it should be able to solve any problem that any other
 * homotopy solver could solve.
 *
 * Hazen 8/12
 *
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall homotopy_general.c
 *  gcc -shared -Wl,-soname,homotopy.so.1 -o homotopy.so.1.0.1 homotopy_general.o -lc -llapack
 *
 * Windows:
 *  gcc -c homotopy_general.c
 *  gcc -shared -o homotopy_general.dll homotopy_general.o -llapack
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "homotopy_common.h"
#include "homotopy_storm.h"

/* Define */
#define MAXGAMMA 1.0e+6
#define PRECISION 1.0e-6
#define VERBOSE 0
#define VERYVERBOSE 0
#define WARNINGS 1

/* Structures */

/* Function Declarations */
void computeC(void);
void computeD(void);
void computeG1G2(void);
void computeG3(void);
int getL1FLTSize(void);
double l2Error(double *);
void update(void);
void updateOnIndices(void);

/* LAPACK Functions */
/* Cholesky solver */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);

/* Global Variables */
static int background_term;
static int cholesky_failure;
static int max_c_index;
static int max_non_zero;
static int ncols;
static int nrows;
static int number_on;
static int positive_only;

static int *am_off;
static int *on_indices;
static int *visited;

static double gamma_s;
static double lambda;

static double *a_mat;
static double *a_y_vec;
static double *c_vec;
static double *d_vec;
static double *g_mat;
static double *g1_vec;
static double *g2_vec;
static double *g3_vec;
static double *x_vec;
static double *y_vec;

static double *work1;
static double *work2;


/*
 * cleanup()
 *
 * Frees all the allocated storage.
 */
void cleanup(void)
{
  freeCommon();

  free(am_off);
  free(on_indices);
  free(visited);

  free(a_mat);
  free(a_y_vec);
  free(c_vec);
  free(d_vec);
  free(g_mat);
  free(g1_vec);
  free(g2_vec);
  free(g3_vec);
  free(x_vec);
  free(y_vec);

  free(work1);
  free(work2);
}

/*
 * computeC()
 *
 * Compute c vector, lambda and the location of the
 * maximum element in the c vector given current x vector.
 */
void computeC(void)
{
  int i,j,k,l;
  
  for(i=0;i<ncols;i++){
    work1[i] = 0.0;
  }
  
  for(i=0;i<number_on;i++){
    j = on_indices[i];
    k = j*ncols;
    for(l=0;l<ncols;l++){
      // using the fact that g_mat is symmetric (this goes across a row).
      work1[l] += g_mat[k+l]*x_vec[j];
    }
  }

  lambda = 0.0;
  max_c_index = 0;
  if(positive_only){
    for(i=0;i<ncols;i++){
      c_vec[i] = a_y_vec[i] - work1[i];
      if(c_vec[i]>lambda){
	lambda = c_vec[i];
	max_c_index = i;
      }
    }
  }
  else{
    for(i=0;i<ncols;i++){
      c_vec[i] = a_y_vec[i] - work1[i];
      if(fabs(c_vec[i])>lambda){
	lambda = fabs(c_vec[i]);
	max_c_index = i;
      }
    }
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
  double tmp;

  /* 
   * Compute A'[S] * A[S].
   * Since we already calculate A' * A, we just need to
   * pull the appropriate elements out of the G matrix.
   */
  for(i=0;i<number_on;i++){
    for(j=0;j<number_on;j++){
      work1[i*number_on+j] = g_mat[on_indices[i]*ncols+on_indices[j]];
    }
  }

  /* compute sign of c vector */
  for(i=0;i<number_on;i++){
    tmp = c_vec[on_indices[i]];
    if(tmp>0.0){
      d_vec[i] = 1.0;
    }
    else if(tmp<0.0){
      d_vec[i] = -1.0;
    }
    else{
      d_vec[i] = 0.0;
    }
  }

  if(number_on==1){
    d_vec[0] = d_vec[0]/work1[0];
  }
  else{
    /* use LAPACK cholesky solver to determine d */
    nrhs = 1;
    n = lda = ldb = number_on;
    dposv_( "Lower", &n, &nrhs, work1, &lda, d_vec, &ldb, &info );
    if((info!=0)&&(WARNINGS)){
      printf("  Cholesky solver failed with error: %d\n", info);
      cholesky_failure = 1;
    }
  }
}

/*
 * computeG1G2()
 *
 * Compute gamma1 and gamma2 vectors at current lambda value.
 */
void computeG1G2(void)
{
  int i,j,k;
  double w;

  /* compute v */
  for(i=0;i<nrows;i++){
    work1[i] = 0.0;
  }

  for(i=0;i<number_on;i++){
    j = on_indices[i];
    for(k=0;k<nrows;k++){
      work1[k] += a_mat[k*ncols+j] * d_vec[i];
    }
  }

  /*
  for(i=0;i<nrows;i++){
    work1[i] = 0.0;
    j = i*ncols;
    for(k=0;k<number_on;k++){
      l = on_indices[k];
      work1[i] += a_mat[j+l]*d_vec[k];
    }
  }
  */

  /* compute w, g1_vec, g2_vec */
  if(positive_only){
    for(i=0;i<ncols;i++){
      g1_vec[i] = w = 0.0;
      if(am_off[i]){
	for(j=0;j<nrows;j++){
	  w += a_mat[j*ncols+i]*work1[j];
	}
	g1_vec[i] = (lambda - c_vec[i])/(1.0 - w);
      }
    }
  }
  else{
    for(i=0;i<ncols;i++){
      g1_vec[i] = g2_vec[i] = w = 0.0;
      if(am_off[i]){
	for(j=0;j<nrows;j++){
	  w += a_mat[j*ncols+i]*work1[j];
	}
	g1_vec[i] = (lambda - c_vec[i])/(1.0 - w);
	g2_vec[i] = (lambda + c_vec[i])/(1.0 + w);
      }
    }
  }
}

/*
 * computeG3()
 *
 * Compute gamma3 vector.
 */
void computeG3(void)
{
  int i;

  for(i=0;i<number_on;i++){
    g3_vec[i] = -1.0*x_vec[on_indices[i]]/d_vec[i];
  }

  if(background_term){
    g3_vec[number_on-1] = 0.0;
  }
}

/*
 * getVisited()
 *
 * Returns the current visited vector in user supplied storage.
 *
 * vis - user supplied storage for the visited vector.
 */
void getVisited(int *vis)
{
  int i;

  for(i=0;i<ncols;i++){
    vis[i] = visited[i];
  }
}

/*
 * getXVector()
 *
 * Returns the current x vector in user supplied storage.
 *
 * xvec - user supplied storage for the x vector, assumed initialized to zero.
 */
void getXVector(double *xvec)
{
  int i,j;

  for(i=0;i<number_on;i++){
    j = on_indices[i];
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
 * bgterm - the last column is a background term (i.e. not subject to compression, not removable).
 * pos_only - only add positive elements to the active set.
 * number_nonzero - the maximum number of non zero elements in the final solution.
 *
 * Since the problem is presumably under-determined rows <= cols.
 */
void initialize(double *A, int rows, int cols, int bgterm, int pos_only, int number_nonzero)
{
  int i,j,k,l;
  double sum;

  /* Initialize that which is common to all solvers.. */
  initCommon();

  background_term = bgterm;
  max_non_zero = number_nonzero;
  ncols = cols;
  nrows = rows;
  positive_only = pos_only;

  /* Allocate memory. */
  am_off = (int *)malloc(sizeof(int)*ncols);
  on_indices = (int *)malloc(sizeof(int)*max_non_zero);
  visited = (int *)malloc(sizeof(int)*ncols);

  a_mat = (double *)malloc(sizeof(double)*nrows*ncols);
  a_y_vec = (double *)malloc(sizeof(double)*ncols);
  c_vec = (double *)malloc(sizeof(double)*ncols);
  d_vec = (double *)malloc(sizeof(double)*max_non_zero);
  g_mat = (double *)malloc(sizeof(double)*ncols*ncols);
  g1_vec = (double *)malloc(sizeof(double)*ncols);
  g2_vec = (double *)malloc(sizeof(double)*ncols);
  g3_vec = (double *)malloc(sizeof(double)*max_non_zero);
  x_vec = (double *)malloc(sizeof(double)*ncols);
  y_vec = (double *)malloc(sizeof(double)*nrows);

  if((max_non_zero*max_non_zero)>ncols){
    work1 = (double *)malloc(sizeof(double)*max_non_zero*max_non_zero);
    work2 = (double *)malloc(sizeof(double)*max_non_zero*max_non_zero);
  }
  else{
    work1 = (double *)malloc(sizeof(double)*ncols);
    work2 = (double *)malloc(sizeof(double)*ncols);
  }

  /* Initializations and preliminary calculations. */
  for(i=0;i<(nrows*ncols);i++){
    a_mat[i] = A[i];
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
double l2Error(double *xvec)
{
  int i,j,l;
  double sum,tmp;

  /* Calculate Ax */
  for(i=0;i<nrows;i++){
    work1[i] = 0.0;
  }

  for(i=0;i<number_on;i++){
    l = on_indices[i];
    tmp = xvec[l];
    for(j=0;j<nrows;j++){
      work1[j] += tmp*a_mat[j*ncols+l];
    }
  }

  /* Calculate error */
  sum = 0.0;
  for(i=0;i<nrows;i++){
    tmp = work1[i] - y_vec[i];
    sum += tmp*tmp;
  }

  return sum;
}

/*
 * newYVector()
 *
 * Setup problem with a new y vector.
 *
 * yvec - the new y vector (assumed of size nrows).
 */
void newYVector(double *yvec)
{
  int i,j;
  double sum;

  /* copy into y_vec */
  for(i=0;i<nrows;i++){
    y_vec[i] = yvec[i];
  }

  /* calculate a_y_vec */
  for(i=0;i<ncols;i++){
    sum = 0.0;
    for(j=0;j<nrows;j++){
      sum += a_mat[j*ncols+i]*y_vec[j];
    }
    a_y_vec[i] = sum;
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
double solve(double epsilon, int max_iters)
{
  int i,is_bad,iters,j;
  double cur_error,low,high,mid;

  /* Initialization. */
  is_bad = 0;
  cholesky_failure = 0;
  epsilon = epsilon*epsilon;

  if(VERBOSE){
    printf(" l2 error target: %f\n\n", epsilon);
  }

  number_on = 0;

  for(i=0;i<ncols;i++){
    am_off[i] = 1;
    visited[i] = 0;
    x_vec[i] = 0.0;
  }

  for(i=0;i<max_non_zero;i++){
    on_indices[i] = 0;
  }

  /* Find starting lambda value. */
  computeC();

  /* Check if we need to do anything. */
  cur_error = l2Error(x_vec);
  if(cur_error<epsilon){
    return lambda;
  }

  /* Add first point to the active set. */
  am_off[max_c_index] = 0;
  visited[max_c_index] = 1;
  number_on = 1;
  on_indices[0] = max_c_index;
  
  /* Follow solution until l2 error is less than epsilon. */
  iters = 0;
  while(((cur_error-epsilon)>PRECISION)&&(number_on<max_non_zero)&&(iters<max_iters)&&(!cholesky_failure)){

    if(VERYVERBOSE){
      printf("iteration %d start (%d %.14f %.14f)\n",iters,max_c_index,lambda,l2Error(x_vec));
    }

    computeD();
    computeG1G2();
    computeG3();
    update();

    cur_error = l2Error(x_vec);

    /* Report current state. */
    if(VERYVERBOSE){
      printf(" iteration %d end (%d, %.14f, %.14f)\n",iters,number_on,lambda,l2Error(x_vec));

      printf("  On set: <index> <on_index> <x> <d_vec>\n");
      for(i=0;i<number_on;i++){
	printf("   %d %d %.14f %.14f\n",i,on_indices[i],x_vec[on_indices[i]],d_vec[i]);
      }
      printf("\n");
    }

    if((cur_error-epsilon)>PRECISION){
      updateOnIndices();
      computeC();
    }
    
    iters++;
  }

  /* Print warning message in the event of too many non-zero elements. */
  if(number_on>=max_non_zero){
    if(WARNINGS){
      printf("Max non zero value reached, solution is not optimal.\n");
      printf("  %d elements in the on set\n", number_on);
      printf("  %d iterations performed\n", iters);
      printf("  current lambda is %f\n", lambda);
      printf("  current error is %f (%f target)\n\n", cur_error, epsilon);
    }
    is_bad = 1;
  }

  /* Print warning message in the event that the 
     maximum number of iterations was reached. */
  if(iters==max_iters){
    if(WARNINGS){
      printf("Maximum number of iterations reached, solution is not optimal.\n");
      printf("  %d elements in the on set\n", number_on);
      printf("  %d iterations performed\n", iters);
      printf("  current lambda is %f\n", lambda);
      printf("  current error is %f (%f target)\n\n", cur_error, epsilon);
    }
    is_bad = 1;
  }

  /* Print cholesky failure warning */
  if(cholesky_failure){
    if(WARNINGS){
      printf(" Cholesky failure\n");
    }
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
      for(i=0;i<number_on;i++){
	j = on_indices[i];
	work2[j] = x_vec[j] - d_vec[i]*mid;
      }
      if(l2Error(work2)>epsilon){
	high = mid;
      }
      else{
	low = mid;
      }
    }

    /* Updata x vector and calculate corresponding lambda value. */
    for(i=0;i<number_on;i++){
      j = on_indices[i];
      x_vec[j] = work2[j];
    }
    computeC();
  }
  else{
    number_on = 0;
    lambda = 0.0;
  }

  /* Print solution information. */
  if(VERBOSE){
    if(is_bad){
      printf(" Failed to find solution\n\n");
    }
    else{
      printf(" (General) Solved in %d iterations, %d non-zero, lambda = %f, l2_error = %f\n", iters, number_on, lambda, l2Error(x_vec));
      printf(" On set: <index> <on_index> <x>\n");
      for(i=0;i<number_on;i++){
	printf("  %d %d %f\n",i,on_indices[i],x_vec[on_indices[i]]);
      }
      printf("\n");
    }
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

  index = 0;
  which_g = 0;
  gamma_s = MAXGAMMA;

  /* find minimum value of g1,g2 */
  if(positive_only){
    for(i=0;i<ncols;i++){
      if(am_off[i]){
	if((g1_vec[i]>0.0)&&(g1_vec[i]<gamma_s)){
	  index = i;
	  which_g = 1;
	  gamma_s = g1_vec[i];
	}
      }
    }
  }
  else{
    for(i=0;i<ncols;i++){
      if(am_off[i]){
	if((g1_vec[i]>0.0)&&(g1_vec[i]<gamma_s)){
	  index = i;
	  which_g = 1;
	  gamma_s = g1_vec[i];
	}
	if((g2_vec[i]>0.0)&&(g2_vec[i]<gamma_s)){
	  index = i;
	  which_g = 2;
	  gamma_s = g2_vec[i];
	}
      }
    }
  }

  /* find minimum value of g3 */
  for(i=0;i<number_on;i++){
    if((g3_vec[i]>0.0)&&(g3_vec[i]<gamma_s)){
      index = on_indices[i];
      which_g = 3;
      gamma_s = g3_vec[i];
    }
  }

  /* update x vector */
  for(i=0;i<number_on;i++){
    x_vec[on_indices[i]] += gamma_s*d_vec[i];
  }

  /* update am_off */
  if(which_g==3){
    if(VERYVERBOSE){
      printf(" Removing: %d %.14f\n",index,gamma_s);
    }
    am_off[index] = 1;
  }
  else{
    if(VERYVERBOSE){
      printf(" Adding: %d %.14f\n",index,gamma_s);
    }
    am_off[index] = 0;
    visited[index] = 1;
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

  number_on = 0;
  for(i=0;i<ncols;i++){
    if((!am_off[i])&&(number_on<max_non_zero)){
      on_indices[number_on] = i;
      number_on++;
    }
  }
}

/*
 * The MIT License
 *
 * Copyright (c) 2012 Zhuang Lab, Harvard University
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
