/*
 * C library for (F)ISTA calculations.
 *
 * Hazen 07/13
 *
 * Compilation instructions:
 *
 * Windows:
 *  gcc -c fista_lib.c
 *  gcc -shared -o fista_lib.dll fista_lib.o
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <windows.h>

/* Define */
#define ISZERO 1.0
#define MAXITERS 100000

/* Structures */

/* Function Declarations */
void cleanup(void);
void computePlK(double *, double, double);
__int64 getClock(void);
void getXVector(double *);
void initialize(double *, int, int, int);
void iterateFISTA(double, double, int);
int iterateFISTAToL0Target(double, double, int);
void iterateISTA(double, double, int);
int l0Norm(void);
double l1Error(void);
double l2Error(void);
void newBVector(double *);
void printProfilingData(void);

/* Global Variables */
static int ncols;
static int nrows;
static int positive_only;

static __int64 solver_time;
static double clock_freq;

static double tk;

static double *a_mat;
static double *a_mat_trans;
static double *a_mat_x;
static double *b_vec;
static double *work1;
static double *work2;
static double *x_vec;
static double *x_vec_old;
static double *y_vec;


/*
 * cleanup()
 *
 * Frees all the allocated storage.
 */
void cleanup(void)
{
  free(a_mat);
  free(a_mat_trans);
  free(a_mat_x);
  free(b_vec);
  free(work1);
  free(work2);
  free(x_vec);
  free(x_vec_old);
  free(y_vec);
}

/*
 * computePlK()
 *
 * Performs a single cycle of update & shrinkage, stores
 * the result in work1.
 *
 * a_x_vec - x vector to use.
 * lambda - lambda.
 * time_step - time step.
 */
void computePlK(double *a_x_vec, double lambda, double step_size)
{
  int i,j,k;
  double lt,sum,t1,t2;

  lt = lambda*step_size;
  t1 = 2.0*step_size;

  /* Compute Ax. */
  for(i=0;i<nrows;i++){
    j = i*ncols;
    a_mat_x[i] = 0.0;
    for(k=0;k<ncols;k++){
      a_mat_x[i] += a_mat[j+k]*a_x_vec[k];
    }
  }

  /* Compute Ax - b. */
  for(i=0;i<nrows;i++){
    work2[i] = a_mat_x[i] - b_vec[i];
  }

  /* Compute work1. */
  for(i=0;i<ncols;i++){
    j = i*nrows;
    sum = 0.0;
    for(k=0;k<nrows;k++){
      sum += a_mat_trans[j+k] * work2[k];
    }
    work1[i] = a_x_vec[i] - t1*sum;
  }

  /* Shrink work1. */
  if(positive_only){
    for(i=0;i<ncols;i++){
      t2 = work1[i] - lt;
      if(t2 < 0.0){
	work1[i] = 0.0;
      }
      else{
	work1[i] = t2;
      }
    }
  }
  else{
    for(i=0;i<ncols;i++){
      t2 = fabs(work1[i]) - lt;
      if(t2<0.0){
	work1[i] = 0.0;
      }
      else{
	work1[i] = (work1[i]>0.0) ? t2 : -t2;
      }
    }
  }
}

/*
 * getClock()
 *
 * Returns the current system clock time.
 */
__int64 getClock(void)
{
  LARGE_INTEGER li;

  QueryPerformanceCounter(&li);
  return ((__int64)li.QuadPart);
}

/*
 * getXVector()
 *
 * Return the current x vector in user supplied storage.
 *
 * xvec - pre-allocated space for the current x vector.
 */
void getXVector(double *xvec)
{
  int i;

  for(i=0;i<ncols;i++){
    xvec[i] = x_vec[i];
  }
}

/*
 * initialize()
 *
 * Initialize (F)ISTA solver.
 *
 * A - the A matrix.
 * rows - the number of rows in the A matrix.
 * cols - the number of columns in the A matrix.
 * pos_only - solutions are positive only.
 */
void initialize(double *A, int rows, int cols, int pos_only)
{
  int i,j;
  LARGE_INTEGER li;

  /* Set global variables. */
  nrows = rows;
  ncols = cols;
  positive_only = pos_only;

  /* Allocate storage. */
  a_mat = (double *)malloc(sizeof(double)*nrows*ncols);
  a_mat_trans = (double *)malloc(sizeof(double)*nrows*ncols);
  a_mat_x = (double *)malloc(sizeof(double)*nrows);
  b_vec = (double *)malloc(sizeof(double)*nrows);
  work1 = (double *)malloc(sizeof(double)*ncols);
  work2 = (double *)malloc(sizeof(double)*nrows);
  x_vec = (double *)malloc(sizeof(double)*ncols);
  x_vec_old = (double *)malloc(sizeof(double)*ncols);
  y_vec = (double *)malloc(sizeof(double)*ncols);

  /* Copy A matrix into local storage. */
  for(i=0;i<nrows;i++){
    for(j=0;j<ncols;j++){
      a_mat[i*ncols+j] = A[i*ncols+j];
      a_mat_trans[j*nrows+i] = A[i*ncols+j];
    }
  }

  /* Profiling code. */
  QueryPerformanceFrequency(&li);
  clock_freq = ((double)li.QuadPart);
 
  solver_time = 0;
}

/*
 * iterateFISTA()
 *
 * Perform requested number of steps of FISTA optimization.
 *
 * lambda - lambda value to use.
 * step_size - step size to use.
 * iters - number of iterations to perform.
 */
void iterateFISTA(double lambda, double step_size, int iters)
{
  int i,j;
  double new_tk,t1;

  for(i=0;i<iters;i++){

    for(j=0;j<ncols;j++){
      x_vec_old[j] = x_vec[j];
    }

    computePlK(y_vec, lambda, step_size);
    for(j=0;j<ncols;j++){
      x_vec[j] = work1[j];
    }

    new_tk = 0.5*(1.0 + sqrt(1.0 + 4.0*tk*tk));
    t1 = (tk - 1.0)/new_tk;

    for(j=0;j<ncols;j++){
      y_vec[j] = x_vec[j] + t1*(x_vec[j] - x_vec_old[j]);
    }

    tk = new_tk;
  }
}

/*
 * iterateFISTAToL0Target()
 *
 * Perform FISTA optimization until l0 target is reached.
 *
 * lambda - lambda value to use.
 * step_size - step size to use.
 * l0_target - number of iterations to perform.
 */
int iterateFISTAToL0Target(double lambda, double step_size, int l0_target)
{
  int cur_l0,iters,total;
  __int64 start;

  start = getClock();

  iters = 1000;
  total = 0;
  cur_l0 = l0_target + 1;
  while((cur_l0 > l0_target)&&(total<MAXITERS)){
    iterateFISTA(lambda, step_size, iters);
    total += iters;
    cur_l0 = l0Norm();
  }

  solver_time = getClock() - start;

  return total;
}

/*
 * iterateISTA()
 *
 * Perform requested number of steps of ISTA optimization.
 *
 * lambda - lambda value to use.
 * step_size - step size to use.
 * iters - number of iterations to perform.
 */
void iterateISTA(double lambda, double step_size, int iters)
{
  int i,j;

  for(i=0;i<iters;i++){
    computePlK(x_vec, lambda, step_size);
    for(j=0;j<ncols;j++){
      x_vec[j] = work1[j];
    }
  }
}

/*
 * l0Norm()
 *
 * Return the l0norm of the current solution.
 */
int l0Norm(void)
{
  int cnts,i;

  /*
  double max;

  max = 0.0;
  for(i=0;i<ncols;i++){
    if(fabs(x_vec[i]) > max){
      max = x_vec[i];
    }
  }
  max = ISZERO * max;
  */

  cnts = 0;
  for(i=0;i<ncols;i++){
    if(fabs(x_vec[i]) > ISZERO){
      cnts++;
    }
  }

  printf("  l0norm %d\n", cnts);
  return cnts;
}

/*
 * l1Error()
 *
 * Return the l1Error of the current solution.
 */
double l1Error(void)
{
  int i;
  double sum;
  
  sum = 0.0;
  for(i=0;i<ncols;i++){
    sum += fabs(x_vec[i]);
  }

  return sum;
}

/*
 * l2Error()
 *
 * Return the (square of the) l2Error of the current solution.
 */
double l2Error(void)
{
  int i,j,k;
  double sum,temp;

  /* Compute Ax. */
  for(i=0;i<nrows;i++){
    j = i*ncols;
    a_mat_x[i] = 0.0;
    for(k=0;k<ncols;k++){
      a_mat_x[i] += a_mat[j+k]*x_vec[k];
    }
  }
  
  /* Compute l2 error. */
  sum = 0.0;
  for(i=0;i<nrows;i++){
    temp = a_mat_x[i] - y_vec[i];
    sum += temp*temp;
  }

  return sum;
}

/*
 * newBVector()
 *
 * Reset solver for a new b vector.
 *
 * bvec - new b vector.
 */
void newBVector(double *bvec)
{
  int i;

  /* Copy b vector, zero Ax */
  for(i=0;i<nrows;i++){
    b_vec[i] = bvec[i];
    a_mat_x[i] = 0.0;
  }

  /* Reset solver. */
  tk = 1.0;
  for(i=0;i<ncols;i++){
    x_vec[i] = 0.0;
    x_vec_old[i] = 0.0;
    y_vec[i] = 0.0;
  }
}

/*
 * printProfilingData()
 *
 * Prints out the profiling information.
 */
void printProfilingData(void)
{
  printf("FISTA profiling data:\n");
  if (0){
    printf(" solver: %lld ticks\n", solver_time);
  }
  else{
    printf(" solver: %.4f seconds\n", ((double)solver_time)/clock_freq);
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
