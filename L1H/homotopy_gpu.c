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
 * Hazen 8/12
 *
 * This version has been cleaned up a bit and optimized for analysis of STORM images.
 * The assumptions here are:
 *   1. x is positive only.
 *   2. The last term of x is assumed to be the background term.
 *   3. The computation of G1 is shortcut based on values in C.
 *
 * Hazen 1/13
 *
 * This version uses a GPU to do some of the calculations. Since the math is done in
 * floating point rather than double you may get slightly different answers..
 *
 * Hazen 1/13
 *
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <CL/opencl.h>

#include "homotopy_gpu.h"

/* Define */
#define DEBUG 0
#define MAXGAMMA 1.0e+6
#define MINCFACTOR 0.5
#define PRECISION 1.0e-6
#define USE_GPU 0
#define VERBOSE 0
#define VERYVERBOSE 0
#define WARNINGS 1

#define ONSET 0
#define OFFSET 1

#define NTIMING 7

/* Structures */

/* Function Declarations */
void computeC(void);
void computeD(void);
void computeG1(void);
void computeG3(void);
size_t gpuCalcGlobalSize(int);
void gpuErrorCheck(int, const char *);
float l2Error(float *);
void printProfilingData(void);
void update(void);
void updateOnIndices(void);
void updateProfilingData(int);

/* LAPACK Functions */
/* Cholesky solver */
extern void dposv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
		   double* b, int* ldb, int* info);

/*
 * Global Variables 
 */

/* homotopy related */
static int cholesky_failure;
static int max_c_index;
static int max_non_zero;
static int ncols;
static int nrows;
static int n_onset;
static int which_g1;
static int which_g3;
static int work_size;

static int *offset;
static int *onset_indices;

static float gamma_g1;
static float gamma_g3;
static float gamma_s;
static float lambda;

static float *a_mat;
static float *a_mat_trans;
static float *a_y_vec;
static float *c_vec;
static float *d_vec;
static float *g_mat;
static float *g1_vec;
static float *x_vec;
static float *y_vec;

static float *work1;
static float *work2;

static double *d_d_vec;
static double *d_work1;

/* GPU related */
int gpu_err_num;

size_t gpu_global_work_size;
size_t gpu_local_work_size;

cl_command_queue gpu_command_queue;
cl_kernel gpu_compute_c_vec_kernel;
cl_kernel gpu_compute_g1_vec_kernel;
cl_kernel gpu_compute_work1_kernel;
cl_context gpu_context;
cl_device_id* gpu_device;
cl_event gpu_event;
cl_platform_id* gpu_platform;
cl_program gpu_program;

cl_mem in_gpu_a_matrix;
cl_mem in_gpu_a_matrix_trans;
cl_mem in_gpu_a_y_vec;
cl_mem in_gpu_c_vec;
cl_mem in_gpu_d_vec;
cl_mem in_gpu_g_matrix;
cl_mem in_gpu_on_indices;
cl_mem in_gpu_work1_vec;
cl_mem in_gpu_x_vec;

cl_mem out_gpu_c_vec;
cl_mem out_gpu_g1_vec;
cl_mem out_gpu_work1_vec;

cl_ulong *timing_data;

/*
 * cleanup()
 *
 * Frees all the allocated storage.
 */
void cleanup(void)
{
  free(offset);
  free(onset_indices);

  free(a_mat);
  free(a_mat_trans);
  free(a_y_vec);
  free(c_vec);
  free(d_vec);
  free(g_mat);
  free(g1_vec);
  free(x_vec);
  free(y_vec);

  free(work1);
  free(work2);

  free(d_d_vec);
  free(d_work1);

  clReleaseCommandQueue(gpu_command_queue);
  clReleaseKernel(gpu_compute_c_vec_kernel);
  clReleaseKernel(gpu_compute_g1_vec_kernel);
  clReleaseKernel(gpu_compute_work1_kernel);
  clReleaseContext(gpu_context);
  free(gpu_device);
  clReleaseEvent(gpu_event);
  free(gpu_platform);
  clReleaseProgram(gpu_program);

  clReleaseMemObject(in_gpu_a_matrix);
  clReleaseMemObject(in_gpu_a_matrix_trans);
  clReleaseMemObject(in_gpu_a_y_vec);
  clReleaseMemObject(in_gpu_c_vec);
  clReleaseMemObject(in_gpu_d_vec);
  clReleaseMemObject(in_gpu_g_matrix);
  clReleaseMemObject(in_gpu_on_indices);
  clReleaseMemObject(in_gpu_work1_vec);
  clReleaseMemObject(in_gpu_x_vec);

  clReleaseMemObject(out_gpu_c_vec);
  clReleaseMemObject(out_gpu_g1_vec);
  clReleaseMemObject(out_gpu_work1_vec);

  free(timing_data);
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

  if (DEBUG){
    printf("(computeC) x_vec: ");
    for(i=0;i<ncols;i++){
      printf("%f ", x_vec[i]);
    }
    printf("\n");
  }

  if (USE_GPU){
    /* Update vectors & variables that have changed. */
    gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_x_vec, CL_TRUE, 0, ncols*sizeof(float), x_vec, 0, NULL, &gpu_event),
		  "clEnqueueWriteBuffer");
    updateProfilingData(0);
    gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_on_indices, CL_TRUE, 0, n_onset*sizeof(int), onset_indices, 0, NULL, &gpu_event),
		"clEnqueueWriteBuffer");
    updateProfilingData(0);
    gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 5, sizeof(cl_int), (void*) &n_onset), "clSetKernelArg");

    /* Calculate new c vector. */
    gpu_global_work_size = gpuCalcGlobalSize(ncols);
    gpuErrorCheck(clEnqueueNDRangeKernel(gpu_command_queue, gpu_compute_c_vec_kernel, 1, NULL, &gpu_global_work_size, &gpu_local_work_size, 0, NULL, &gpu_event), 
		  "clEnqueueNDRangeKernel");
    updateProfilingData(1);

    /* Get new c vector. */
    gpuErrorCheck(clEnqueueReadBuffer(gpu_command_queue, out_gpu_c_vec, CL_TRUE, 0, ncols*sizeof(float), c_vec, 0, NULL, &gpu_event),
		  "clEnqueueReadBuffer");
    updateProfilingData(0);

    /* Find the maximum of the c vector & its index. */
    lambda = 0.0;
    max_c_index = 0;
    for(i=0;i<ncols;i++){
      if(c_vec[i]>lambda){
	lambda = c_vec[i];
	max_c_index = i;
      }
    }
  }
  else{
    for(i=0;i<ncols;i++){
      work1[i] = 0.0;
    }
  
    for(i=0;i<n_onset;i++){
      j = onset_indices[i];
      k = j*ncols;
      for(l=0;l<ncols;l++){
	// using the fact that g_mat is symmetric (this goes across a row).
	work1[l] += g_mat[k+l]*x_vec[j];
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
  }

  if (DEBUG){
    printf("(computeC) onset size: %d\n", n_onset);
    for(i=0;i<ncols;i++){
      printf("%f ", c_vec[i]);
    }
    printf("\n");
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
  float tmp;

  /* 
   * Compute A'[S] * A[S].
   * Since we already calculate A' * A, we just need to
   * pull the appropriate elements out of the G matrix.
   */
  for(i=0;i<n_onset;i++){
    for(j=0;j<n_onset;j++){
      d_work1[i*n_onset+j] = (double)g_mat[onset_indices[i]*ncols+onset_indices[j]];
    }
  }

  /* compute sign of c vector */
  for(i=0;i<n_onset;i++){
    tmp = c_vec[onset_indices[i]];
    if(tmp>0.0){
      d_d_vec[i] = 1.0;
    }
    else if(tmp<0.0){
      d_d_vec[i] = -1.0;
    }
    else{
      d_d_vec[i] = 0.0;
    }
  }

  if(n_onset==1){
    d_d_vec[0] = d_d_vec[0]/d_work1[0];
  }
  else{
    /* use LAPACK cholesky solver to determine d */
    nrhs = 1;
    n = lda = ldb = n_onset;
    dposv_( "Lower", &n, &nrhs, d_work1, &lda, d_d_vec, &ldb, &info );
    if((info!=0)&&(WARNINGS)){
      printf("  Cholesky solver failed with error: %d\n", info);
      cholesky_failure = 1;
    }
  }

  for(i=0;i<n_onset;i++){
    d_vec[i] = (float)d_d_vec[i];
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

}

/*
 * computeG1()
 *
 * Compute gamma1 at current lambda value.
 */
void computeG1(void)
{
  int i,j,k,l;
  double w,tmp;

  /* 
   * Compute v 
   */

  if (USE_GPU){
    /* Update vectors & variables that have changed. */
    gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_d_vec, CL_TRUE, 0, n_onset*sizeof(float), d_vec, 0, NULL, &gpu_event),
		  "clEnqueueWriteBuffer");
    updateProfilingData(2);
    gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_on_indices, CL_TRUE, 0, n_onset*sizeof(int), onset_indices, 0, NULL, &gpu_event),
		  "clEnqueueWriteBuffer");
    updateProfilingData(2);
    gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 5, sizeof(cl_int), (void*) &n_onset), "clSetKernelArg");

    /* Calculate new work1 vector. */
    gpu_global_work_size = gpuCalcGlobalSize(nrows);
    gpuErrorCheck(clEnqueueNDRangeKernel(gpu_command_queue, gpu_compute_work1_kernel, 1, NULL, &gpu_global_work_size, &gpu_local_work_size, 0, NULL, &gpu_event),
		  "clEnqueueNDRangeKernel");
    updateProfilingData(3);

    /* Get new work1 vector. */
    gpuErrorCheck(clEnqueueReadBuffer(gpu_command_queue, out_gpu_work1_vec, CL_TRUE, 0, nrows*sizeof(float), work1, 0, NULL, &gpu_event),
		  "clEnqueueReadBuffer");
    updateProfilingData(2);
  }
  else{  
    for(i=0;i<nrows;i++){
      work1[i] = 0.0;
      j = i*ncols;
      for(k=0;k<n_onset;k++){
	l = onset_indices[k];
	work1[i] += a_mat[j+l]*d_vec[k];
      }
    }
  }

  if (DEBUG){
    printf("(computeG1) work1: ");
    for(i=0;i<nrows;i++){
      printf("%f ", work1[i]);
    }
    printf("\n");
  }

  /* 
   * Compute g1 
   */

  if (USE_GPU){
    /* Update vectors & variables that have changed. */
    gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_c_vec, CL_TRUE, 0, ncols*sizeof(float), c_vec, 0, NULL, &gpu_event),
		  "clEnqueueWriteBuffer");
    updateProfilingData(4);
    gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_work1_vec, CL_TRUE, 0, nrows*sizeof(int), work1, 0, NULL, &gpu_event),
		  "clEnqueueWriteBuffer");
    updateProfilingData(4);
    gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 5, sizeof(cl_float), (void*) &lambda), "clSetKernelArg");

    /* Calculate new work1 vector. */
    gpu_global_work_size = gpuCalcGlobalSize(ncols);
    gpuErrorCheck(clEnqueueNDRangeKernel(gpu_command_queue, gpu_compute_g1_vec_kernel, 1, NULL, &gpu_global_work_size, &gpu_local_work_size, 0, NULL, &gpu_event),
		  "clEnqueueNDRangeKernel");
    updateProfilingData(5);

    /* Get new g1 vector. */
    gpuErrorCheck(clEnqueueReadBuffer(gpu_command_queue, out_gpu_g1_vec, CL_TRUE, 0, ncols*sizeof(float), g1_vec, 0, NULL, &gpu_event),
		  "clEnqueueReadBuffer");
    updateProfilingData(4);
    
    /* Find the maximum of the g1 vector & its index. */
    gamma_g1 = MAXGAMMA;
    which_g1 = 0;
    for(i=0;i<ncols;i++){
      if(offset[i]==OFFSET){
	if((g1_vec[i]>0.0)&&(g1_vec[i]<gamma_g1)){
	  gamma_g1 = g1_vec[i];
	  which_g1 = i;
	}
      }
    }
  }

  else{
    gamma_g1 = MAXGAMMA;
    for(i=0;i<ncols;i++){
      if((offset[i]==OFFSET)&&(c_vec[i]>0.0)){
	w = 0.0;
	for(j=0;j<nrows;j++){
	  w += a_mat[j*ncols+i]*work1[j];
	}
	tmp = (lambda - c_vec[i])/(1.0 - w);
	if((tmp>0.0)&&(tmp<gamma_g1)){
	  gamma_g1 = tmp;
	  which_g1 = i;
	}
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
  int i,j;
  double tmp;

  gamma_g3 = MAXGAMMA;
  for(i=0;i<(n_onset-1);i++){
    j = onset_indices[i];
    tmp = -1.0*x_vec[j]/d_vec[i];
    if((tmp>0.0)&&(tmp<gamma_g3)){
      gamma_g3 = tmp;
      which_g3 = j;
    }
  }
}

/*
 * getXVector()
 *
 * Returns the current x vector in user supplied storage.
 *
 * xvec - user supplied storage for the x vector, assumed initialized to zero.
 */
void getXVector(float *xvec)
{
  int i,j;

  for(i=0;i<n_onset;i++){
    j = onset_indices[i];
    xvec[j] = x_vec[j];
  }
}

/*
 * gpuCalcGlobalSize()
 *
 * Figure out correct global size given problem size.
 *
 * problem_size - number of elements in the problem.
 */
size_t gpuCalcGlobalSize(int problem_size)
{
  size_t tmp;

  tmp = ((size_t)problem_size)/gpu_local_work_size;
  tmp = gpu_local_work_size*(tmp+1);
  
  return tmp;
}

/*
 * gpuErrorCheck()
 *
 * Check the return code of GPU function calls for indications of errors.
 *
 * error_code - the error code.
 * fn_name - the name of the function that was called.
 */
void gpuErrorCheck(int error_code, const char *fn_name)
{
  if(error_code){
    printf("OpenCL error: %s - %d\n", fn_name, error_code);
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
void initialize(float *A, int rows, int cols, int bgterm, int pos_only, int number_nonzero)
{
  int i,j,k,l;
  double sum;

  max_non_zero = number_nonzero;
  ncols = cols;
  nrows = rows;

  /* Allocate memory. */
  offset = (int *)malloc(sizeof(int)*ncols);
  onset_indices = (int *)malloc(sizeof(int)*(max_non_zero+1));

  a_mat = (float *)malloc(sizeof(float)*nrows*ncols);
  a_mat_trans = (float *)malloc(sizeof(float)*nrows*ncols);
  a_y_vec = (float *)malloc(sizeof(float)*ncols);
  c_vec = (float *)malloc(sizeof(float)*ncols);
  d_vec = (float *)malloc(sizeof(float)*max_non_zero);
  g_mat = (float *)malloc(sizeof(float)*ncols*ncols);
  g1_vec = (float *)malloc(sizeof(float)*ncols);
  x_vec = (float *)malloc(sizeof(float)*ncols);
  y_vec = (float *)malloc(sizeof(float)*nrows);

  d_d_vec = (double *)malloc(sizeof(double)*max_non_zero);

  if((max_non_zero*max_non_zero)>ncols){
    work_size = max_non_zero*max_non_zero;
  }
  else{
    work_size = ncols;
  }
  work1 = (float *)malloc(sizeof(float)*work_size);
  work2 = (float *)malloc(sizeof(float)*work_size);
  d_work1 = (double *)malloc(sizeof(double)*work_size);

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

  /* 
   * GPU memory initialization. 
   */

  /* Allocate OpenCL buffer memory objects. */

  /* Input */
  in_gpu_a_matrix = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, nrows*ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_a_matrix_trans = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, nrows*ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_a_y_vec = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_c_vec = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_d_vec = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, max_non_zero*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_g_matrix = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, ncols*ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_on_indices = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, (max_non_zero+1)*sizeof(int), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_work1_vec = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, work_size*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  in_gpu_x_vec = clCreateBuffer(gpu_context, CL_MEM_READ_ONLY, ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  /* Output */
  out_gpu_c_vec = clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY, ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  out_gpu_g1_vec = clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY, ncols*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  out_gpu_work1_vec = clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY, work_size*sizeof(float), NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateBuffer");

  /* Fill memory objects. */
  gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_a_matrix, CL_FALSE, 0, nrows*ncols*sizeof(float), a_mat, 0, NULL, NULL),
		"clEnqueueWriteBuffer");
  gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_a_matrix_trans, CL_FALSE, 0, nrows*ncols*sizeof(float), a_mat_trans, 0, NULL, NULL),
		"clEnqueueWriteBuffer");
  gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_g_matrix, CL_FALSE, 0, ncols*ncols*sizeof(float), g_mat, 0, NULL, NULL),
		"clEnqueueWriteBuffer");

  /* Set kernel argument values. */
  gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 0, sizeof(cl_mem), (void*) &in_gpu_g_matrix), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 1, sizeof(cl_mem), (void*) &in_gpu_a_y_vec), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 2, sizeof(cl_mem), (void*) &in_gpu_x_vec), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 3, sizeof(cl_mem), (void*) &in_gpu_on_indices), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 4, sizeof(cl_int), (void*) &ncols), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 5, sizeof(cl_int), (void*) &n_onset), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_c_vec_kernel, 6, sizeof(cl_mem), (void*) &out_gpu_c_vec), "clSetKernelArg");

  gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 0, sizeof(cl_mem), (void*) &in_gpu_a_matrix), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 1, sizeof(cl_mem), (void*) &in_gpu_d_vec), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 2, sizeof(cl_mem), (void*) &in_gpu_on_indices), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 3, sizeof(cl_int), (void*) &nrows), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 4, sizeof(cl_int), (void*) &ncols), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 5, sizeof(cl_int), (void*) &n_onset), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_work1_kernel, 6, sizeof(cl_mem), (void*) &out_gpu_work1_vec), "clSetKernelArg");

  gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 0, sizeof(cl_mem), (void*) &in_gpu_a_matrix_trans), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 1, sizeof(cl_mem), (void*) &in_gpu_c_vec), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 2, sizeof(cl_mem), (void*) &in_gpu_work1_vec), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 3, sizeof(cl_int), (void*) &nrows), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 4, sizeof(cl_int), (void*) &ncols), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 5, sizeof(cl_float), (void*) &lambda), "clSetKernelArg");
  gpuErrorCheck(clSetKernelArg(gpu_compute_g1_vec_kernel, 6, sizeof(cl_mem), (void*) &out_gpu_g1_vec), "clSetKernelArg");

  /* Allocate space for timing data storage. */
  timing_data = (cl_ulong *)malloc(NTIMING*sizeof(cl_ulong));

  for(i=0;i<NTIMING;i++){
    timing_data[i] = 0;
  }
}

/*
 * initializeGPU()
 *
 * Perform initial setup for GPU based calculations. The
 * program uses the first device on the specified platform.
 *
 * gpu_code - cl code to be compiled into GPU programs.
 * which_platform - the number of the GPU platform to use.
 * which_device - the number of the GPU device to use.
 * local_work_size - sub problem size for GPU calculations?
 *
 * FIXME: which_platform can only be 0??
 */
void initializeGPU(char *gpu_code, int which_platform, int which_device, int local_work_size)
{
  char *log;
  cl_uint num_devices;
  int error_code;
  size_t gpu_code_size, log_size;

  num_devices = 1;
  gpu_local_work_size = (size_t)local_work_size;

  /* Connect to GPU. */
  gpu_platform = (cl_platform_id *)malloc(sizeof(cl_platform_id));
  gpu_device = (cl_device_id *)malloc(sizeof(cl_device_id));
  gpuErrorCheck(clGetPlatformIDs(1, gpu_platform, NULL),
		"clGetPlatformIDs");
  gpuErrorCheck(clGetDeviceIDs(gpu_platform[which_platform], CL_DEVICE_TYPE_GPU, num_devices, gpu_device, &num_devices),
		"clGetDeviceIDs");
  gpu_context = clCreateContext(0, 1, gpu_device, NULL, NULL, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateContext");
  gpu_command_queue = clCreateCommandQueue(gpu_context, gpu_device[which_device], CL_QUEUE_PROFILING_ENABLE, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateCommandQueue");

  /* Compile GPU code. */
  gpu_code_size = strlen(gpu_code);
  gpu_program = clCreateProgramWithSource(gpu_context, 1, (const char **)&gpu_code, &gpu_code_size, &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateProgramWithSource");
  error_code = clBuildProgram(gpu_program, 1, gpu_device, "-cl-fast-relaxed-math", NULL, NULL);
  gpuErrorCheck(error_code, "clBuildProgram");
  if (error_code == CL_BUILD_PROGRAM_FAILURE){
    gpuErrorCheck(clGetProgramBuildInfo(gpu_program, gpu_device[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size), "clGetProgramBuildInfo");
    log = (char *) malloc(log_size+1);
    gpuErrorCheck(clGetProgramBuildInfo(gpu_program, gpu_device[0], CL_PROGRAM_BUILD_LOG, log_size, log, NULL), "clGetProgramBuildInfo");
    log[log_size] = '\0';
    printf("------- Begin GPU Build Error ----------\n");
    printf("%s\n", log); 
    printf("------- End GPU Build Error ----------\n");
  }

  /* Create kernels */
  gpu_compute_c_vec_kernel = clCreateKernel(gpu_program, "computeCVec", &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateKernel");
  gpu_compute_work1_kernel = clCreateKernel(gpu_program, "computeWork1", &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateKernel");
  gpu_compute_g1_vec_kernel = clCreateKernel(gpu_program, "computeG1Vec", &gpu_err_num);
  gpuErrorCheck(gpu_err_num, "clCreateKernel");

}

/*
 * l2Error()
 *
 * Calculate the l2Error given a x vector.
 *
 * xvec - a x_vector
 */
float l2Error(float *xvec)
{
  int i,j,l;
  float sum,tmp;

  /* Calculate Ax */
  for(i=0;i<nrows;i++){
    work1[i] = 0.0;
  }

  for(i=0;i<n_onset;i++){
    l = onset_indices[i];
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
 * Setup problem with a new y vector. This does the following:
 *  1. Copy the y vector into local storage.
 *  2. Compute the product Ay and save this (this is also the initial C vector).
 *  3. Compute the maximum element of Ay and its location.
 *  4. Initialize offset vectors.
 *
 * yvec - the new y vector (assumed of size nrows).
 */
void newYVector(float *yvec)
{
  int i,j;
  float sum;

  /* 1. Copy into y_vec. */
  for(i=0;i<nrows;i++){
    y_vec[i] = yvec[i];
  }

  /* 2/3. Calculate a_y_vec, c_vec, find maximum. */
  lambda = 0.0;
  max_c_index = 0;
  for(i=0;i<ncols;i++){
    sum = 0.0;
    for(j=0;j<nrows;j++){
      sum += a_mat[j*ncols+i]*y_vec[j];
    }
    a_y_vec[i] = sum;
    c_vec[i] = sum;
    if(sum>lambda){
      lambda = sum;
      max_c_index = i;
    }
  }

  /* Update a_y_vec in GPU memory. */
  if (USE_GPU){
    gpuErrorCheck(clEnqueueWriteBuffer(gpu_command_queue, in_gpu_a_y_vec, CL_FALSE, 0, ncols*sizeof(float), a_y_vec, 0, NULL, &gpu_event),
		  "clEnqueueWriteBuffer");
    updateProfilingData(6);
  }

  /* 4. Initialize offset vector. */
  for(i=0;i<ncols;i++){
    offset[i] = OFFSET;
  }
  offset[max_c_index] = ONSET;

}

/*
 * printProfilingData()
 *
 * Print GPU profiling information.
 */
void printProfilingData(void)
{
  int i;
  cl_uint sum;

  sum = 0;
  printf("\nGPU profiling results (nanoseconds):\n");
  for(i=0;i<NTIMING;i++){
    sum += timing_data[i];
    printf(" %ld) %d (%.f)\n", i, timing_data[i], (1.0e-6 * (double)timing_data[i]));
  }
  printf(" Total: %ld (%.f)\n", sum, (1.0e-6 * (double)sum));
  printf("\n");
}

/*
 * solve()
 *
 * Follow solution path until l2 error is equal to epsilon.
 *
 * epsilon - desired l2 error.
 * max_iters - maximum number of iterations before giving up.
 */
float solve(float epsilon, int max_iters)
{
  int i,is_bad,iters,j;
  float cur_error,low,high,mid;
  double temp;

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
      printf("iteration %d start (%d %.14f %.14f)\n",iters,max_c_index,lambda,l2Error(x_vec));
    }

    computeD();
    computeG1();
    computeG3();
    update();

    cur_error = l2Error(x_vec);

    /* Report current state. */
    if(VERYVERBOSE){
      printf(" iteration %d end (%d, %.14f, %.14f)\n",iters,n_onset,lambda,l2Error(x_vec));

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

  /* Print warning message in the event of too many non-zero elements. */
  if(n_onset>=max_non_zero){
    if(WARNINGS){
      printf("Max non zero value reached, solution is not optimal.\n");
      printf("  %d elements in the on set\n", n_onset);
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
      printf("  %d elements in the on set\n", n_onset);
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

    /* Updata x vector and calculate corresponding lambda value. */
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
      printf(" Removing: %d %.14f\n",index,gamma_s);
    }
    offset[index] = OFFSET;
  }
  else{
    if(VERYVERBOSE){
      printf(" Adding: %d %.14f\n",index,gamma_s);
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
 * updateProfilingData()
 *
 * Updates the profiling data.
 *
 */
void updateProfilingData(int event_type)
{
  cl_ulong end, start;

  //gpuErrorCheck(clWaitForEvents(1, &gpu_event), "clWaitForEvents");
  gpuErrorCheck(clFinish(gpu_command_queue), "clFinish");
  gpuErrorCheck(clGetEventProfilingInfo(gpu_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL), "clGetEventProfilingInfo");
  gpuErrorCheck(clGetEventProfilingInfo(gpu_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL), "clGetEventProfilingInfo");
  //gpuErrorCheck(clGetEventProfilingInfo(gpu_event, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &start, NULL), "clGetEventProfilingInfo");

  //printf(" %d %ld %ld\n", event_type, start, end);
  timing_data[event_type] += (end - start);
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
