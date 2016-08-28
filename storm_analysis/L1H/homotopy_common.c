/*
 * Functions that are common to all the homotopy solvers.
 *
 * Hazen 08/13
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#include "homotopy_common.h"

static int total_iterations;
#ifdef _WIN32
static __int64 start;
#else
static struct timespec start;
#endif
static double clock_freq;

static int *failure_counter;
static __int64 *profile_counter;

/*
 * freeCommon()
 *
 * Free the common global variables.
 */
void freeCommon(void)
{
  free(failure_counter);
  free(profile_counter);
}

/*
 * getClock()
 *
 * get profiling clock time.
 */
__int64 getClock(void)
{
  #ifdef _WIN32
  LARGE_INTEGER li;

  QueryPerformanceCounter(&li);
  return ((__int64)li.QuadPart);
  #endif

  return 0;
}

/*
 * getL1FLTSize()
 *
 * Returns the size of the L1FLT type for the purpose of 
 * being able to query the nature of L1FLT.
 */
int getL1FLTSize(void)
{
  return ((int)sizeof(L1FLT));
}

/*
 * initCommon()
 *
 * Initialize (common) global variables.
 */
void initCommon(void)
{
  int i;
  #ifdef _WIN32
  LARGE_INTEGER li;
  #endif

  failure_counter = (int *)malloc(sizeof(int)*4);
  profile_counter = (__int64 *)malloc(sizeof(__int64)*10);

  /* Failure mode logging. */
  for(i=0;i<4;i++){
    failure_counter[i] = 0;
  }

  clock_freq = 1.0e9;

  #ifdef _WIN32
  QueryPerformanceFrequency(&li);
  clock_freq = ((double)li.QuadPart);
  #endif

  for(i=0;i<10;i++){
    profile_counter[i] = 0;
  }

  total_iterations = 0;
}

/*
 * printFailureCounter()
 *
 * Prints out failure counter data.
 */
void printFailureCounter(void)
{
  printf(" Total blocks: %d, Failures: Max Non Zero %d, Max Iterations %d, Cholesky %d\n", failure_counter[0], failure_counter[1], failure_counter[2], failure_counter[3]);
}

/*
 * printProfilingData()
 *
 * Prints out the profiling information.
 */
void printProfilingData(void)
{
  int i;
  __int64 sum;

  sum = 0;
  for(i=0;i<5;i++){
    sum += profile_counter[i];
  }

  printf("Total iterations: %d\n", total_iterations);
  printf("Profiling data:\n");
  if (sum == 0){
    printf(" No profiling data\n");
  }
  else if ((((double)sum)/clock_freq) < 1.0){
    #ifdef _WIN32
    printf("   C: %lld ticks\n", profile_counter[0]);
    printf("   D: %lld ticks\n", profile_counter[1]);
    printf("  G1: %lld ticks\n", profile_counter[2]);
    printf("  G3: %lld ticks\n", profile_counter[3]);
    printf("  L2: %lld ticks\n", profile_counter[4]);
    printf("  NY: %lld ticks\n", profile_counter[5]);
    printf("   S: %lld ticks\n", profile_counter[6]);
    printf(" sum: %lld ticks\n\n", sum);
    #else
    printf("   C: %ld ticks\n", profile_counter[0]);
    printf("   D: %ld ticks\n", profile_counter[1]);
    printf("  G1: %ld ticks\n", profile_counter[2]);
    printf("  G3: %ld ticks\n", profile_counter[3]);
    printf("  L2: %ld ticks\n", profile_counter[4]);
    printf("  NY: %ld ticks\n", profile_counter[5]);
    printf("   S: %ld ticks\n", profile_counter[6]);
    printf(" sum: %ld ticks\n\n", sum);
    #endif
  }
  else {
    printf("   C: %.4f seconds\n", ((double)profile_counter[0])/clock_freq);
    printf("   D: %.4f seconds\n", ((double)profile_counter[1])/clock_freq);
    printf("  G1: %.4f seconds\n", ((double)profile_counter[2])/clock_freq);
    printf("  G3: %.4f seconds\n", ((double)profile_counter[3])/clock_freq);
    printf("  L2: %.4f seconds\n", ((double)profile_counter[4])/clock_freq);
    printf("  NY: %.4f seconds\n", ((double)profile_counter[5])/clock_freq);
    printf("   S: %.4f seconds\n", ((double)profile_counter[6])/clock_freq);
    printf(" sum: %.4f seconds\n\n", ((double)sum)/clock_freq);
  }
}

/*
 * resetFailureCounter()
 *
 * Resets the failure counter.
 */
void resetFailureCounter(void)
{
  int i;

  for(i=0;i<4;i++){
    failure_counter[i] = 0;
  }
}

/*
 * startClock()
 *
 * Start the profiling clock.
 */
void startClock(void)
{
#ifdef _WIN32
  start = getClock();
#else
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
#endif
}

/*
 * stopClock()
 *
 * Stop the profiling clock.
 */
void stopClock(int which)
{
#ifdef _WIN32
  profile_counter[which] += getClock() - start;
#else
  struct timespec end;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
  profile_counter[which] += end.tv_nsec - start.tv_nsec;
#endif
}

/*
 * updateFailureCounter()
 *
 * Update the failure counter.
 */
void updateFailureCounter(int which)
{
  failure_counter[which] += 1;
}

/*
 * updateIterations()
 *
 * Update the iteration counter.
 */
void updateIterations(int iters)
{
  total_iterations += iters;
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
