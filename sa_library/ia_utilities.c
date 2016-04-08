/*
 * Utility functions for image analysis.
 *
 * Hazen
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall ia_utilities.c
 *  gcc -shared -Wl,-soname,ia_utilities.so.1 -o ia_utilities.so.1.0.1 ia_utilities.o -lc
 *
 * Windows:
 *  gcc -c ia_utilities.c
 *  gcc -shared -o ia_utilities.dll ia_utilities.o
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "multi_fit.h"


/* Function Declarations */
int findLocalMaxima(double *, int *, double *, double, double, int, int, int, int);
int getBackgroundIndex();
int getErrorIndex();
int getHeightIndex();
int getNPeakPar();
int getNResultsPar();
int getStatusIndex();
int getXCenterIndex();
int getXWidthIndex();
int getYCenterIndex();
int getYWidthIndex();
int getZCenterindex();
void initializePeaks(double *, double *, double *, double, double, int, int);
int mergeNewPeaks(double *, double *, double *, double, double, int, int);
void peakToPeakDist(double *, double *, double *, double *, double *, int, int);
void peakToPeakIndex(double *, double *, double *, double *, int *, int, int);
int removeClosePeaks(double *, double *, double, double, int);
int removeNeighbors(double *, double *, double, int);
void smoothImage(double *, int);


/* Functions */

/*
 * findLocalMaxima(image, taken, peaks, threshold, radius, image_size_x, image_size_y, margin, peak_size)
 *
 * Finds the locations of all the local maxima in an image with
 * intensity greater than threshold. Adds them to the list if
 * that location has not already been used.
 *
 * image - the image to analyze, assumed square.
 * taken - spots in the image where peaks have already been
 * peaks - pre-allocated storage for peak data.
 * threshold - minumum peak intensity.
 * radius - circle in which the peak is maximal.
 * image_size_x - size of the image in x (fast axis).
 * image_size_y - size of the image in y (slow axis).
 * margin - number of pixels around the edge to ignore.
 * peak_size - size of the peaks array.
 *
 * Returns the number of peaks found.
 */
int findLocalMaxima(double *image, int *taken, double *peaks, double threshold, double radius, int image_size_x, int image_size_y, int margin, int peak_size)
{
  int cnts,i,j,k,l,m,max,r,t;
  double tmp;

  cnts = 0;
  r = (int)(radius+0.5);
  radius = radius*radius;
  for(i=margin;i<(image_size_y-margin);i++){
    for(j=margin;j<(image_size_x-margin);j++){
      m = i*image_size_x+j;
      tmp = image[m];
      if(tmp>threshold){
	max = 1;
	k = -r;
	while((k<=r)&&max){
	  t = m+k*image_size_x;
	  for(l=-r;l<=r;l++){
	    if((k*k+l*l)<radius){
	      if((k<=0)&&(l<=0)){
		if(tmp<image[t+l]){
		  max=0;
		}
	      }
	      else{
		if(tmp<=image[t+l]){
		  max=0;
		}
	      }
	    }
	  }
	  k++;
	}
	if(max){
	  if((taken[m]<2)&&(cnts<peak_size)){
	    peaks[cnts*NRESULTSPAR+XCENTER] = j;
	    peaks[cnts*NRESULTSPAR+YCENTER] = i;
	    peaks[cnts*NRESULTSPAR+ZCENTER] = 0.0;
	    taken[m] += 1;
	    cnts++;
	  }
	}
      }
    }
  }

  return cnts;
}


/*
 * getBackgroundIndex.
 *
 * Returns the index of the background.
 */
int getBackgroundIndex()
{
  return BACKGROUND;
}


/*
 * getErrorIndex.
 *
 * Returns the index of the peak fit error.
 */
int getErrorIndex()
{
  return IERROR;
}


/*
 * getHeightIndex.
 *
 * Returns the index of the peak height parameter.
 */
int getHeightIndex()
{
  return HEIGHT;
}


/*
 * getNPeakPar.
 *
 * Returns the number of parameters in the peak fit array.
 */
int getNPeakPar()
{
  return NPEAKPAR;
}


/*
 * getNResultsPar.
 *
 * Returns the number of parameters in the results array.
 */
int getNResultsPar()
{
  return NRESULTSPAR;
}


/*
 * getStatusIndex.
 *
 * Returns the index of the status flag.
 */
int getStatusIndex()
{
  return STATUS;
}


/*
 * getXCenterIndex.
 *
 * Returns the index of the center in x.
 */
int getXCenterIndex()
{
  return XCENTER;
}

/*
 * getXWidthIndex.
 *
 * Returns the index of the peak width in x.
 */
int getXWidthIndex()
{
  return XWIDTH;
}

/*
 * getYCenterIndex.
 *
 * Returns the index of the center in y.
 */
int getYCenterIndex()
{
  return YCENTER;
}

/*
 * getYWidthIndex.
 *
 * Returns the index of the center in y.
 */
int getYWidthIndex()
{
  return YWIDTH;
}


/*
 * getZCenterIndex.
 *
 * Returns the index of the center in z.
 */
int getZCenterIndex()
{
  return ZCENTER;
}

/*
 * initializePeaks(peaks, image, background, sigma, n_peaks, x_size)
 *
 * Initialize peak values with the best guess of their height,
 * background and sigma.
 *
 * peaks - the list of peaks to initialize.
 * image - the image.
 * background - the estimated background.
 * sigma - the starting value for sigma.
 * zvalue - the starting value for the center z position.
 * n_peaks - the number of peaks.
 * x_size - the size of the image and background arrays in x (fast axis).
 */
void initializePeaks(double *peaks, double *image, double *background, double sigma, double zvalue, int n_peaks, int x_size)
{
  int i,j,k;

  for (i=0;i<n_peaks;i++){
    j = i*NRESULTSPAR;    
    k = peaks[j+YCENTER]*x_size + peaks[j+XCENTER];
    
    peaks[j+HEIGHT] = image[k] - background[k];
    peaks[j+XWIDTH] = sigma;
    peaks[j+YWIDTH] = sigma;
    peaks[j+BACKGROUND] = background[k];
    peaks[j+ZCENTER] = zvalue;
    peaks[j+STATUS] = RUNNING;
    peaks[j+IERROR] = 0.0;
  }
}

 /*
 * mergeNewPeaks(in_peaks, new_peaks, out_peaks, radius, num_in_peaks, num_new_peaks)
 *
 * Merge the current peaks list (in_peaks) with a new peak list
 * (new_peaks), putting the results in out_peaks (pre-allocated). 
 * Peaks from new_peaks that are at least a distance radius from 
 * peaks in in_peaks are added to the end of the out_peaks list.
 * The total number of peaks is returned. Old peaks that are near
 * newly added peaks are marked as running.
 *
 * in_peaks - the list of "good" peaks.
 * new_peaks - potential new peaks to add to the list of "good" peaks.
 * out_peaks - pre-allocated storage for the merged peak list.
 * radius - cut-off radius, new peaks that are closer than this
 *          distance to current peaks are ignored.
 * neighborhood - cut-off radius for flagging neighbors as running.
 * num_in_peaks - length of the in_peaks array (number of peaks).
 * num_new_peaks - length of the new_peaks array (number of peaks).
 *
 * Returns the number of peaks added to the out_peaks list.
 */
int mergeNewPeaks(double *in_peaks, double *new_peaks, double *out_peaks, double radius, double neighborhood, int num_in_peaks, int num_new_peaks)
{
  int bad,i,j,k;
  double dx,dy,rad,x,y;

  // first, copy in_peaks to out_peaks.
  for(i=0;i<num_in_peaks;i++){
    for(j=0;j<NRESULTSPAR;j++){
      out_peaks[i*NRESULTSPAR+j] = in_peaks[i*NRESULTSPAR+j];
    }
  }

  // then check new peaks & add if they are ok.
  radius = radius*radius;
  neighborhood = neighborhood*neighborhood;
  k = num_in_peaks;
  for(i=0;i<num_new_peaks;i++){
    x = new_peaks[i*NRESULTSPAR+XCENTER];
    y = new_peaks[i*NRESULTSPAR+YCENTER];
    bad = 0;
    j = 0;
    while((j<num_in_peaks)&&(!bad)){
      dx = x - in_peaks[j*NRESULTSPAR+XCENTER];
      dy = y - in_peaks[j*NRESULTSPAR+YCENTER];
      rad = dx*dx+dy*dy;
      if(rad<radius){
	bad = 1;
      }
      // FIXME: This could mark as running peaks that are 
      //   close to a bad peak which won't get added anyway.
      else if(rad<neighborhood){
	out_peaks[j*NRESULTSPAR+ZCENTER] = 0.0;
	out_peaks[j*NRESULTSPAR+STATUS] = RUNNING;
      }
      j++;
    }
    if(!bad){
      for(j=0;j<NRESULTSPAR;j++){
	out_peaks[k*NRESULTSPAR+j] = new_peaks[i*NRESULTSPAR+j];
      }
      k++;
    }
  }
  
  k -= num_in_peaks;
  return k;
}


/*
 * peakToPeakDist(x1, y1, x2, y2, dist1, n1, n2)
 *
 * Calculates the distance from each peak in 1 to the nearest peak in 2.
 * This doesn't have much to do with fitting, but it is useful for
 * evaluating the results.
 *
 * x1 - x locations of 1 peaks.
 * y1 - y locations of 1 peaks.
 * x2 - x locations of 2 peaks.
 * y2 - y locations of 2 peaks.
 * dist1 - distance to nearest peak in 2.
 * n1 - number of 1 peaks.
 * n2 - number of 2 peaks.
 */
void peakToPeakDist(double *x1, double *y1, double *x2, double *y2, double *dist1, int n1, int n2)
{
  int i,j;
  double x,y,best_d,d;
  
  for(i=0;i<n1;i++){
    x = x1[i];
    y = y1[i];
    best_d = (x-x2[0])*(x-x2[0])+(y-y2[0])*(y-y2[0]);
    for(j=1;j<n2;j++){
      d = (x-x2[j])*(x-x2[j])+(y-y2[j])*(y-y2[j]);
      if(d<best_d){
	best_d = d;
      }
    }
    dist1[i] = sqrt(best_d);
  }
}


/*
 * peakToPeakIndex(x1, y1, x2, y2, dist1, n1, n2)
 *
 * Returns the index of the closest peak in (x2, y2) to
 * each peak in (x1, y1).
 *
 * x1 - x locations of 1 peaks.
 * y1 - y locations of 1 peaks.
 * x2 - x locations of 2 peaks.
 * y2 - y locations of 2 peaks.
 * index - index of the nearest peak in 2.
 * n1 - number of 1 peaks.
 * n2 - number of 2 peaks.
 */
void peakToPeakIndex(double *x1, double *y1, double *x2, double *y2, int *index, int n1, int n2)
{
  int i,j,best_j;
  double x,y,best_d,d;
  
  for(i=0;i<n1;i++){
    x = x1[i];
    y = y1[i];
    best_j = 0;
    best_d = (x-x2[0])*(x-x2[0])+(y-y2[0])*(y-y2[0]);
    for(j=1;j<n2;j++){
      d = (x-x2[j])*(x-x2[j])+(y-y2[j])*(y-y2[j]);
      if(d<best_d){
	best_d = d;
	best_j = j;
      }
    }
    index[i] = best_j;
  }
}


/*
 * removeClosePeaks(in_peaks, out_peaks, radius, num_in_peaks)
 *
 * Construct a new list of peaks, out_peaks, which only includes
 * members of in_peaks that are at least radius away from their
 * brighter neighbors.
 *
 * in_peaks - the list peaks to process.
 * out_peaks - pre-allocated storage for the refined peak list.
 * radius - cut-off radius, peaks that are less than radius
 *          from a brighter neighbor are removed.
 * neighborhood - cut-off radius for flagging neighbors as running.
 * 
 * num_in_peaks - length of the in_peaks array (number of peaks).
 *
 * Returns the number of peaks added to the out_peaks list.
 *
 * FIXME: Intensity scaled radial cutoff?
 */
int removeClosePeaks(double *in_peaks, double *out_peaks, double radius, double neighborhood, int num_in_peaks)
{
  int bad,i,j,k;
  double dx,dy,h,x,y;

  radius = radius*radius;
  neighborhood = neighborhood*neighborhood;

  // 1. Flag the peaks to be removed.
  for(i=0;i<num_in_peaks;i++){
    x = in_peaks[i*NRESULTSPAR+XCENTER];
    y = in_peaks[i*NRESULTSPAR+YCENTER];
    h = in_peaks[i*NRESULTSPAR+HEIGHT];
    bad = 0;
    j = 0;
    while((j<num_in_peaks)&&(!bad)){
      if(j!=i){
	dx = x - in_peaks[j*NRESULTSPAR+XCENTER];
	dy = y - in_peaks[j*NRESULTSPAR+YCENTER];
	if(((dx*dx+dy*dy)<radius)&&(in_peaks[j*NRESULTSPAR+HEIGHT]>h)){
	  bad = 1;
	}
      }
      j++;
    }
    if(bad){
      in_peaks[i*NRESULTSPAR+STATUS] = BADPEAK;
    }
  }

  // 2. Flag (non-bad) neighbors of bad peaks as running.
  for(i=0;i<num_in_peaks;i++){
    if(in_peaks[i*NRESULTSPAR+STATUS]==BADPEAK){
      x = in_peaks[i*NRESULTSPAR+XCENTER];
      y = in_peaks[i*NRESULTSPAR+YCENTER];
      for(j=0;j<num_in_peaks;j++){
	if(j!=i){
	  dx = x - in_peaks[j*NRESULTSPAR+XCENTER];
	  dy = y - in_peaks[j*NRESULTSPAR+YCENTER];
	  if(((dx*dx+dy*dy)<neighborhood)&&(in_peaks[j*NRESULTSPAR+STATUS]!=BADPEAK)){
	    in_peaks[j*NRESULTSPAR+STATUS] = RUNNING;
	  }
	}
      }
    }
  }

  // 3. Create a new list with the bad peaks removed.
  
  k = 0;
  for(i=0;i<num_in_peaks;i++){

    // if the peak is not bad then add it to the list of out peaks.
    if(in_peaks[i*NRESULTSPAR+STATUS]!=BADPEAK){
      for(j=0;j<NRESULTSPAR;j++){
	out_peaks[k*NRESULTSPAR+j] = in_peaks[i*NRESULTSPAR+j];
      }
      k++;
    }
  }
  
  return k;
}


/*
 * removeNeighbors(in_peaks, out_peaks, radius, num_in_peaks)
 *
 * Construct a new list of peaks, out_peaks, which only includes
 * members of in_peaks that are at least radius away from all
 * of their neighbors.
 *
 * in_peaks - the list peaks to process.
 * out_peaks - pre-allocated storage for the refined peak list.
 * radius - cut-off radius, peaks that are less than radius
 *          from a brighter neighbor are removed.
 * num_in_peaks - length of the in_peaks array (number of peaks).
 *
 * Returns the number of peaks added to the out_peaks list.
 *
 */
int removeNeighbors(double *in_peaks, double *out_peaks, double radius, int num_in_peaks)
{
  int i,j,k,bad;
  double x,y,dx,dy;

  radius = radius*radius;
  k = 0;
  for(i=0;i<num_in_peaks;i++){
    x = in_peaks[i*NRESULTSPAR+XCENTER];
    y = in_peaks[i*NRESULTSPAR+YCENTER];
    bad = 0;
    j = 0;
    while((j<num_in_peaks)&&(!bad)){
      if(j!=i){
	dx = x - in_peaks[j*NRESULTSPAR+XCENTER];
	dy = y - in_peaks[j*NRESULTSPAR+YCENTER];
	if((dx*dx+dy*dy)<radius){
	  bad = 1;
	}
      }
      j++;
    }
    if(!bad){
      for(j=0;j<NRESULTSPAR;j++){
	out_peaks[k*NRESULTSPAR+j] = in_peaks[i*NRESULTSPAR+j];
      }
      k++;
    }
  }
  
  return k;
}


/*
 * smoothImage(image, image_size)
 *
 * Smooths an image by convolving with a gaussian.
 *
 * image - the image to analyze, assumed square.
 * image_size - size of the image (assumed square).
 *
 */
void smoothImage(double *image, int image_size)
{
  int i,j,k,l;
  double sum;

  // sigma = 0.5
  double gauss[] = {0.01134, 0.08382, 0.01134,
		    0.08382, 0.61935, 0.08382,
		    0.01134, 0.08382, 0.01134};

  for(i=1;i<(image_size-1);i++){
    for(j=1;j<(image_size-1);j++){
      sum = 0.0;
      for(k=0;k<3;k++){
	for(l=0;l<3;l++){
	  sum += image[(i+k-1)*image_size+(j+l-1)]*gauss[k*3+l];
	}
      }
      image[i*image_size+j] = sum;
    }
  }
}


/*
 * The MIT License
 *
 * Copyright (c) 2016 Zhuang Lab, Harvard University
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
