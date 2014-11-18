/*
 * That which is common to all homotopy based image analyzers.
 *
 * Hazen 04/13
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "homotopy_common.h"
#include "homotopy_imagea_common.h"

/* Define */
#define MAXITERS 800
#define NUMBERNONZERO 50

/* Structures */

/* Function Declarations */
void removePeak(double *, int *, int, int, double, double, double *, double *, double *);

/*
 * closeFile()
 *
 * Close the file for saving the high resolution data.
 */
void closeFile(void)
{
  fclose(hres_fp);
}

/*
 * finishUp()
 *
 * Free space allocated for image analysis.
 */
void finishUp(void)
{
  /* Homotopy library cleanup. */
  cleanup();
 
  free(xvec);
  free(yvec);
}

/*
 * getPeaks()
 *
 * Return all the peaks in a hres image. Zeros out the
 * hres image in the process of doing this.
 *
 * hres - the high resolution image.
 * peak_x - pre-allocated storage for peak x positions.
 * peak_y - pre-allocated storage for peak y positions.
 * peak_i - pre-allocated storage for peak intensity.
 * peak_c - pre-allocated storage for the number of pixels in the peak.
 * max_peaks - size of the peak_x, peak_y, peak_i arrays.
 *
 * returns - number of peaks that were found.
 */
int getPeaks(double *hres, double *peak_x, double *peak_y, double *peak_i, int *peak_c, int max_peaks)
{
  int i,j,n,o_size,t1;
  double xt,yt,t,inv_scale;

  inv_scale = 1.0/((double)scale);
  n = 0;
  for(i=0;i<hres_y;i++){
    t1 = i*hres_x;
    for(j=0;j<hres_x;j++){
      if(hres[t1+j]>0.0){
	o_size = 0;
	xt = 0.0;
	yt = 0.0;
	t = 0.0;
	removePeak(hres, &o_size, i, j, 0.0, 0.0, &xt, &yt, &t);
	if((n<max_peaks)&&(t>0.0)){
	  peak_x[n] = inv_scale*((double)i + yt/t)+0.5+0.5*inv_scale;
	  peak_y[n] = inv_scale*((double)j + xt/t)+0.5+0.5*inv_scale;
	  peak_i[n] = t;
	  peak_c[n] = o_size;
	  n++;
	}
      }
    }
  }
  return n;
}

/*
 * openFile()
 *
 * Opens file for hres data, writes header & header data. This should
 * be called after the image parameters have been set.
 *
 * file_name - The name of the file to create for saving the data.
 *
 * Returns the last frame in the file, or zero if the file does not exist.
 */
int openFile(char *file_name)
{
  unsigned char data;
  int i,last_frame;

  last_frame = 0;
  if(access(file_name, F_OK) != -1){

    /* open the file for appending. */
    hres_fp = fopen(file_name, "rb+");

    /* determine frame number of last object added. */
    /*
     * FIXME: should probably use a 64 bit version of fseek.
     */
    fseek(hres_fp,-12,SEEK_END);
    fread(&last_frame, sizeof(int), 1, hres_fp);

    /* set fp to the end of the file. */
    fseek(hres_fp,0,SEEK_END);
  }
  else{
    /* open file */
    hres_fp = fopen(file_name, "wb");

    /*
     * write header, 100 bytes, mostly blank for now.. 
     */

    /* x size of the hres picture */
    fwrite(&hres_x,sizeof(int),1,hres_fp);

    /* y size of the hres picture */
    fwrite(&hres_y,sizeof(int),1,hres_fp);
    
    /* a bunch of blank space for possible future changes to the format */
    data = 0;
    for(i=0;i<92;i++){
      fwrite(&data,1,1,hres_fp);
    }
  }

  return last_frame;
}

/*
 * removePeak()
 *
 * Find the center of mass of a peak (all 8-connected pixels)
 * and remove the peak from the image.
 *
 * hres - the high resolution image.
 * o_size - pointer to storage for the number of elements in the object
 *          that were above the threshold.
 * i - current "x" location in the image.
 * j - current "y" location in the image.
 * dx - the current offset in x from point zero.
 * dy - the current offset in y from point zero.
 * xt - running sum of the intensity * dx.
 * yt - running sum of the intensity * dy.
 * t - running sum of the intensity.
 */
void removePeak(double *hres, int *o_size, int i, int j, double dx, double dy, double *xt, double *yt, double *t)
{
  double temp;

  temp = hres[i * hres_x + j];
  *t += temp;
  *xt += dx * temp;
  *yt += dy * temp;
  *o_size += 1;
  hres[i * hres_x + j] = 0.0;
  /* top row */
  if (i>0){
    if (j>0){
      if (hres[(i-1) * hres_x + (j-1)] > 0.0){
	removePeak(hres, o_size, i-1, j-1, dx-1.0, dy-1.0, xt, yt, t);
      }
    }
    if (hres[(i-1) * hres_x + j] > 0.0){
      removePeak(hres, o_size, i-1, j, dx, dy-1.0, xt, yt, t);
    }
    if (j<(hres_x-1)){
      if (hres[(i-1) * hres_x + (j+1)] > 0.0){
	removePeak(hres, o_size, i-1, j+1, dx+1.0, dy-1.0, xt, yt, t);
      }
    }
  }
  /* middle */
  if (j>0){
    if (hres[i * hres_x + (j-1)] > 0.0){
      removePeak(hres, o_size, i, j-1, dx-1.0, dy, xt, yt, t);
    }
  }
  if (j<(hres_x-1)){
    if (hres[i * hres_x + (j+1)] > 0.0){
      removePeak(hres, o_size, i, j+1, dx+1.0, dy, xt, yt, t);
    }
  }
  /* bottom row */
  if (i<(hres_y-1)){
    if (j>0){
      if (hres[(i+1) * hres_x + (j-1)] > 0.0){
	removePeak(hres, o_size, i+1, j-1, dx-1.0, dy+1.0, xt, yt, t);
      }
    }
    if (hres[(i+1) * hres_x + j] > 0.0){
      removePeak(hres, o_size, i+1, j, dx, dy+1.0, xt, yt, t);
    }
    if (j<(hres_x-1)){
      if (hres[(i+1) * hres_x + (j+1)] > 0.0){
	removePeak(hres, o_size, i+1, j+1, dx+1.0, dy+1.0, xt, yt, t);
      }
    }
  }

}

/*
 * saveHighRes()
 *
 * Save the high resolution image in high-resolution format.
 *
 * hres - the high resolution image
 * frame - the associated frame number
 *
 * Returns the number of non-zero elements.
 */
int saveHighRes(double *hres, int frame)
{
  int cnt,i;
  float tmp;

  cnt = 0;
  for(i=0;i<(hres_x*hres_y);i++){
    if(hres[i]>0.0){
      cnt++;
      tmp = (float)hres[i];
      fwrite(&frame,sizeof(int),1,hres_fp);
      fwrite(&i,sizeof(int),1,hres_fp);
      fwrite(&tmp,sizeof(float),1,hres_fp);
    }
  }

  return cnt;
}

/*
 * See the accompanying license.txt file.
 *
 * Copyright (c) 2013 Zhuang Lab, Harvard University
 */
