/*
 * Utility functions for image analysis.
 *
 * Hazen 10/17
 */

/* Include */
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "multi_fit.h"

/*
 * This is matched in Python. A description of the fields is
 * provided with the Python version.
 */
typedef struct flmData
{
  int margin;
  int n_peaks;
  int radius;
  int z_range;
  
  int xsize;
  int ysize;
  int zsize;

  double threshold;
  
  double *z_values;

  int32_t **taken;
  double **images;
} flmData;


/* Function Declarations */
int calcMaxPeaks(flmData *);
void findLocalMaxima(flmData *, double *, double *, double *);
int isLocalMaxima(flmData *, double, int, int, int, int, int, int, int, int);

/*
 * calcMaxPeaks()
 *
 * Return the maximum number of peaks that could be in an image stack. This 
 * is just the number of pixels above threshold.
 */
int calcMaxPeaks(flmData *flm_data)
{
  int np,xi,yi,zi;

  np = 0;
  for(zi=0;zi<flm_data->zsize;zi++){
    for(yi=flm_data->margin;yi<(flm_data->ysize - flm_data->margin);yi++){
      for(xi=flm_data->margin;xi<(flm_data->xsize - flm_data->margin);xi++){
	if(flm_data->images[zi][yi*flm_data->xsize+xi]>flm_data->threshold){
	  if(flm_data->taken[zi][yi*flm_data->xsize+xi]==0){
	    np++;
	  }
	}
      }
    }
  }
  return np;
}

/*
 * findLocalMaxima()
 *
 * Finds the locations of all the local maxima in a stack of images with
 * intensity greater than threshold. Adds them to the list if that location 
 * has not already been used.
 *
 * Note: This destructively modifies the images.
 */
void findLocalMaxima(flmData *flm_data, double *z, double *y, double *x)
{
  int np,xi,yi,zi;
  int ex,ey,ez,sx,sy,sz;

  np = 0;
  for(zi=0;zi<flm_data->zsize;zi++){

    /* Set z search range. */
    sz = zi - flm_data->z_range;
    if(sz<0){ sz = 0;}
    ez = zi + flm_data->z_range;
    if(ez>=flm_data->zsize){ ez = flm_data->zsize-1; }
    
    for(yi=flm_data->margin;yi<(flm_data->ysize - flm_data->margin);yi++){

      /* Set y search range. */
      sy = yi - flm_data->radius;
      if(sy<0){ sy = 0; }
      ey = yi + flm_data->radius;
      if(ey>=flm_data->ysize){ ey = flm_data->ysize-1; }

      for(xi=flm_data->margin;xi<(flm_data->xsize - flm_data->margin);xi++){
	if(flm_data->images[zi][yi*flm_data->xsize+xi]>flm_data->threshold){
	  if(flm_data->taken[zi][yi*flm_data->xsize+xi]==0){

	    /* Set x search range. */
	    sx = xi - flm_data->radius;
	    if(sx<0){ sx = 0; }
	    ex = xi + flm_data->radius;
	    if(ex>=flm_data->xsize){ ex = flm_data->xsize-1; }

	    if(isLocalMaxima(flm_data, flm_data->images[zi][yi*flm_data->xsize+xi], sz, ez, sy, yi, ey, sx, xi, ex)){
	      flm_data->taken[zi][yi*flm_data->xsize+xi]++;
	      z[np] = flm_data->z_values[zi];
	      y[np] = yi;
	      x[np] = xi;
	      np++;
	    }

	    if (np >= flm_data->n_peaks){
	      printf("Warning! Found maximum number of peaks!\n");
	      return;
	    }
	  }
	}
      }
    }
  }

  flm_data->n_peaks = np;
}

/*
 * isLocalMaxima()
 *
 * Does a local search to check if the current pixel is a maximum. The search area
 * is a cylinder with it's axis pointed along the z axis.
 */
int isLocalMaxima(flmData *flm_data, double cur, int sz, int ez, int sy, int cy, int ey, int sx, int cx, int ex)
{
  int dx,dy,rr,xi,yi,zi;

  rr = flm_data->radius * flm_data->radius;
  
  for(zi=sz;zi<=ez;zi++){
    for(yi=sy;yi<=ey;yi++){
      dy = (yi - cy)*(yi - cy);
      for(xi=sx;xi<=sx;xi++){
	dx = (xi - cx)*(xi - cx);
	if((dx+dy)<rr){
	  if(flm_data->images[zi][yi*flm_data->xsize+xi]>cur){
	    return 0;
	  }
	  else{
	    flm_data->images[zi][yi*flm_data->xsize+xi] = flm_data->threshold - 0.1;
	  }
	}
      }
    }
  }
  return 1;
}


/*
 * The MIT License
 *
 * Copyright (c) 2017 Zhuang Lab, Harvard University
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
