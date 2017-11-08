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

#include "../dbscan/kdtree.h"
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
struct kdtree *createKDTree(double *, double *, int);
void findLocalMaxima(flmData *, double *, double *, double *, double *);
void freeKDTree(struct kdtree *);
int isLocalMaxima(flmData *, double, int, int, int, int, int, int, int, int);
int markDimmerPeaks(double *, double *, double *, int32_t *, double, double, int);
void runningIfHasNeighbors(double *, double *, double *, double *, int32_t *, double, int, int);


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
	  if(flm_data->taken[zi][yi*flm_data->xsize+xi]<1){
	    np++;
	  }
	}
      }
    }
  }
  return np;
}

/*
 * createKDTree()
 *
 * Create a KD tree from two arrays of positions.
 */
struct kdtree *createKDTree(double *x, double *y, int n)
{
  int i;
  double pos[2];
  struct kdtree *kd;

  kd = kd_create(2);
  for(i=0;i<n;i++){
    pos[0] = x[i];
    pos[1] = y[i];
    kd_insert(kd, pos, (void *)(intptr_t)i);
  }

  return kd;
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
void findLocalMaxima(flmData *flm_data, double *z, double *y, double *x, double *h)
{
  int np,xi,yi,zi;
  int ex,ey,ez,sx,sy,sz;
  double cur;

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
	  if(flm_data->taken[zi][yi*flm_data->xsize+xi]<1){

	    /* Set x search range. */
	    sx = xi - flm_data->radius;
	    if(sx<0){ sx = 0; }
	    ex = xi + flm_data->radius;
	    if(ex>=flm_data->xsize){ ex = flm_data->xsize-1; }

	    cur = flm_data->images[zi][yi*flm_data->xsize+xi];
	    if(isLocalMaxima(flm_data, cur, sz, ez, sy, yi, ey, sx, xi, ex)){
	      flm_data->taken[zi][yi*flm_data->xsize+xi]++;
	      z[np] = flm_data->z_values[zi];
	      y[np] = yi;
	      x[np] = xi;
	      h[np] = cur;
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
 * freeKDTree()
 *
 * Frees an existing kdtree.
 */
void freeKDTree(struct kdtree *kd)
{
  kd_clear(kd);
  kd_free(kd);
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
      for(xi=sx;xi<=ex;xi++){
	dx = (xi - cx)*(xi - cx);
	if((dx+dy)<=rr){
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
 * markDimmerPeaks()
 *
 * For each peak, check if it has a brighter neighbor within radius, and if it
 * does mark the peak for removal (by setting the status to ERROR) and the 
 * neighbors as running.
 */
int markDimmerPeaks(double *x, double *y, double *h, int32_t *status, double r_removal, double r_neighbors, int np)
{
  int i,j,k;
  int is_dimmer, removed;
  double pos[2];
  struct kdres *set_r, *set_n;
  struct kdtree *kd;

  removed = 0;
  kd = createKDTree(x, y, np);

  for(i=0;i<np;i++){

    /* Skip error peaks. */
    if(status[i] == ERROR){
      continue;
    }

    /* Check for neighbors within r_removal. */
    pos[0] = x[i];
    pos[1] = y[i];
    set_r = kd_nearest_range(kd, pos, r_removal);

    /* Every point will have at least itself as a neighbor. */
    if(kd_res_size(set_r) < 2){
      kd_res_free(set_r);
      continue;
    }

    /* Check for brighter neighbors. */
    is_dimmer = 0;
    for(j=0;j<kd_res_size(set_r);j++){
      k = (intptr_t)kd_res_item_data(set_r);
      if(h[k] > h[i]){
	is_dimmer = 1;
	break;
      }
      kd_res_next(set_r);
    }
    kd_res_free(set_r);

    if(is_dimmer){
      removed++;
      status[i] = ERROR;

      /* Check for neighbors within r_neighbors. */
      set_n = kd_nearest_range(kd, pos, r_neighbors);
      for(j=0;j<kd_res_size(set_n);j++){
	k = (intptr_t)kd_res_item_data(set_n);
	if (status[k] == CONVERGED){
	  status[k] = RUNNING;
	}
	kd_res_next(set_n);
      }
      kd_res_free(set_n);
    }    
  }

  freeKDTree(kd);

  return removed;
}

/*
 * runningIfHasNeighbors()
 *
 * Update status based on proximity of new peaks (n_x, n_y) to current peaks (c_x, c_y).
 *
 * This works the simplest way by making a KD tree from the new peaks then comparing
 * the old peaks against this tree. However this might not be the fastest way given
 * that there will likely be a lot more current peaks then new peaks.
 */
void runningIfHasNeighbors(double *c_x, double *c_y, double *n_x, double *n_y, int32_t *status, double radius, int nc, int nn)
{
  int i;
  double pos[2];
  struct kdres *set;
  struct kdtree *kd;

  kd = createKDTree(n_x, n_y, nn);
  
  for(i=0;i<nc;i++){

    /* Skip RUNNING and ERROR peaks. */
    if((status[i] == RUNNING) || (status[i] == ERROR)){
      continue;
    }
    
    /* Check for new neighbors within radius. */
    pos[0] = c_x[i];
    pos[1] = c_y[i];
    set = kd_nearest_range(kd, pos, radius);

    if (kd_res_size(set) > 0){
      status[i] = RUNNING;
    }
    
    kd_res_free(set);
  }
  
  freeKDTree(kd);
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
