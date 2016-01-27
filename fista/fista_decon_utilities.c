/*
 * 01/16
 *
 * Utility functions for processing 3D images that have been deconvolved using FISTA.
 *
 * Note that these routines are designed to work on double precision images where
 * we don't have to worry about two adjacent pixels having identical values.
 *
 * Hazen
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall fista_decon_utilities.c
 *  gcc -shared -Wl,-soname,fista_decon_utilities.so.1 -o fista_decon_utilities.so.1.0.1 fista_decon_utilities.o -lc
 *  ln -s fista_decon_utilities.so.1.0.1 fista_decon_utilities.so
 *
 * Windows:
 *  gcc -c ia_utilities.c
 *  gcc -shared -o ia_utilities.dll ia_utilities.o
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>

/* Function Declarations */
int isMaxima(double *, int, int, int, int, int, int);
int label(double *, int *, double, int, int, int, int);
void labelImage(double *, int *, double, int, int, int, int, int, int, int);
void moments(double *, double *, int *, int, int, int, int);

/* Functions */

/*
 * isMaxima(image, ci, cj, ck, z_size, y_size, x_size)
 *
 * Returns 1 if the image at (ci,cj,ck) is greater than
 * all of the surrounding pixels, otherwise 0 is returned.
 *
 * image - The image.
 * ci - x index.
 * cj - y index.
 * ck - z index.
 * x_size - Image size in x.
 * y_size - Image size in y.
 * z_size - Image size in z.
 *
 * Returns 1/0 if the point is a maxima, or not.
 */
int isMaxima(double *image, int ci, int cj, int ck, int x_size, int y_size, int z_size)
{
  int i,j,k,ti,tj,tk;
  double cur;

  cur = image[ci * (y_size * z_size) + cj * z_size + ck];
  for(i=-1;i<2;i++){
    ti = ci + i;
    if ((ti >= 0) && (ti < x_size)){
      for(j=-1;j<2;j++){
	tj = cj + j;
	if ((tj >= 0) && (tj < y_size)){
	  for(k=-1;k<2;k++){
	    tk = ck + k;
	    if ((tk >= 0) && (tk < z_size)){
	      if ((ti != 0) || (tj != 0) || (tk != 0)){
		if (cur < image[ti * (y_size * z_size) + tj * z_size + tk]){
		  return 0;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  return 1;
}

/*
 * label(image, labels, threshold, margin, z_size, y_size, x_size)
 *
 * Finds the locations of all the local maxima in an image and
 * gives them and any surrounding pixels of lesser value, but
 * greater than threshold the same label.
 *
 * image - The image to label.
 * labels - An array with the same dimensions as image to record the labels.
 * threshold - Minimum intensity threshold.
 * margin - Number of pixels around the edge (in x, y) to ignore.
 * x_size - Image size in x (slowest axis).
 * y_size - Image size in y (second slowest axis).
 * z_size - Image size in z (fast axis).
 *
 * Returns the number of unique labels.
 */
int label(double *image, int *labels, double threshold, int margin, int x_size, int y_size, int z_size)
{
  int i,j,k,t;
  int cur_label;

  cur_label = 1;
  for(i=margin;i<(x_size-margin);i++){    
    for(j=margin;j<(y_size-margin);j++){
      for(k=0;k<z_size;k++){
	t = i*(y_size*z_size)+j*z_size+k;
	if(image[t] > threshold){
	  if(isMaxima(image,i,j,k,x_size,y_size,z_size)){
	    labelImage(image,labels,threshold,cur_label,i,j,k,x_size,y_size,z_size);
	    cur_label += 1;
	  }
	}
      }
    }
  }

  cur_label -= 1;
  return cur_label;
}

/*
 * labelImage(image, labels, threshold, cur_label, ci, cj, ck, z_size, y_size, x_size)
 *
 * This is called recursively to do the actual labeling.
 *
 * image - The image to label.
 * labels - An array with the same dimensions as image to record the labels.
 * threshold - Minimum intensity threshold.
 * cur_label - The current label.
 * ci - x index.
 * cj - y index.
 * ck - z index.
 * x_size - Image size in x.
 * y_size - Image size in y.
 * z_size - Image size in z.
 */
void labelImage(double *image, int *labels, double threshold, int cur_label, int ci, int cj, int ck, int x_size, int y_size, int z_size)
{
  int i,j,k,t,ti,tj,tk;
  double cur;

  // Label current position.
  t = ci*(y_size*z_size) + cj*z_size + ck;
  labels[t] = cur_label;
  cur = image[t];

  // Check surrounding positions.
  for(i=-1;i<2;i++){
    ti = ci + i;
    if ((ti >= 0) && (ti < x_size)){
      for(j=-1;j<2;j++){
	tj = cj + j;
	if ((tj >= 0) && (tj < y_size)){
	  for(k=-1;k<2;k++){
	    tk = ck + k;
	    if ((tk >= 0) && (tk < z_size)){
	      t = ti*(y_size*z_size) + tj*z_size + tk;
	      if ((image[t] > threshold) && (image[t] < cur) && (labels[t] == 0)){
		labelImage(image,labels,threshold,cur_label,ti,tj,tk,x_size,y_size,z_size);
	      }
	    }
	  }
	}
      }
    }
  }
}

/*
 * moments(image, peaks, labels, counts, z_size, y_size, x_size)
 *
 * Compute the zeroth and first moments of each labeled section
 * of the image, i.e. sum and center of mass. The results are
 * stored in the peak array as [total1, cz1, cy1, cz1, total2, cz2, ..]
 *
 * image - The image to label.
 * peaks - The array to hold the results in, of size [number of labels, 4]
 * labels - An array with the same dimensions as image to record the labels.
 * counts - The number of labels.
 * x_size - Image size in x.
 * y_size - Image size in y.
 * z_size - Image size in z.
 */
void moments(double *image, double *peaks, int *labels, int counts, int x_size, int y_size, int z_size)
{
  int i,j,k,l,t;

  for(i=0;i<x_size;i++){
    for(j=0;j<y_size;j++){
      for(k=0;k<z_size;k++){
	t = i*(y_size*z_size) + j*z_size + k;
	if(labels[t]>0){
	  l = labels[t]-1;
	  peaks[4*l] += image[t];
	  peaks[4*l+1] += image[t]*(double)i;
	  peaks[4*l+2] += image[t]*(double)j;
	  peaks[4*l+3] += image[t]*(double)k;
	}
      }
    }
  }

  for(i=0;i<counts;i++){
    if(peaks[4*i]>0){
      peaks[4*i+1] = peaks[4*i+1]/peaks[4*i];
      peaks[4*i+2] = peaks[4*i+2]/peaks[4*i];
      peaks[4*i+3] = peaks[4*i+3]/peaks[4*i];
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
