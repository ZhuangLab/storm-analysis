/*
 * Draw gaussians on an array.
 * 
 * Hazen 03/13
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Define */
#define NPARAMS 5
#define XC 0
#define YC 1
#define HEIGHT 2
#define XW 3
#define YW 4

/* Function definitions */
void drawGaussians(double *, double *, int, int, int, int);

/* Functions */
void drawGaussians(double *image, double *gaussian_params, int image_x, int image_y, int number_gaussians, int resolution)
{
  int i, j, k;
  int awidth, fx, fy, sx, sy;
  double x, y, xw, yw;
  double dx, dy, px, py, sgx, sgy;
  double intens, sum;
  double norm, step_size, start;
  
  norm = 1.0/((double)(resolution*resolution));
  step_size = 1.0/((double)resolution);
  start = -0.5 + 0.5/((double)resolution);

  for(i=0; i < number_gaussians; i++){
    intens = gaussian_params[i*NPARAMS+HEIGHT];
    if(intens > 0.0){
      px = gaussian_params[i*NPARAMS+XC];
      py = gaussian_params[i*NPARAMS+YC];
      xw = gaussian_params[i*NPARAMS+XW];
      yw = gaussian_params[i*NPARAMS+YW];

      if(xw > yw){
	awidth = (int)(5.0 * xw);
      } else {
	awidth = (int)(5.0 * yw);
      }

      if(awidth < 1){
	awidth = 1;
      }
      
      sgx = 1.0/(2.0 * xw * xw);
      sgy = 1.0/(2.0 * yw * yw);

      sx = (int)px - awidth;
      sy = (int)py - awidth;
      fx = (int)px + awidth + 1;
      fy = (int)py + awidth + 1;
      if(sx < 0)       { sx = 0; }
      if(fx > image_y) { fx = image_y; }
      if(sy < 0)       { sy = 0; }
      if(fy > image_x) { fy = image_x; }

      for(j=sx;j<fx;j++){
	for(k=sy;k<fy;k++){
	  sum = 0.0;
	  for(dx=start;dx<0.5;dx+=step_size){
	    x = (((double)j)-px+dx)*(((double)j)-px+dx)*sgx;
	    for(dy=start;dy<0.5;dy+=step_size){
	      y = (((double)k)-py+dy)*(((double)k)-py+dy)*sgy;
	      sum += norm * intens * exp(-1.0 * (x + y));
	    }
	  }
	  image[j*image_x+k] += sum;
	}
      }
    }
  }
}

