/*
 * C library for the affine transform.
 *
 * Hazen 05/17
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall affine_transform.c
 *  gcc -shared -Wl,-soname,affine_transform.so.1 -o affine_transform.so.1.0.1 affine_transform.o -lc
 *  ln -s affine_transform.so.1.0.1 affine_transform.so
 *
 * Windows:
 *  gcc -c -O3 affine_transform.c
 *  gcc -shared -o affine_transform.dll affine_transform.o
 */

#include <stdlib.h>
#include <stdio.h>


struct atrans_struct {
  double xt[3];
  double yt[3];
};
typedef struct atrans_struct atrans;


void cleanup(atrans *);
atrans *initialize(double *, double *);
void transform(atrans *, double *, double *, int, int);


void cleanup(atrans *at)
{
  free(at);
}


atrans *initialize(double *xt_in, double *yt_in)
{
  int i;
  atrans *at;

  at = (atrans *)malloc(sizeof(atrans));

  for(i=0;i<3;i++){
    at->xt[i] = xt_in[i];
    at->yt[i] = yt_in[i];
  }
  
  return at;
}

void transform(atrans *at, double *im, double *im_trans, int sx, int sy)
{
  int i,j,xi,yi;
  double dx,dy,im1,im2,im3,im4,tr,xf,yf;
  double *xt, *yt;

  xt = at->xt;
  yt = at->yt;
  for(i=0;i<sx;i++){
    for(j=0;j<sy;j++){
      xf = xt[0] + xt[1]*i + xt[2]*j;
      yf = yt[0] + yt[1]*i + yt[2]*j;
      xi = (int)xf;
      yi = (int)yf;
      
      if ((xi>0) && (xi < (sx-1)) && (yi>0) && (yi < (sy-1))){
	dx = xf - xi;
	dy = yf - yi;
	
	im1 = im[xi*sy+yi];
	im2 = im[xi*sy+yi+1];
	im3 = im[(xi+1)*sy+yi];
	im4 = im[(xi+1)*sy+yi+1];

	tr = (1.0-dx)*(1.0-dy)*im1 + (1.0-dx)*dy*im2 + dx*(1.0-dy)*im3 + dx*dy*im4;

	im_trans[i*sy+j] = tr;
      }
      else{
	im_trans[i*sy+j] = 0.0;
      }
    }
  }
}
