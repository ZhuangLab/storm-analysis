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


atrans *initialize(double *xt, double *yt)
{
  int i;
  atrans *at;

  at = (atrans *)malloc(sizeof(atrans));

  for(i=0;i<3;i++){
    at->xt[i] = xt[i];
    at->yt[i] = yt[i];
  }
  
  return at;
}

void transform(atrans *at, double *im, double *im_trans, int sy, int sx)
{
  int i,j,xi,yi;
  double dx,dy,im1,im2,im3,im4,tr,xf,yf;
  double *xt, *yt;

  xt = at->xt;
  yt = at->yt;
  for(i=0;i<sy;i++){
    for(j=0;j<sx;j++){
      xf = xt[0] + xt[1]*j + xt[2]*i;
      yf = yt[0] + yt[1]*j + yt[2]*i;
      xi = (int)xf;
      yi = (int)yf;
      
      if ((xi>0) && (xi < (sx-1)) && (yi>0) && (yi < (sy-1))){
	dx = xf - xi;
	dy = yf - yi;
	
	im1 = im[yi*sx+xi];
	im2 = im[yi*sx+xi+1];
	im3 = im[(yi+1)*sx+xi];
	im4 = im[(yi+1)*sx+xi+1];

	tr = (1.0-dy)*(1.0-dx)*im1 + (1.0-dy)*dx*im2 + dy*(1.0-dx)*im3 + dy*dx*im4;

	im_trans[i*sx+j] = tr;
      }
      else{
	im_trans[i*sx+j] = 0.0;
      }
    }
  }
}
