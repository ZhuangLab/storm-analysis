/*
 * Calculate cubic spline values and derivatives given spline coefficients.
 * 
 * Notes:
 *  1. This uses static global variables to store the spline coefficients
 *     and other information about the current computation state. You cannot 
 *     have both a 2D and 3D spline active at the same time as they replace 
 *     each other at the init step.
 *
 *  2. Nominally xsize, ysize and zsize could be different, but this whether
 *     or not this works correctly has not been tested.
 *
 *  3. Unlike the Python version, this will not return the correct value if
 *     x is exactly at the maximum value that the spline covers.
 *
 *  4. This library is not thread safe.
 *
 * Hazen 12/13
 *
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall cubic_spline.c
 *  gcc -shared -Wl,-soname,cubic_spline.so.1 -o cubic_spline.so.1.0.1 cubic_spline.o -lc -lm
 *
 * Windows:
 *  gcc -c cubic_spline.c
 *  gcc -shared -o cubic_spline.dll cubic_spline.o
 */


/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cubic_spline.h"

/* Define */
#define RANGECHECK 1
#define RANGEWARN 0

/* Structures */

/* (Local) function declarations */

double dot(double *, double *, int);

double rangeCheckD(const char *, double, double, double);
int rangeCheckI(const char *, int, int, int);

/* Global variables */
static int xsize;
static int ysize;
static int zsize;
static double *aij;
static double *delta_f;
static double *delta_dxf;
static double *delta_dyf;
static double *delta_dzf;

/*
 * computeDelta2D()
 *
 * This computes the 16 cross-product values that you need to for 
 * 2D splines and their derivatives.
 * i.e. 1,x,y,xx,xy,yy,..
 *
 * y_delta - The delta in y (0.0 - 1.0).
 * x_delta - The delta in x (0.0 - 1.0).
 */
void computeDelta2D(double y_delta, double x_delta)
{
  int i,j;
  double cx,cy;

  if(RANGECHECK){
    x_delta = rangeCheckD("computeDelta2D,x_delta", x_delta, 0.0, 1.0);
    y_delta = rangeCheckD("computeDelta2D,y_delta", y_delta, 0.0, 1.0);
  }

  cx = 1.0;
  for(i=0;i<4;i++){
    cy = 1.0;
    for(j=0;j<4;j++){
      delta_f[4*i+j] = cx * cy;
      if(i<3){
	delta_dxf[4*(i+1)+j] = ((double)i+1) * cx * cy;
      }
      if(j<3){
	delta_dyf[4*i+j+1] = ((double)j+1) * cx * cy;
      }
      cy = cy * y_delta;
    }
    cx = cx * x_delta;
  }
}

/*
 * computeDelta3D()
 *
 * This computes the 64 cross-product values that you need to for 
 * 3D splines and their derivatives.
 * i.e. 1,x,y,z,xx,xy,yy,zx,zy,zz,..
 *
 * z_delta - The delta in z (0.0 - 1.0).
 * y_delta - The delta in y (0.0 - 1.0).
 * x_delta - The delta in x (0.0 - 1.0).
 */
void computeDelta3D(double z_delta, double y_delta, double x_delta)
{
  int i,j,k;
  double cx,cy,cz;

  if(RANGECHECK){
    x_delta = rangeCheckD("computeDelta3D,x_delta", x_delta, 0.0, 1.0);
    y_delta = rangeCheckD("computeDelta3D,y_delta", y_delta, 0.0, 1.0);
    z_delta = rangeCheckD("computeDelta3D,z_delta", z_delta, 0.0, 1.0);
  }

  cx = 1.0;
  for(i=0;i<4;i++){
    cy = 1.0;
    for(j=0;j<4;j++){
      cz = 1.0;
      for(k=0;k<4;k++){
	delta_f[i*16+j*4+k] = cx * cy * cz;
	if(i<3){
	  delta_dxf[(i+1)*16+j*4+k] = ((double)i+1) * cx * cy * cz;
	}
	if(j<3){
	  delta_dyf[i*16+(j+1)*4+k] = ((double)j+1) * cx * cy * cz;
	}
	if(k<3){
	  delta_dzf[i*16+j*4+k+1] = ((double)k+1) * cx * cy * cz;
	}
	cz = cz * z_delta;
      }
      cy = cy * y_delta;
    }
    cx = cx * x_delta;
  }
}

/*
 * dot()
 *
 * Compute the dot product of two vectors.
 *
 * v1 - The first vector.
 * v2 - The second vector.
 * len - The number of elements.
 *
 * Returns the dot product.
 */
double dot(double *v1, double *v2, int len)
{
  int i;
  double pd;

  pd = 0.0;
  for (i=0;i<len;i++){
    pd += v1[i]*v2[i];
  }

  return pd;
}

/*
 * dxfAt2D()
 *
 * Compute the derivative of the spline in x at x,y coordinate (integer). 
 * In order for this to work correctly computeDelta2D should have already 
 * been called.
 *
 * yc - The y coordinate (integer).
 * xc - The x coordinate (integer).
 *
 * Return - The derivative in x.
 */
double dxfAt2D(int yc, int xc)
{
  double yv;

  if(RANGECHECK){
    xc = rangeCheckI("dxfAt2D,xc", xc, 0, xsize);
    yc = rangeCheckI("dxfAt2D,yc", yc, 0, ysize);
  }

  yv = dot(&(aij[(xc*ysize+yc)*16]), delta_dxf, 16);

  return yv;
}

/*
 * dxfAt3D()
 *
 * Compute the derivative of the spline in x at x,y,z coordinate (integer). 
 * In order for this to work correctly computeDelta3D should have already 
 * been called.
 *
 * zc - The z coordinate (integer).
 * yc - The y coordinate (integer).
 * xc - The x coordinate (integer).
 *
 * Return - The derivative in y.
 */
double dxfAt3D(int zc, int yc, int xc)
{
  double yv;

  xc = rangeCheckI("dxfAt3D,xc", xc, 0, xsize);
  yc = rangeCheckI("dxfAt3D,yc", yc, 0, ysize);
  zc = rangeCheckI("dxfAt3D,zc", zc, 0, zsize);

  yv = dot(&(aij[(xc*(ysize*zsize)+yc*zsize+zc)*64]), delta_dxf, 64);

  return yv;
}

/*
 * dxfSpline2D()
 *
 * Return the derivative in x of the spline at coordinate x,y.
 *
 * y - The y coordinate (double).
 * x - The x coordinate (double).
 *
 * Return the derivative in x.
 */
double dxfSpline2D(double y, double x)
{
  int xc,yc;
  double x_delta,y_delta,yv;

  if(RANGECHECK){
    x = rangeCheckD("dxfSpline2D,x", x, 0.0, (double)(xsize));
    y = rangeCheckD("dxfSpline2D,y", y, 0.0, (double)(ysize));
  }

  xc = (int)x;
  x_delta = x - (double)xc;

  yc = (int)y;
  y_delta = y - (double)yc;

  computeDelta2D(y_delta, x_delta);
  yv = dxfAt2D(yc, xc);
  return yv;
}

/*
 * dxfSpline3D()
 *
 * Return the derivative in y of the spline at coordinate x,y,z.
 *
 * z - The x coordinate (double).
 * y - The y coordinate (double).
 * x - The x coordinate (double).
 */
double dxfSpline3D(double z, double y, double x)
{
  int xc,yc,zc;
  double x_delta,y_delta,z_delta,yv;

  if(RANGECHECK){
    x = rangeCheckD("dxfSpline3D,x", x, 0.0, (double)(xsize));
    y = rangeCheckD("dxfSpline3D,y", y, 0.0, (double)(ysize));
    z = rangeCheckD("dxfSpline3D,z", z, 0.0, (double)(zsize));
  }

  xc = (int)x;
  x_delta = x - (double)xc;

  yc = (int)y;
  y_delta = y - (double)yc;

  zc = (int)z;
  z_delta = z - (double)zc;

  computeDelta3D(z_delta, y_delta, x_delta);
  yv = dxfAt3D(zc, yc, xc);
  return yv;
}

/*
 * dyfAt2D()
 *
 * Compute the derivative of the spline in y at x,y coordinate (integer). 
 * In order for this to work correctly computeDelta2D should have already 
 * been called.
 *
 * yc - The y coordinate (integer).
 * xc - The x coordinate (integer).
 *
 * Return - The derivative in y.
 */
double dyfAt2D(int yc, int xc)
{
  double yv;

  if(RANGECHECK){
    xc = rangeCheckI("dyfAt2D,xc", xc, 0, xsize);
    yc = rangeCheckI("dyfAt2D,yc", yc, 0, ysize);
  }

  yv = dot(&(aij[(xc*ysize+yc)*16]), delta_dyf, 16);

  return yv;
}

/*
 * dyfAt3D()
 *
 * Compute the derivative of the in spline x at x,y,z coordinate (integer). 
 * In order for this to work correctly computeDelta3D should have already 
 * been called.
 *
 * zc - The z coordinate (integer).
 * yc - The y coordinate (integer).
 * xc - The x coordinate (integer).
 *
 * Return - The derivative in y.
 */
double dyfAt3D(int zc, int yc, int xc)
{
  double yv;

  if(RANGECHECK){
    xc = rangeCheckI("dyfAt3D,xc", xc, 0, xsize);
    yc = rangeCheckI("dyfAt3D,yc", yc, 0, ysize);
    zc = rangeCheckI("dyfAt3D,zc", zc, 0, zsize);
  }

  yv = dot(&(aij[(xc*(ysize*zsize)+yc*zsize+zc)*64]), delta_dyf, 64);

  return yv;
}

/*
 * dyfSpline2D()
 *
 * Return the derivative in y of the spline at coordinate x,y.
 *
 * y - The y coordinate (double).
 * x - The x coordinate (double).
 *
 * Return the derivative in y.
 */
double dyfSpline2D(double y, double x)
{
  int xc,yc;
  double x_delta,y_delta,yv;

  if(RANGECHECK){
    x = rangeCheckD("dyfSpline2D,x", x, 0.0, (double)(xsize));
    y = rangeCheckD("dyfSpline2D,y", y, 0.0, (double)(ysize));
  }

  xc = (int)x;
  x_delta = x - (double)xc;

  yc = (int)y;
  y_delta = y - (double)yc;

  computeDelta2D(y_delta, x_delta);
  yv = dyfAt2D(yc, xc);
  return yv;
}

/*
 * dyfSpline3D()
 *
 * Return the derivative in y of the spline at coordinate x,y,z.
 *
 * z - The x coordinate (double).
 * y - The y coordinate (double).
 * x - The x coordinate (double).
 */
double dyfSpline3D(double z, double y, double x)
{
  int xc,yc,zc;
  double x_delta,y_delta,z_delta,yv;

  if(RANGECHECK){
    x = rangeCheckD("dyfSpline3D,x", x, 0.0, (double)(xsize));
    y = rangeCheckD("dyfSpline3D,y", y, 0.0, (double)(ysize));
    z = rangeCheckD("dyfSpline3D,z", z, 0.0, (double)(zsize));
  }

  xc = (int)x;
  x_delta = x - (double)xc;

  yc = (int)y;
  y_delta = y - (double)yc;

  zc = (int)z;
  z_delta = z - (double)zc;

  computeDelta3D(z_delta, y_delta, x_delta);
  yv = dyfAt3D(zc, yc, xc);
  return yv;
}

/*
 * dzfAt3D()
 *
 * Compute the derivative of the in spline z at x,y,z coordinate (integer). 
 * In order for this to work correctly computeDelta3D should have already 
 * been called.
 *
 * zc - The z coordinate (integer).
 * yc - The y coordinate (integer).
 * xc - The x coordinate (integer).
 *
 * Return - The derivative in y.
 */
double dzfAt3D(int zc, int yc, int xc)
{
  double yv;

  if(RANGECHECK){
    xc = rangeCheckI("dzfAt3D,xc", xc, 0, xsize);
    yc = rangeCheckI("dzfAt3D,yc", yc, 0, ysize);
    zc = rangeCheckI("dzfAt3D,zc", zc, 0, zsize);
  }

  yv = dot(&(aij[(xc*(ysize*zsize)+yc*zsize+zc)*64]), delta_dzf, 64);

  return yv;
}

/*
 * dzfSpline3D()
 *
 * Return the derivative in z of the spline at coordinate x,y,z.
 *
 * z - The x coordinate (double).
 * y - The y coordinate (double).
 * x - The x coordinate (double).
 */
double dzfSpline3D(double z, double y, double x)
{
  int xc,yc,zc;
  double x_delta,y_delta,z_delta,yv;

  if(RANGECHECK){
    x = rangeCheckD("dzfSpline3D,x", x, 0.0, (double)(xsize));
    y = rangeCheckD("dzfSpline3D,y", y, 0.0, (double)(ysize));
    z = rangeCheckD("dzfSpline3D,z", z, 0.0, (double)(zsize));
  }

  xc = (int)x;
  x_delta = x - (double)xc;

  yc = (int)y;
  y_delta = y - (double)yc;

  zc = (int)z;
  z_delta = z - (double)zc;

  computeDelta3D(z_delta, y_delta, x_delta);
  yv = dzfAt3D(zc, yc, xc);
  return yv;
}

/*
 * fAt2D()
 *
 * Compute the spline at x,y coordinate (integer). In order
 * for this to work correctly computeDelta2D should have
 * already been called.
 *
 * yc - The y coordinate (integer).
 * xc - The x coordinate (integer).
 */
double fAt2D(int yc, int xc)
{
  double yv;

  if(RANGECHECK){
    xc = rangeCheckI("fAt2D,xc", xc, 0, xsize);
    yc = rangeCheckI("fAt2D,yc", yc, 0, ysize);
  }

  yv = dot(&(aij[(xc*ysize+yc)*16]), delta_f, 16);

  return yv;
}

/*
 * fAt3D()
 *
 * Compute the spline at x,y,z coordinate (integer). In order
 * for this to work correctly computeDelta3D should have
 * already been called.
 *
 * zc - The z coordinate (integer).
 * yc - The y coordinate (integer).
 * xc - The x coordinate (integer).
 */
double fAt3D(int zc, int yc, int xc)
{
  double yv;

  if(RANGECHECK){
    xc = rangeCheckI("fAt3D,xc", xc, 0, xsize);
    yc = rangeCheckI("fAt3D,yc", yc, 0, ysize);
    zc = rangeCheckI("fAt3D,zc", zc, 0, zsize);
  }

  yv = dot(&(aij[(xc*(ysize*zsize)+yc*zsize+zc)*64]), delta_f, 64);

  return yv;
}

/*
 * fSpline2D()
 *
 * Return the spline value at coordinate x,y.
 *
 * y - The y coordinate (double).
 * x - The x coordinate (double).
 */
double fSpline2D(double y, double x)
{
  int xc,yc;
  double x_delta,y_delta,yv;

  if(RANGECHECK){
    x = rangeCheckD("fSpline2D,x", x, 0.0, (double)(xsize));
    y = rangeCheckD("fSpline2D,y", y, 0.0, (double)(ysize));
  }

  xc = (int)x;
  x_delta = x - (double)xc;

  yc = (int)y;
  y_delta = y - (double)yc;

  computeDelta2D(y_delta, x_delta);
  yv = fAt2D(yc, xc);
  return yv;
}

/*
 * fSpline3D()
 *
 * Return the spline value at coordinate x,y,z.
 *
 * z - The z coordinate (double).
 * y - The y coordinate (double).
 * x - The x coordinate (double).
 */
double fSpline3D(double z, double y, double x)
{
  int xc,yc,zc;
  double x_delta,y_delta,z_delta,yv;

  if(RANGECHECK){
    x = rangeCheckD("fSpline3D,x", x, 0.0, (double)(xsize));
    y = rangeCheckD("fSpline3D,y", y, 0.0, (double)(ysize));
    z = rangeCheckD("fSpline3D,z", z, 0.0, (double)(zsize));
  }

  xc = (int)x;
  x_delta = x - (double)xc;

  yc = (int)y;
  y_delta = y - (double)yc;

  zc = (int)z;
  z_delta = z - (double)zc;

  computeDelta3D(z_delta, y_delta, x_delta);
  yv = fAt3D(zc, yc, xc);
  return yv;
}

/*
 * getXSize()
 *
 * Returns the x size of the current spline.
 */
int getXSize(void)
{
  return xsize;
}

/*
 * getYSize()
 *
 * Returns the y size of the current spline.
 */
int getYSize(void)
{
  return ysize;
}

/*
 * getZSize()
 *
 * Returns the z size of the current spline.
 */
int getZSize(void)
{
  return zsize;
}

/*
 * initSpline2D()
 *
 * Initialize the spline coefficient values.
 *
 * new_aij - The spline values.
 * new_xsize - The number of cells in x.
 * new_ysize - The number of cells in y.
 */
void initSpline2D(double *new_aij, int new_xsize, int new_ysize)
{
  int i;

  xsize = new_xsize;
  ysize = new_ysize;
  zsize = 0;

  /* Delete previous spline data. */
  splineCleanup();

  /* 
   * Allocate storage for new spline data.
   *
   * delta_dfz is not actually used, but we allocated it anyway to
   * keep the cleanup() function simple.
   */
  aij = (double *)malloc(sizeof(double)*xsize*ysize*16);
  delta_f = (double *)malloc(sizeof(double)*16);
  delta_dxf = (double *)malloc(sizeof(double)*16);
  delta_dyf = (double *)malloc(sizeof(double)*16);
  delta_dzf = (double *)malloc(sizeof(double));

  /* Copy spline coefficients. */
  for (i=0;i<(xsize*ysize*16);i++){
    aij[i] = new_aij[i];
  }

  /* Initialize delta arrays. */
  for (i=0;i<16;i++){
    delta_f[i] = 0.0;
    delta_dxf[i] = 0.0;
    delta_dyf[i] = 0.0;
  }
}

/*
 * initSpline3D()
 *
 * Initialize the spline coefficient values.
 *
 * new_aij - The spline values.
 * new_xsize - The number of cells in x.
 * new_ysize - The number of cells in y.
 * new_zsize - The number of cells in z.
 */
void initSpline3D(double *new_aij, int new_xsize, int new_ysize, int new_zsize)
{
  int i;

  xsize = new_xsize;
  ysize = new_ysize;
  zsize = new_zsize;

  /* Delete previous spline data. */
  splineCleanup();

  /* Allocate storage for new spline data. */
  aij = (double *)malloc(sizeof(double)*xsize*ysize*zsize*64);
  delta_f = (double *)malloc(sizeof(double)*64);
  delta_dxf = (double *)malloc(sizeof(double)*64);
  delta_dyf = (double *)malloc(sizeof(double)*64);
  delta_dzf = (double *)malloc(sizeof(double)*64);

  /* Copy spline coefficients. */
  for (i=0;i<(xsize*ysize*zsize*64);i++){
    aij[i] = new_aij[i];
  }

  /* Initialize delta arrays. */
  for (i=0;i<64;i++){
    delta_f[i] = 0.0;
    delta_dxf[i] = 0.0;
    delta_dyf[i] = 0.0;
    delta_dzf[i] = 0.0;
  }
}

/*
 * rangeCheckD()
 *
 * Checks whether the variable is in the expected range.
 * Prints a warning if it is not (double version).
 *
 * fname - The name of the current function.
 * value - The variable value.
 * vmin - The minimum expected variable value.
 * vmax - The maximum expected variable value.
 */
double rangeCheckD(const char *fname, double value, double vmin, double vmax)
{
  if (value < vmin){
    if (RANGEWARN){
      printf("Value out of range %f (%f) in %s\n", value, vmin, fname);
    }
    value = vmin;
  }
  if (value > vmax){
    if (RANGEWARN){
      printf("Value out of range %f (%f) in %s\n", value, vmax, fname);
    }
    value = vmax;
  }

  return value;
}

/*
 * rangeCheckI()
 *
 * Checks whether the variable is in the expected range.
 * Prints a warning if it is not (integer version).
 *
 * fname - The name of the current function.
 * value - The variable value.
 * vmin - The minimum expected variable value.
 * vmax - The maximum expected variable value.
 */
int rangeCheckI(const char *fname, int value, int vmin, int vmax)
{
  if (value < vmin){
    if (RANGEWARN){
      printf("Value out of range %d (%d) in %s\n", value, vmin, fname);
    }
    value = vmin;
  }
  if (value >= vmax){
    if (RANGEWARN){
      printf("Value out of range %d (%d) in %s\n", value, vmax, fname);
    }
    value = vmax - 1;
  }

  return value;
}

/*
 * splineCleanup()
 *
 * Free the allocated storage.
 */
void splineCleanup(void)
{
  if (aij != NULL){
    free(aij);
  }
  if (delta_f != NULL){
    free(delta_f);
    free(delta_dxf);
    free(delta_dyf);
    free(delta_dzf);
  }
}

