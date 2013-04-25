/*
 * 12/10
 *
 * Image deconvolution following Lingenfelter 2009.
 * In the simplest case this is just Richardson-Lucy
 * deconvolution. The algorithm attempts to find the
 * background at the same time it is analyzing the
 * foreground.
 * 
 * "Forward" is calculating the fit from the deconvolution
 * coefficients. "Backward" is calculating how to update
 * the deconvolution coefficients.
 *
 * Static variable wickedly non-thread safe implementation!
 *
 *
 * 01/11
 *
 * Added cull of lowest coefficient in a neighborhood. The
 * thinking is that this is a way to approach a L0 norm.
 *
 *
 * 03/11
 *
 * Switched to using structs instead of static global
 * variables to facilitate being able to run multiple
 * deconvolutions at the same time (with some limitations).
 *
 *
 * 05/11
 *
 * Added start/stop technology so that I can do updates
 * following the CrowdSTORM approach.
 *
 *
 * 06/11
 *
 * Large scale code overhaul to make it (in theory anyway)
 * easier to implement arbritrary PSF deconvolution.
 * Mostly this is for the purpose of being able to 3D decon.
 *
 * Hazen
 *
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall mlem_sparse.c
 *  gcc -shared -Wl,-soname,mlem_sparse.so.1 -o mlem_sparse.so.1.0.1 mlem_sparse.o -lc
 *
 * Windows:
 *  gcc -c mlem_sparse.c
 *  gcc -shared -o mlem_sparse.dll mlem_sparse.o
 *
 */

/* Define */
#define PSFSIZE 4.0

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* Structures */
typedef struct
{
  int half_x_size;
  int half_y_size;
  int x_size;
  int y_size;
  double xc;
  double yc;
  double zc;
  double width_x;
  double width_y;
  double *coeff;
} gauss;

typedef struct
{
  int p_num_i;       /* number of pixels covered by the psf. */
  int p_loc;         /* location of the psf (in the high-res image). */
  double p_xj;       /* magnitude of the psf. */
  double p_comp;     /* compression factor of the psf (for variable compression). */
  int *p_i;          /* index to pixels in the image. */
  double *p_coeff;   /* magnitude scale factor associated with each pixel. */
} psf;

typedef struct
{
  int nbackground_psf;  /* number of background psfs. */
  int nforeground_psf;  /* number of foreground psfs. */

  /* data & fit */
  double *y_data;       /* image data. */
  double *y_fit;        /* current best fit. */

  psf **background_psf; /* array of pointers to background psf structures. */
  psf **foreground_psf; /* array of pointers to foreground psf structures. */
} decon;


/*
 * Function Declarations.
 */

/* Internal */
void _backwardFixedCompression(psf *, double *, double *, double);
void _backwardVariableCompression(psf *, double *, double *); 
void _calculateGaussian(gauss *);
void _calculateGaussianSigma(gauss *, double, double, double);
void _calculateGaussianZ(gauss *, double, double, double);
void _cleanupDecon(decon *);
void _cleanupDeconBack(decon *);
void _cleanupDeconFore(decon *);
void _cleanupGauss(gauss *);
void _cleanupPsf(psf *, int);
void _forward(psf *, double *);
gauss* _getGauss(int, int, int);
void _initializeGauss2D(double);
void _normalizePsf(psf *);

/* API */
void backward(decon *);
void backwardCompressed(decon *, double);
void backwardCompressedFixedBg(decon *, double);
void backwardVarCompressionFixedBg(decon *);
void cleanup(decon *);
void cull(decon *, double);
void forward(decon *);
double fractionLow(decon *, double);
void getAPsf(decon *, psf *, int);
void getBackground(decon *, double *);
double getDiff(decon *);
void getFit(decon *, double *);
void getForeground(decon *, double *);
void getGauss(gauss *, int, int, int);
double *getPeaks(decon *, int *, double, int, int);
void newBackground(decon *, double *, int);
void newForeground(decon *, double *, double);
void setCompression(decon *, double *);
//void setForeground(decon *, double *);
decon *setup2D(double, int, int);


/* Global Variables */

/* general */
int decon_size;
int image_size; /* assumed square */
int margin;
int scale;
int zsteps;
gauss **gpsf;


/*** Internal Functions ***/

/*
 * _backwardFixedCompression()
 *
 * Calculate ej and xj for a psf structure.
 */
void _backwardFixedCompression(psf *a_psf, double *data, double *fit, double compression)
{
  int i;
  double x_ej,x_xj;
  int *x_i;
  double *x_coeff;
  
  x_xj = a_psf->p_xj;
  x_i = a_psf->p_i;
  x_coeff = a_psf->p_coeff;

  x_ej = 0.0;
  for(i=0;i<(a_psf->p_num_i);i++){
    x_ej += x_coeff[i]*data[x_i[i]]/fit[x_i[i]];
  }

  x_xj = x_xj*x_ej*compression;
  if(x_xj<0.0){
    x_xj = 0.0;
  }

  a_psf->p_xj = x_xj;
}


/*
 * _backwardVariableCompression()
 *
 * Calculate ej and xj for a psf structure.
 */
void _backwardVariableCompression(psf *a_psf, double *data, double *fit)
{
  int i;
  double x_ej, x_xj;
  int *x_i;
  double *x_coeff;
  
  x_xj = a_psf->p_xj;
  x_i = a_psf->p_i;
  x_coeff = a_psf->p_coeff;

  x_ej = 0.0;

  for(i=0;i<(a_psf->p_num_i);i++){
    x_ej += x_coeff[i] * data[x_i[i]]/fit[x_i[i]];
  }

  x_xj = x_xj*x_ej*a_psf->p_comp;
  if(x_xj<0.0){
    x_xj = 0.0;
  }

  a_psf->p_xj = x_xj;
}


/*
 * _calculateGaussian()
 *
 * This allocates space for & calculates a normalized gaussian.
 *
 * gaussian - Pointer gaussian structure.
 */
void _calculateGaussian(gauss *gaussian)
{
  int i, j;
  double total, sum, sgx, sgy, dx, dy, x, y, cx, cy;

  gaussian->coeff = (double *)malloc(sizeof(double)*gaussian->x_size*gaussian->y_size);

  total = 0.0;
  cx = gaussian->xc;
  cy = gaussian->yc;
  sgx = 2.0/(gaussian->width_x*gaussian->width_x);
  sgy = 2.0/(gaussian->width_y*gaussian->width_y);
  for(i=0;i<gaussian->y_size;i++){
    for(j=0;j<gaussian->x_size;j++){
      sum = 0.0;
      for(dy=-0.40;dy<0.5;dy+=0.2){
	y = (i-cy+dy)*(i-cy+dy)*sgy;
	for(dx=-0.40;dx<0.5;dx+=0.2){
	  x = (j-cx+dx)*(j-cx+dx)*sgx;
	  sum += 0.04*exp(-1.0*(x+y));
	}
      }
      gaussian->coeff[i*gaussian->x_size+j] = sum;
      total += sum;
    }
  }

  total = 1.0/total;
  for(i=0;i<(gaussian->x_size*gaussian->y_size);i++){
    gaussian->coeff[i] = gaussian->coeff[i]*total;
  }
}


/*
 * _calculateGaussianSigma()
 *
 * calculate gaussian from x-center, y-center and sigma.
 *
 * gaussian - Pointer gaussian structure.
 * cx - center in x.
 * cy - center in y.
 * sigma - sigma.
 */
void _calculateGaussianSigma(gauss *gaussian, double cx, double cy, double sigma)
{
  int offset,size;
  
  offset = (int)(PSFSIZE*sigma);
  size = 2 * offset + 1;
  gaussian->half_x_size = offset;
  gaussian->half_y_size = offset;
  gaussian->x_size = size;
  gaussian->y_size = size;
  gaussian->xc = cx + (double)offset;
  gaussian->yc = cy + (double)offset;
  gaussian->zc = 0.0;
  gaussian->width_x = 2.0*sigma;
  gaussian->width_y = 2.0*sigma;
  _calculateGaussian(gaussian);
}


/*
 * _calculateGaussianSigma()
 *
 * calculate gaussian from x-center, y-center and z-center.
 *
 * gaussian - Pointer gaussian structure.
 * cx - center in x.
 * cy - center in y.
 * cz - center in z.
 */
void _calculateGaussianZ(gauss *gaussian, double cx, double cy, double cz)
{
  gaussian->x_size = 1;
  gaussian->y_size = 1;
  gaussian->xc = cx;
  gaussian->yc = cy;
  gaussian->zc = cz;
  gaussian->width_x = 1.0;
  gaussian->width_y = 1.0;
  _calculateGaussian(gaussian);
}

/*
 * _cleanupDecon()
 *
 * Free space allocated for a decon struct.
 */
void _cleanupDecon(decon *dec)
{
  _cleanupDeconBack(dec);
  _cleanupDeconFore(dec);

  free(dec->y_fit);
  free(dec);
}



/*
 * _cleanupDeconBack()
 *
 * Free space allocated for the background psfs.
 */
void _cleanupDeconBack(decon *dec)
{
  int i;

  for(i=0;i<dec->nbackground_psf;i++){
    _cleanupPsf(dec->background_psf[i], 1);
  }
  free(dec->background_psf);
  dec->nbackground_psf = 0;
}


/*
 * _cleanupDeconFore()
 *
 * Free space allocated for the foreground psfs.
 */
void _cleanupDeconFore(decon *dec)
{
  int i;

  for(i=0;i<dec->nforeground_psf;i++){
    _cleanupPsf(dec->foreground_psf[i], 0);
  }
  free(dec->foreground_psf);
  dec->nforeground_psf = 0;
}


/*
 * _cleanupGauss()
 *
 * Free space allocated for a gauss struct.
 */
void _cleanupGauss(gauss *a_gauss)
{
  free(a_gauss->coeff);
  free(a_gauss);
}


/*
 * _cleanupPsf()
 *
 * Free space allocated for a psf struct.
 */
void _cleanupPsf(psf *a_psf, int free_coeff)
{
  free(a_psf->p_i);
  if(free_coeff){
    free(a_psf->p_coeff);
  }
  free(a_psf);
}


/*
 * _forward()
 *
 * Calculate yi(x) from psf.
 */
void _forward(psf *a_psf, double *fit)
{
  int i;
  double x_xj;
  int *x_i;
  double *x_coeff;

  x_xj = a_psf->p_xj;
  x_i = a_psf->p_i;
  x_coeff = a_psf->p_coeff;

  for(i=0;i<(a_psf->p_num_i);i++){
    fit[x_i[i]] += x_coeff[i] * x_xj;
  }
}


/*
 * _getGauss
 *
 * Returns the appropriate gauss structure.
 *
 * gx - gaussian x index
 * gy - gaussian y index
 * gz - gaussian z index
 *
 */
gauss* _getGauss(int gx, int gy, int gz)
{
  return gpsf[gz*scale*scale + gy*scale + gx];
}


/*
 * _initializeGauss2D()
 *
 * Setups up the arrays for the pre-calculated gaussian.
 *
 * sigma - gaussian sigma
 */
void _initializeGauss2D(double sigma)
{
  int i,j;
  double x,y,step_size,start_p;

  margin = (int)(PSFSIZE*sigma);
  zsteps = 1;
  gpsf = (gauss **)malloc(sizeof(gauss *)*scale*scale);
  step_size = 1.0/((double)scale);
  start_p = -0.5+0.5*step_size;

  for(i=0;i<scale;i++){
    for(j=0;j<scale;j++){
      y = ((double)i)*step_size + start_p;
      x = ((double)j)*step_size + start_p;
      gpsf[i*scale+j] = (gauss *)malloc(sizeof(gauss));
      _calculateGaussianSigma(gpsf[i*scale+j], x, y, sigma);
    }
  }
}


/*
 * _normalizePsf()
 *
 * Normalize a psf.
 */
void _normalizePsf(psf *a_psf)
{
  int i;
  double sum;

  for(i=0;i<a_psf->p_num_i;i++){
    sum += a_psf->p_coeff[i];
  }
  sum = 1.0/sum;
  for(i=0;i<a_psf->p_num_i;i++){
    a_psf->p_coeff[i] = a_psf->p_coeff[i]*sum;
  }
}



/*** API Functions ***/


/*
 * backward()
 *
 * Calculate ej and xj.
 */
void backward(decon *dec)
{
  int i;
  psf *a_psf;
  
  for(i=0;i<dec->nbackground_psf;i++){
    _backwardFixedCompression(dec->background_psf[i], dec->y_data, dec->y_fit, 1.0);    
  }

  for(i=0;i<dec->nforeground_psf;i++){
    a_psf = dec->foreground_psf[i];
    if(a_psf->p_xj > 0.0){
      _backwardFixedCompression(a_psf, dec->y_data, dec->y_fit, 1.0);    
    }
  }
}


/*
 * backwardCompressed()
 *
 * Calculate ej and xj w/ compression parameter 
 * following Lingenfelter L1 norm.
 *
 * cmpr - compression parameter
 */
void backwardCompressed(decon *dec, double compression)
{
  int i;
  double cmpr;
  psf *a_psf;

  cmpr = 1.0/(1.0 + compression);

  for(i=0;i<dec->nbackground_psf;i++){
    _backwardFixedCompression(dec->background_psf[i], dec->y_data, dec->y_fit, 1.0);    
  }

  for(i=0;i<dec->nforeground_psf;i++){
    a_psf = dec->foreground_psf[i];
    if(a_psf->p_xj > 0.0){
      _backwardFixedCompression(a_psf, dec->y_data, dec->y_fit, cmpr);
    }
  }
}


/*
 * backwardCompressedFixedBg()
 *
 * Calculate ej and xj w/ compression parameter 
 * following Lingenfelter L1 norm, but do not
 * update the background.
 *
 * cmpr - compression parameter
 */
void backwardCompressedFixedBg(decon *dec, double compression)
{
  int i;
  double cmpr;
  psf *a_psf;

  cmpr = 1.0/(1.0 + compression);

  for(i=0;i<dec->nforeground_psf;i++){
    a_psf = dec->foreground_psf[i];
    if(a_psf->p_xj > 0.0){
      _backwardFixedCompression(a_psf, dec->y_data, dec->y_fit, cmpr);
    }
  }
}


/*
 * backwardVarCompressionFixedBg()
 *
 * Calculate ej and xj w/ spatial variable
 * compression parameter.
 *
 */
void backwardVarCompressionFixedBg(decon *dec)
{
  int i;
  psf *a_psf;

  for(i=0;i<dec->nforeground_psf;i++){
    a_psf = dec->foreground_psf[i];
    if(a_psf->p_xj > 0.0){
      _backwardVariableCompression(a_psf, dec->y_data, dec->y_fit);
    }
  }
}


/*
 * cleanup()
 * 
 * Free space allocated for global variables and
 * the decon structure, should be called when
 * finished with the analysis.
 */
void cleanup(decon *dec)
{
  int i;

  for(i=0;i<(scale*scale*zsteps);i++){
    _cleanupGauss(gpsf[i]);
  }
  free(gpsf);

  _cleanupDecon(dec);
}


/*
 * cull(threshold)
 *
 * Zero out foreground psfs that are less than threshold.
 *
 * threshold - all psfs that are less than this value will be 
 *             zeroed out (and subsequently ignored).
 */
void cull(decon *dec, double threshold)
{
  int i;
  psf *a_psf;

  for(i=0;i<dec->nforeground_psf;i++){
    a_psf = dec->foreground_psf[i];
    if(a_psf->p_xj < threshold){
      a_psf->p_xj = 0.0;
    }
  }
}


/*
 * forward(threshold)
 * 
 * Calculate y_fit from fit coefficients.
 *
 */
void forward(decon *dec)
{
  int i;
  double *fit;
  psf *a_psf;

  fit = dec->y_fit;

  for(i=0;i<(image_size*image_size);i++){
    fit[i] = 0.0;
  }

  for(i=0;i<dec->nbackground_psf;i++){
    _forward(dec->background_psf[i], fit);
  }

  for(i=0;i<dec->nforeground_psf;i++){
    a_psf = dec->foreground_psf[i];
    if(a_psf->p_xj > 0.0){
      _forward(a_psf, fit);
    }
  }
}

/*
 * fractionLow(threshold)
 * 
 * Return the fraction of (foreground) xj that are above threshold.
 *
 * threshold - minimum xj value.
 */
double fractionLow(decon *dec, double threshold)
{
  int i, cnts;
  double tmp;

  /* count xj below threshold. */
  cnts = 0;
  for(i=0;i<(dec->nforeground_psf);i++){
    if(dec->foreground_psf[i]->p_xj < threshold){
      cnts++;
    }
  }

  tmp = 0.0;
  if(dec->nforeground_psf>0){
    tmp = (float)cnts/(float)dec->nforeground_psf;
  }
  return tmp;
}


/*
 * getAPsf
 *
 * Returns values of a psf in a user supplied structure.
 *
 * which - index of psf of interest.
 */
void getAPsf(decon *dec, psf *a_psf, int which)
{
  int i;
  psf *which_psf;

  which_psf = dec->foreground_psf[which];
  a_psf->p_num_i = which_psf->p_num_i;
  a_psf->p_loc = which_psf->p_loc;
  a_psf->p_xj = which_psf->p_xj;
  a_psf->p_comp = which_psf->p_comp;

  for(i=0;i<which_psf->p_num_i;i++){
    a_psf->p_i[i] = which_psf->p_i[i];
    a_psf->p_coeff[i] = which_psf->p_coeff[i];
  }
}


/*
 * getBackground
 *
 * Returns background in a user supplied array.
 *
 * back - pointer to storage for the background data (of type double).
 */
void getBackground(decon *dec, double *back)
{
  int i;

  for(i=0;i<(image_size*image_size);i++){
    back[i] = 0.0;
  }

  for(i=0;i<dec->nbackground_psf;i++){
    _forward(dec->background_psf[i], back);
  }
}


/*
 * getDiff
 *
 * Returns the (absolute) difference between the image and the fit.
 */
double getDiff(decon *dec)
{
  int i;
  double diff;
  
  diff = 0.0;
  for(i=0;i<(image_size*image_size);i++){
    diff += fabs(dec->y_data[i]-dec->y_fit[i]);
  }
  return diff;
}


/*
 * getFit
 *
 * Returns the fit image in a user supplied array.
 *
 * fit - pointer to storage for the fit data (of type double).
 */
void getFit(decon *dec, double *fit)
{
  int i;

  for(i=0;i<(image_size*image_size);i++){
    fit[i] = dec->y_fit[i];
  }
}


/*
 * getForeground
 *
 * Returns foreground in a user supplied array.
 *
 * fore - pointer to storage for the foreground data (of type double).
 */
void getForeground(decon *dec, double *fore)
{
  int i;
  psf *a_psf;

  for(i=0;i<(decon_size*decon_size);i++){
    fore[i] = 0.0;
  }

  for(i=0;i<(dec->nforeground_psf);i++){
    a_psf = dec->foreground_psf[i];
    fore[a_psf->p_loc] += a_psf->p_xj;
  }
}


/*
 * getGauss
 *
 * Returns values of a deconvolution gaussian 
 * in a user supplied structure.
 *
 * gx - gaussian x index
 * gy - gaussian y index
 * gz - gaussian z index
 *
 */
void getGauss(gauss *store, int gx, int gy, int gz)
{
  int i;
  gauss *t_gauss;

  t_gauss = _getGauss(gx,gy,gz);
  store->x_size = t_gauss->x_size;
  store->y_size = t_gauss->y_size;
  store->xc = t_gauss->xc;
  store->yc = t_gauss->yc;
  store->zc = t_gauss->zc;
  store->width_x = t_gauss->width_x;
  store->width_y = t_gauss->width_y;
  for(i=0;i<(t_gauss->x_size*t_gauss->y_size);i++){
    store->coeff[i] = t_gauss->coeff[i];
  }
}


/*
 * getPeaks
 *
 * Returns stats of all the peaks in the deconvolved image.
 *
 * n - pointer to integer, will hold the number of peaks found.
 * threshold - minimum peak height of interest.
 * guard_size - don't return peaks closer than this to the edge of the image.
 * n_size - size of neighborhood to average over for peak location.
 *
 * returns - [height1, xc1, yc1, bg1, height2, xc2, yc2, bg2, ...]
 */

#define MAXPEAKS 2000
double *getPeaks(decon *dec, int *n, double threshold, int guard_size, int n_size)
{
  int i,j,k,l,np;
  double tmp,xsum,ysum,sum,inv_scale,offset;
  int *xci,*yci;
  double *fore,*peaks;

  np = 0;
  xci = (int *)malloc(sizeof(int)*MAXPEAKS);
  yci = (int *)malloc(sizeof(int)*MAXPEAKS);
  fore = (double *)malloc(sizeof(double)*decon_size*decon_size);
  getForeground(dec, fore);

  // find local maxima above threshold
  for(i=(scale*guard_size);i<(decon_size-scale*guard_size);i++){
    l = i*decon_size;
    for(j=(scale*guard_size);j<(decon_size-scale*guard_size);j++){
      if(fore[l+j]>threshold){
	tmp = fore[l+j];
	if((tmp>fore[(l-decon_size)+j-1])&&
	   (tmp>fore[(l-decon_size)+j])&&
	   (tmp>fore[(l-decon_size)+j+1])&&
	   (tmp>fore[l+j-1])&&
	   (tmp>fore[l+j+1])&&
	   (tmp>fore[(l+decon_size)+j-1])&&
	   (tmp>fore[(l+decon_size)+j])&&
	   (tmp>fore[(l+decon_size)+j+1])){
	  yci[np] = i;
	  xci[np] = j;
	  np++;
	}
      }
    }
  }
  *n = np;

  // compute stats over the neighborhood
  inv_scale = 1.0/(double)scale;
  // offset = 0.0;
  offset = 0.5-0.5*inv_scale;
  peaks = (double *)malloc(sizeof(double)*4*np);
  for(i=0;i<np;i++){
    xsum = 0.0;
    ysum = 0.0;
    sum = 0.0;
    for(j=-n_size;j<=n_size;j++){
      l = (yci[i]+j)*decon_size;
      for(k=-n_size;k<=n_size;k++){
	tmp = fore[l+xci[i]+k];
	sum += tmp;
	xsum += tmp*(double)k;
	ysum += tmp*(double)j;
      }
    }
    tmp = 1.0/sum;
    //peaks[4*i] = sum/(sigma*sigma*2.0*3.14159);
    peaks[4*i] = sum/(2.0*3.14159);
    peaks[4*i+1] = (xsum*tmp + (double)xci[i])*inv_scale - offset;
    peaks[4*i+2] = (ysum*tmp + (double)yci[i])*inv_scale - offset;
    peaks[4*i+3] = 0.0;
    // printf("%f %f %f\n",peaks[4*i],peaks[4*i+1],peaks[4*i+2]);
  }

  free(xci);
  free(yci);
  free(fore);

  return peaks;
}


/*
 * newBackground()
 *
 * Setups up all the arrays for MLEM fitting of the
 * background.
 *
 * Background is fit with a piecewise linear function.
 * Basically a bunch of pyramids centered on grid points.
 *
 * gridsize must be chosen so that image_size % GRIDSIZE = 1
 */

void newBackground(decon *dec, double *background, int gridsize)
{
  int i,j,k,wi,wj,b_size,temp;
  double dx,dy;
  psf *a_psf;

  if(dec->nbackground_psf>0){
    _cleanupDeconBack(dec);
  }

  b_size = image_size/gridsize+1;
  dec->nbackground_psf = b_size*b_size;
  dec->background_psf = (psf **)malloc(sizeof(psf *)*dec->nbackground_psf);
  dec->y_data = background;

  for(i=0;i<b_size;i++){
    for(j=0;j<b_size;j++){
      a_psf = (psf *)malloc(sizeof(psf));
      a_psf->p_num_i = 0;      
      a_psf->p_i = (int *)malloc(sizeof(int)*4*gridsize*gridsize);
      a_psf->p_coeff = (double *)malloc(sizeof(double)*4*gridsize*gridsize);
      dec->background_psf[i*b_size+j] = a_psf;
    }
  }

  for(i=0;i<image_size;i++){
    temp = i*image_size;
    wi = i/gridsize;
    dy = (double)(i - wi * gridsize)/(double)gridsize;
    for(j=0;j<image_size;j++){
      wj = j/gridsize;
      dx = (double)(j - wj * gridsize)/(double)gridsize;

      // a
      a_psf = dec->background_psf[wi*b_size+wj];
      k = a_psf->p_num_i;
      a_psf->p_i[k] = temp+j;
      a_psf->p_coeff[k] = (1.0-dx)*(1.0-dy);
      a_psf->p_num_i++;

      // b
      if((wj<b_size)&&(dx>0.0)){
	a_psf = dec->background_psf[wi*b_size+(wj+1)];
	k = a_psf->p_num_i;
	a_psf->p_i[k] = temp+j;
	a_psf->p_coeff[k] = dx*(1.0-dy);
	a_psf->p_num_i++;
      }

      // c
      if((wi<b_size)&&(dy>0.0)){
	a_psf = dec->background_psf[(wi+1)*b_size+wj];
	k = a_psf->p_num_i;
	a_psf->p_i[k] = temp+j;
	a_psf->p_coeff[k] = (1.0-dx)*dy;
	a_psf->p_num_i++;
      }

      // d
      if ((wi<b_size)&&(wj<b_size)&&(dx>0.0)&&(dy>0.0)){
	a_psf = dec->background_psf[(wi+1)*b_size+(wj+1)];
	k = a_psf->p_num_i;
	a_psf->p_i[k] = temp+j;
	a_psf->p_coeff[k] = dx*dy;
	a_psf->p_num_i++;
      }
    }
  }

  /* normalize coefficients */
  for(i=0;i<(b_size*b_size);i++){
    _normalizePsf(dec->background_psf[i]);
  }

  /* calculate starting values */
  for(i=0;i<(b_size*b_size);i++){
    a_psf = dec->background_psf[i];
    a_psf->p_xj = 0.0;
    for(j=0;j<a_psf->p_num_i;j++){
      a_psf->p_xj += a_psf->p_coeff[j]*dec->y_data[a_psf->p_i[j]];
    }
  }
}


/*
 * newForeground()
 *
 * Setups up all the arrays for MLEM fitting of the
 * foreground.
 *
 * The foreground is fit with symmetric gaussian peaks.
 *
 * image - pointer to the image data (of type double).
 */
void newForeground(decon *dec, double *image, double threshold)
{
  int i,j,k,l,m,n,o,p,half_x,half_y,f_size,f_len;
  psf *a_psf;
  gauss *t_gauss;

  dec->y_data = image;

  /* 0. Free current storage. */
  if(dec->nforeground_psf>0){
    _cleanupDeconFore(dec);
  }

  /* 1. Determine number of points above threshold. */
  f_len = 0;
  for(i=margin;i<(image_size-margin);i++){
    for(j=margin;j<(image_size-margin);j++){
      if(image[i*image_size+j]>threshold){
	f_len++;
      }
    }
  }

  /* 2. Allocate space. */
  f_size = f_len*scale*scale*zsteps;
  dec->nforeground_psf = f_size;
  dec->foreground_psf = (psf **)malloc(sizeof(psf *)*f_size);

  /* 3. Initialize arrays. */
  p = 0;
  for(i=margin;i<(image_size-margin);i++){
    for(j=margin;j<(image_size-margin);j++){
      if(image[i*image_size+j]>threshold){
	for(k=0;k<scale;k++){
	  for(l=0;l<scale;l++){
	    for(m=0;m<zsteps;m++){
	      //a_psf = dec->foreground_psf[p];
	      a_psf = (psf *)malloc(sizeof(psf));
	      a_psf->p_loc = (i*scale+k)*decon_size+(j*scale+l);
	      a_psf->p_xj = image[i*image_size+j]/(double)scale;
	      a_psf->p_comp = 0.0;

	      t_gauss = _getGauss(l,k,m);
	      a_psf->p_num_i = t_gauss->x_size*t_gauss->y_size;
	      half_x = t_gauss->half_x_size;
	      half_y = t_gauss->half_y_size;
	      a_psf->p_i = (int *)malloc(sizeof(int)*t_gauss->x_size*t_gauss->y_size);
	      for(n=0;n<t_gauss->y_size;n++){
		for(o=0;o<t_gauss->x_size;o++){
		  a_psf->p_i[n*t_gauss->x_size+o] = (i+n-half_y)*image_size+(j+o-half_x);
		}
	      }
	      a_psf->p_coeff = t_gauss->coeff;
	      dec->foreground_psf[p] = a_psf;
	      p++;
	    }
	  }
	}
      }
    }
  }

  if(p!=f_size){
    printf("f_size init error: %d %d\n", p, f_size);
  }
}


/*
 * setCompression(comp)
 *
 * Set compression values to values given in the comp array,
 * which should be of size decon_size x decon_size.
 *
 * comp - compression value array.
 */
void setCompression(decon *dec, double *comp)
{
  int i;
  psf *a_psf;

  for(i=0;i<dec->nforeground_psf;i++){
    a_psf = dec->foreground_psf[i];
    a_psf->p_comp = comp[a_psf->p_loc];
  }
}


/*
 * setForeground(high_res)
 *
 * Set foreground with an existing high resolution image.
 *
 * high_res - high resolution image.
 */

/*
void setForeground(decon *dec, double *high_res)
{
  int i,j,k,l,ri,rj,f_len,f_size;
  double *t_gauss;
  int *x_i;
  int *x_j;
  int *x_loc;
  double *x_coeff;
  double *x_comp;
  double *x_xj;
  double *x_ej;

  / free old decon arrays /
  free(dec->f_i);
  free(dec->f_j);
  free(dec->f_loc);
  free(dec->f_coeff);
  free(dec->f_comp);
  free(dec->f_xj);
  free(dec->f_ej);

  / determine decon arrays size /
  f_size = 0;
  for(i=0;i<decon_size;i++){
    for(j=0;j<decon_size;j++){
      if(high_res[i*decon_size+j]>0.0){
	f_size++;
      }
    }
  }
  f_len = f_size*peak_arr_size*peak_arr_size;

  / allocate space /
  x_i = (int *)malloc(sizeof(int)*f_len);
  x_j = (int *)malloc(sizeof(int)*f_len);
  x_loc = (int *)malloc(sizeof(int)*f_size);
  x_coeff = (double *)malloc(sizeof(double)*f_len);
  x_comp = (double *)malloc(sizeof(double)*f_size);
  x_xj = (double *)malloc(sizeof(double)*f_size);
  x_ej = (double *)malloc(sizeof(double)*f_size);

  / generate decon arrays. /
  f_len = 0;
  f_size = 0;
  for(i=0;i<decon_size;i++){
    for(j=0;j<decon_size;j++){
      if(high_res[i*decon_size+j]>0.0){
	ri = i/scale;
	rj = j/scale;
	t_gauss = gauss[(i%scale)*scale+(j%scale)];
	for(k=0;k<peak_arr_size;k++){
	  for(l=0;l<peak_arr_size;l++){
	    x_i[f_len] = (ri+k-peak_size)*image_size+(rj+l-peak_size);
	    x_j[f_len] = f_size;
	    x_coeff[f_len] = t_gauss[k*peak_arr_size+l];
	    f_len++;	      
	  }
	}
	x_loc[f_size] = i*decon_size+j;
	x_comp[f_size] = 1.0;
	x_xj[f_size] = high_res[i*decon_size+j];
	f_size++;
      }
    }
  }

  / save results in decon struct. /
  dec->f_nij = f_len;
  dec->f_nj = f_size;
  dec->f_i = x_i;
  dec->f_j = x_j;
  dec->f_loc = x_loc;
  dec->f_coeff = x_coeff;
  dec->f_comp = x_comp;
  dec->f_xj = x_xj;
  dec->f_ej = x_ej;
}
*/

/*
 * setup2D
 *
 * Initializes things for 2D deconvolution.
 * This must be called before anything else.
 *
 * sigma_value - sigma of peaks in the image (assumed constant across the image).
 * image_size_value - size of the image (assumed square).
 * scale_value - the degree of "super-resolution" to attemp. (2 = 2x, 4 = 4x, etc..).
 *
 * returns a pointer to a initialized decon struct.
 */
decon *setup2D(double sigma_value, int image_size_value, int scale_value)
{
  decon *dec;

  decon_size = image_size_value*scale_value;
  image_size = image_size_value;
  scale = scale_value;
  _initializeGauss2D(sigma_value);

  dec = (decon *)malloc(sizeof(decon));
  dec->nbackground_psf = 0;
  dec->nforeground_psf = 0;
  dec->y_fit = (double *)malloc(sizeof(double)*image_size*image_size);

  return dec;
}
