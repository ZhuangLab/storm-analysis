/*
 * Routine(s) for attempting to (MLE) fit multiple gaussians to an image.
 * The approach follows Laurence and Chromy, Nature Methods, 2010.
 *
 * Hazen 10/17
 */

/* Include */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "multi_fit.h"

#define MARGIN 10      /* Margin around the edge of the image. This
			  also sets the maximum number of pixels that will
			  be fit.
			  
			  FIXME: This should be adjustable. */

/* Structures */
typedef struct
{
  int wx;                       /* peak width in pixels in x (2 x wx + 1 = size_x). */
  int wy;                       /* peak width in pixels in y (2 x wy + 1 = size_y). */
  int xc;                       /* peak center in x in pixels (xc - wx = xi). */
  int yc;                       /* peak center in y in pixels (yc - wy = yi). */

  double wx_term;
  double wy_term;
  
  double xt[2*MARGIN+1];
  double ext[2*MARGIN+1];
  double yt[2*MARGIN+1];
  double eyt[2*MARGIN+1];
} daoPeak;

typedef struct
{
  int zfit;                     /* This is a flag for the 'Z' fitting model. */
  
  double min_z;                 /* minimum z value. */
  double max_z;                 /* maximum z value. */  

  double wx_z_params[5];        /* x width versus z parameters. */
  double wy_z_params[5];        /* y width versus z parameters. */
} daoFit;


/* Functions */
void daoAddPeak(fitData *);
peakData *daoAllocPeaks(int);
void daoCalcJH2DFixed(fitData *, double *, double *);
void daoCalcJH2D(fitData *, double *, double *);
void daoCalcJH3D(fitData *, double *, double *);
void daoCalcJHZ(fitData *, double *, double *);
void daoCalcLocSize(peakData *);
void daoCalcPeakShape(fitData *);
int daoCalcWidth(fitData *, double, int);
void daoCalcWidthsFromZ(fitData *, peakData *);
int daoCheck(fitData *);
int daoCheck2DFixed(fitData *);
int daoCheck2D(fitData *);
int daoCheck3D(fitData *);
int daoCheckZ(fitData *);
void daoCleanup(fitData *);
void daoCopyPeak(peakData *, peakData *);
void daoFreePeaks(peakData *, int);
fitData* daoInitialize(double *, double *, double, int, int);
void daoInitialize2DFixed(fitData *);
void daoInitialize2D(fitData *);
void daoInitialize3D(fitData *);
void daoInitializeZ(fitData *, double *, double *, double, double);
void daoNewPeaks(fitData *, double *, char *, int);
void daoSubtractPeak(fitData *);
void daoUpdate(peakData *);
void daoUpdate2DFixed(fitData *, double *);
void daoUpdate2D(fitData *, double *);
void daoUpdate3D(fitData *, double *);
void daoUpdateZ(fitData *, double *);


/*
 * daoAddPeak()
 *
 * Add working_peak peak to the foreground and background 
 * data arrays.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoAddPeak(fitData *fit_data)
{
  int j,k,l,m;
  double bg,mag,tmp;
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  peak->added++;

  /* Add peak to foreground and background arrays. */
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  mag = peak->params[HEIGHT];
  for(j=0;j<peak->size_y;j++){
    tmp = dao_peak->eyt[j];
    for(k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      fit_data->f_data[m] += mag * tmp * dao_peak->ext[k];
      fit_data->bg_counts[m] += 1;
      fit_data->bg_data[m] += bg + fit_data->scmos_term[m];
    }
  }
}


/*
 * daoAllocPeaks()
 *
 * Allocate storage for daoPeaks.
 */
struct peakData *daoAllocPeaks(int n_peaks)
{
  int i;
  peakData *new_peaks;

  new_peaks = (peakData *)malloc(sizeof(peakData)*n_peaks);  
  for(i=0;i<n_peaks;i++){
    new_peaks[i].peak_model = (daoPeak *)malloc(sizeof(daoPeak));
  }
  return new_peaks;
}


/* 
 * daoCalcJH2DFixed()
 *
 * Calculate Jacobian and Hessian for the 2DFixed model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DFixed(fitData *fit_data, double *jacobian, double *hessian)
{
  int j,k,l,m;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[4];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize Jacobian and Hessian. */
  for(j=0;j<4;j++){
    jacobian[j] = 0.0;
  }
  for(j=0;j<16;j++){
    hessian[j] = 0.0;
  }
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<peak->size_y;j++){
    yt = dao_peak->yt[j];
    eyt = dao_peak->eyt[j];
    for(k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
      xi = fit_data->x_data[m];
      xt = dao_peak->xt[k];
      ext = dao_peak->ext[k];
      e_t = ext*eyt;

      jt[0] = e_t;
      jt[1] = 2.0*a1*width*xt*e_t;
      jt[2] = 2.0*a1*width*yt*e_t;
      jt[3] = 1.0;
	
      // calculate jacobian
      t1 = 2.0*(1.0 - xi/fi);
      jacobian[0] += t1*jt[0];
      jacobian[1] += t1*jt[1];
      jacobian[2] += t1*jt[2];
      jacobian[3] += t1*jt[3];
      
      // calculate hessian without second derivative terms.
      t2 = 2.0*xi/(fi*fi);

      hessian[0] += t2*jt[0]*jt[0];
      hessian[1] += t2*jt[0]*jt[1];
      hessian[2] += t2*jt[0]*jt[2];
      hessian[3] += t2*jt[0]*jt[3];
	    
      // hessian[4]
      hessian[5] += t2*jt[1]*jt[1];
      hessian[6] += t2*jt[1]*jt[2];
      hessian[7] += t2*jt[1]*jt[3];
	
      // hessian[8]
      // hessian[9]
      hessian[10] += t2*jt[2]*jt[2];
      hessian[11] += t2*jt[2]*jt[3];
	
      // hessian[12]
      // hessian[13]
      // hessian[14]
      hessian[15] += t2*jt[3]*jt[3];
    }
  }
}  


/* 
 * daoCalcJH2D()
 *
 * Calculate Jacobian and Hessian for the 2D model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2D(fitData *fit_data, double *jacobian, double *hessian)
{
  int j,k,l,m;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;
  
  for(j=0;j<5;j++){
    jacobian[j] = 0.0;
  }
  for(j=0;j<25;j++){
    hessian[j] = 0.0;
  }
  l = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<peak->size_y;j++){
    yt = dao_peak->yt[j];
    eyt = dao_peak->eyt[j];
    for(k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
      xi = fit_data->x_data[m];
      xt = dao_peak->xt[k];
      ext = dao_peak->ext[k];
      e_t = ext*eyt;
      
      jt[0] = e_t;
      jt[1] = 2.0*a1*width*xt*e_t;
      jt[2] = 2.0*a1*width*yt*e_t;
      jt[3] = -a1*xt*xt*e_t-a1*yt*yt*e_t;
      jt[4] = 1.0;
	  
      // calculate jacobian
      t1 = 2.0*(1.0 - xi/fi);
      jacobian[0] += t1*jt[0];
      jacobian[1] += t1*jt[1];
      jacobian[2] += t1*jt[2];
      jacobian[3] += t1*jt[3];
      jacobian[4] += t1*jt[4];
      
      // calculate hessian without second derivative terms.
      t2 = 2.0*xi/(fi*fi);

      // calculate hessian without second derivative terms.
      hessian[0] += t2*jt[0]*jt[0];
      hessian[1] += t2*jt[0]*jt[1];
      hessian[2] += t2*jt[0]*jt[2];
      hessian[3] += t2*jt[0]*jt[3];
      hessian[4] += t2*jt[0]*jt[4];
	  
      // hessian[5]
      hessian[6] += t2*jt[1]*jt[1];
      hessian[7] += t2*jt[1]*jt[2];
      hessian[8] += t2*jt[1]*jt[3];
      hessian[9] += t2*jt[1]*jt[4];
	    
      // hessian[10]
      // hessian[11]
      hessian[12] += t2*jt[2]*jt[2];
      hessian[13] += t2*jt[2]*jt[3];
      hessian[14] += t2*jt[2]*jt[4];
	  
      // hessian[15]
      // hessian[16]
      // hessian[17]
      hessian[18] += t2*jt[3]*jt[3];
      hessian[19] += t2*jt[3]*jt[4];
	
      // hessian[20]
      // hessian[21]
      // hessian[22]
      // hessian[23]
      hessian[24] += t2*jt[4]*jt[4];
    }
  }
}


/* 
 * daoCalcJH3D()
 *
 * Calculate Jacobian and Hessian for the 3D model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH3D(fitData *fit_data, double *jacobian, double *hessian)
{
  int j,k,l,m;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,a3,a5;
  double jt[6];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  for(j=0;j<6;j++){
    jacobian[j] = 0.0;
  }
  for(j=0;j<36;j++){
    hessian[j] = 0.0;
  }
  l = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  a3 = peak->params[XWIDTH];
  a5 = peak->params[YWIDTH];
  for(j=0;j<peak->size_y;j++){
    yt = dao_peak->yt[j];
    eyt = dao_peak->eyt[j];
    for(k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
      xi = fit_data->x_data[m];
      xt = dao_peak->xt[k];
      ext = dao_peak->ext[k];
      e_t = ext*eyt;

      jt[0] = e_t;
      jt[1] = 2.0*a1*a3*xt*e_t;
      jt[2] = -a1*xt*xt*e_t;
      jt[3] = 2.0*a1*a5*yt*e_t;
      jt[4] = -a1*yt*yt*e_t;
      jt[5] = 1.0;
      
      // calculate jacobian
      t1 = 2.0*(1.0 - xi/fi);
      jacobian[0] += t1*jt[0];
      jacobian[1] += t1*jt[1];
      jacobian[2] += t1*jt[2];
      jacobian[3] += t1*jt[3];
      jacobian[4] += t1*jt[4];
      jacobian[5] += t1*jt[5];

      // calculate hessian without second derivative terms.
      t2 = 2.0*xi/(fi*fi);

      hessian[0] += t2*jt[0]*jt[0];
      hessian[1] += t2*jt[0]*jt[1];
      hessian[2] += t2*jt[0]*jt[2];
      hessian[3] += t2*jt[0]*jt[3];
      hessian[4] += t2*jt[0]*jt[4];
      hessian[5] += t2*jt[0]*jt[5];
	    
      // hessian[6]
      hessian[7]  += t2*jt[1]*jt[1];
      hessian[8]  += t2*jt[1]*jt[2];
      hessian[9]  += t2*jt[1]*jt[3];
      hessian[10] += t2*jt[1]*jt[4];
      hessian[11] += t2*jt[1]*jt[5];
	    
      // hessian[12]
      // hessian[13]
      hessian[14] += t2*jt[2]*jt[2];
      hessian[15] += t2*jt[2]*jt[3];
      hessian[16] += t2*jt[2]*jt[4];
      hessian[17] += t2*jt[2]*jt[5];
	
      // hessian[18]
      // hessian[19]
      // hessian[20]
      hessian[21] += t2*jt[3]*jt[3];
      hessian[22] += t2*jt[3]*jt[4];
      hessian[23] += t2*jt[3]*jt[5];
	    
      // hessian[24]
      // hessian[25]
      // hessian[26]
      // hessian[27]
      hessian[28] += t2*jt[4]*jt[4];
      hessian[29] += t2*jt[4]*jt[5];

      // hessian[30]
      // hessian[31]
      // hessian[32]
      // hessian[33]
      // hessian[34]
      hessian[35] += t2*jt[5]*jt[5];
    }
  }
}


/* 
 * daoCalcJHZ()
 *
 * Calculate Jacobian and Hessian for the Z model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJHZ(fitData *fit_data, double *jacobian, double *hessian)
{
  int j,k,l,m;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,a3,a5;
  double z0,z1,z2,zt,gx,gy;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;
  daoFit *dao_fit;
  
  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;
  dao_fit = (daoFit *)fit_data->fit_model;

  for(j=0;j<5;j++){
    jacobian[j] = 0.0;
  }
  for(j=0;j<25;j++){
    hessian[j] = 0.0;
  }
  l = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  a3 = peak->params[XWIDTH];
  a5 = peak->params[YWIDTH];
    
  // dwx/dz calcs
  z0 = (peak->params[ZCENTER] - dao_fit->wx_z_params[1]) / dao_fit->wx_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  zt = 2.0*z0 + 3.0*dao_fit->wx_z_params[3]*z1 + 4.0*dao_fit->wx_z_params[4]*z2;
  gx = -2.0*zt/(dao_fit->wx_z_params[0]*dao_peak->wx_term);

  // dwy/dz calcs
  z0 = (peak->params[ZCENTER] - dao_fit->wy_z_params[1]) / dao_fit->wy_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  zt = 2.0*z0 + 3.0*dao_fit->wy_z_params[3]*z1 + 4.0*dao_fit->wy_z_params[4]*z2;
  gy = -2.0*zt/(dao_fit->wy_z_params[0]*dao_peak->wy_term);
  for(j=0;j<peak->size_y;j++){
    yt = dao_peak->yt[j];
    eyt = dao_peak->eyt[j];
    for(k=0;k<peak->size_x;k++){
      m = j*fit_data->image_size_x + k + l;
      fi = fit_data->f_data[m] + fit_data->bg_data[m] / ((double)fit_data->bg_counts[m]);
      xi = fit_data->x_data[m];
      xt = dao_peak->xt[k];
      ext = dao_peak->ext[k];
      e_t = ext*eyt;
	
      // first derivatives
      jt[0] = e_t;
      jt[1] = 2.0*a1*a3*xt*e_t;
      jt[2] = 2.0*a1*a5*yt*e_t;
      jt[3] = -a1*xt*xt*gx*e_t-a1*yt*yt*gy*e_t;
      jt[4] = 1.0;
	  
      // calculate jacobian
      t1 = 2.0*(1.0 - xi/fi);
      jacobian[0] += t1*jt[0];
      jacobian[1] += t1*jt[1];
      jacobian[2] += t1*jt[2];
      jacobian[3] += t1*jt[3];
      jacobian[4] += t1*jt[4];
	
      // calculate hessian without second derivative terms.
      t2 = 2.0*xi/(fi*fi);
	
      hessian[0] += t2*jt[0]*jt[0];
      hessian[1] += t2*jt[0]*jt[1];
      hessian[2] += t2*jt[0]*jt[2];
      hessian[3] += t2*jt[0]*jt[3];
      hessian[4] += t2*jt[0]*jt[4];
	
      // hessian[5]
      hessian[6] += t2*jt[1]*jt[1];
      hessian[7] += t2*jt[1]*jt[2];
      hessian[8] += t2*jt[1]*jt[3];
      hessian[9] += t2*jt[1]*jt[4];
      
      // hessian[10]
      // hessian[11]
      hessian[12] += t2*jt[2]*jt[2];
      hessian[13] += t2*jt[2]*jt[3];
      hessian[14] += t2*jt[2]*jt[4];
	  
      // hessian[15]
      // hessian[16]
      // hessian[17]
      hessian[18] += t2*jt[3]*jt[3];
      hessian[19] += t2*jt[3]*jt[4];

      // hessian[20]
      // hessian[21]
      // hessian[22]
      // hessian[23]
      hessian[24] += t2*jt[4]*jt[4];
    }
  }
}


/*
 * daoCalcLocSize()
 *
 * Given the current center and width of the peak, calculate
 * the location of the start of the fitting area and the size
 * of the fitting area.
 *
 * For DAOSTORM it is easier to use the peak center and peak
 * width because the peak fitting area can change. However to 
 * reduce code redundancy we need the x, y location of the 
 * fitting area and its size for the shared functions.
 *
 * peak - pointer to a peakData structure.
 */
void daoCalcLocSize(peakData *peak)
{
  daoPeak *dao_peak;

  dao_peak = (daoPeak *)peak->peak_model;
  peak->xi = dao_peak->xc - dao_peak->wx;
  peak->yi = dao_peak->yc - dao_peak->wy;
  peak->size_x = 2 * dao_peak->wx + 1;
  peak->size_y = 2 * dao_peak->wy + 1;  
}


/*
 * daoCalcPeakShape()
 *
 * Calculate peak shape.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoCalcPeakShape(fitData *fit_data)
{
  int j,n,wx,wy,xc,yc;
  double xt,yt;
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  peak->added++;
    
  /* Calculate peak shape. */
  wx = dao_peak->wx;
  wy = dao_peak->wy;

  xc = dao_peak->xc;
  yc = dao_peak->yc;

  for(j=(xc-wx);j<=(xc+wx);j++){
    xt = (double)j - peak->params[XCENTER];
    n = j-xc+wx;
    dao_peak->xt[n] = xt;
    dao_peak->ext[n] = exp(-xt*xt*peak->params[XWIDTH]);
  }
  for(j=(yc-wy);j<=(yc+wy);j++){
    yt = (double)j - peak->params[YCENTER];
    n = j-yc+wy;
    dao_peak->yt[n] = yt;
    dao_peak->eyt[n] = exp(-yt*yt*peak->params[YWIDTH]);
  }
}

  
/*
 * daoCalcWidth()
 *
 * Given a peak_width, returns the appropriate 
 * bounding box to use for fitting.
 *
 * fit_data - pointer to a fitData structure.
 * peak_width - the new actual width of the peak.
 * old_w - the old integer width of the peak.
 */
int daoCalcWidth(fitData *fit_data, double peak_width, int old_w)
{
  int new_w;
  double tmp;

  if(peak_width < 0.0){
    if(TESTING){
      printf(" Got negative peak width! %.3f\n", peak_width);
    }
    return 1;
  }
  else{
    new_w = old_w;
    tmp = 4.0*sqrt(1.0/(2.0*peak_width));
    if(fabs(tmp - (double)old_w) > HYSTERESIS){
      new_w = (int)round(tmp);
    }
    if(new_w > fit_data->margin){
      new_w = fit_data->margin;
    }
    return new_w;
  }
}


/*
 * daoCalcWidthsFromZ()
 *
 * Updates wx, wy given z.
 *
 * fit_data - pointer to a fitData structure.
 * peak - Peak data structure to update.
 */
void daoCalcWidthsFromZ(fitData *fit_data, peakData *peak)
{
  double z0,z1,z2,z3,tmp;
  daoFit *dao_fit;
  daoPeak *dao_peak;

  dao_fit = (daoFit *)fit_data->fit_model;
  dao_peak = (daoPeak *)peak->peak_model;

  // wx
  z0 = (peak->params[ZCENTER] - dao_fit->wx_z_params[1]) / dao_fit->wx_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  z3 = z2*z0;
  tmp = 1.0 + z1 + dao_fit->wx_z_params[3] * z2 + dao_fit->wx_z_params[4] * z3;
  dao_peak->wx_term = tmp*tmp;
  peak->params[XWIDTH] = 2.0/(dao_fit->wx_z_params[0] * tmp);

  // wy
  z0 = (peak->params[ZCENTER] - dao_fit->wy_z_params[1]) / dao_fit->wy_z_params[2];
  z1 = z0*z0;
  z2 = z1*z0;
  z3 = z2*z0;
  tmp = 1.0 + z1 + dao_fit->wy_z_params[3] * z2 + dao_fit->wy_z_params[4] * z3;
  dao_peak->wy_term = tmp*tmp;
  peak->params[YWIDTH] = 2.0/(dao_fit->wy_z_params[0] * tmp);
}


/*
 * daoCheck()
 *
 * Check that the parameters of working_peak are still valid.
 *
 * fit_data - pointer to a fitData structure.
 *
 * Return 0 if okay.
 */
int daoCheck(fitData *fit_data)
{
  int margin,xc,yc;  
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /*
   * Check that the peak hasn't moved to close to the 
   * edge of the image. Flag the peak as bad if it has.
   */
  margin = fit_data->margin;
  xc = dao_peak->xc;
  yc = dao_peak->yc;
  if((xc <= margin)||(xc >= (fit_data->image_size_x - margin - 1))||(yc <= margin)||(yc >= (fit_data->image_size_y - margin - 1))){
    fit_data->n_margin++;
    if(TESTING){
      printf("object outside margins, %d %.5f, %.5f\n", peak->index, peak->params[XCENTER], peak->params[YCENTER]);
    }
    return 1;
  }
  
  /*
   * Check for negative height.
   */
  if(peak->params[HEIGHT]<0.0){
    fit_data->n_neg_height++;
    if(TESTING){
      printf("negative height, %d %.5f, %.5f (%.5f, %.5f)\n", peak->index, peak->params[BACKGROUND], peak->params[HEIGHT], peak->params[XCENTER], peak->params[YCENTER]);
    }
    return 1;
  }

  return 0;
}


/*
 * daoCheck2DFixed()
 *
 * Check that the parameters of working_peak are still valid (2DFixed model).
 *
 * fit_data - pointer to a fitData structure.
 *
 * Return 0 if okay.
 */
int daoCheck2DFixed(fitData *fit_data){
  return daoCheck(fit_data);
}


/*
 * daoCheck2D()
 *
 * Check that the parameters of working_peak are still valid (2D model).
 *
 * fit_data - pointer to a fitData structure.
 *
 * Return 0 if okay.
 */
int daoCheck2D(fitData *fit_data){
  peakData *peak;

  peak = fit_data->working_peak;
  
  if(daoCheck(fit_data)){
    return 1;
  }

  /*
   * Check for negative widths
   */
  if((peak->params[XWIDTH]<0.0)||(peak->params[YWIDTH]<0.0)){
    fit_data->n_neg_width++;
    if(TESTING){
      printf("negative widths, %.3f, %.3f (%.3f, %.3f)\n", peak->params[XWIDTH], peak->params[YWIDTH], peak->params[XCENTER], peak->params[YCENTER]);
    }
    return 1;
  }
  
  return 0;
}


/*
 * daoCheck3D()
 *
 * Check that the parameters of working_peak are still valid (3D model).
 *
 * fit_data - pointer to a fitData structure.
 *
 * Return 0 if okay.
 */
int daoCheck3D(fitData *fit_data){
  return daoCheck2D(fit_data);
}


/*
 * daoCheckZ()
 *
 * Check that the parameters of working_peak are still valid (Z model).
 *
 * fit_data - pointer to a fitData structure.
 *
 * Return 0 if okay.
 */
int daoCheckZ(fitData *fit_data){
  return daoCheck(fit_data);
}


/*
 * daoCleanup()
 *
 * Frees the fitData structure.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoCleanup(fitData *fit_data)
{
  int i;

  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->nfit;i++){
      free((daoPeak *)fit_data->fit[i].peak_model);
    }
    free(fit_data->fit);
  }
  free((daoPeak *)fit_data->working_peak->peak_model);
  free((daoFit *)fit_data->fit_model);
  mFitCleanup(fit_data);
}


/*
 * daoCopyPeak()
 *
 * Copies the contents of peak structure into another peak structure.
 *
 * original - pointer to a peakData structure.
 * copy - pointer to a peakData structure.
 */
void daoCopyPeak(peakData *original, peakData *copy)
{
  int i;
  daoPeak *dao_copy, *dao_original;

  dao_copy = (daoPeak *)copy->peak_model;
  dao_original = (daoPeak *)original->peak_model;

  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(original, copy);

  /* Copy the parts that are specific to 3D-DAOSTORM / sCMOS. */
  dao_copy->wx = dao_original->wx;
  dao_copy->wy = dao_original->wy;
  dao_copy->xc = dao_original->xc;
  dao_copy->yc = dao_original->yc;

  dao_copy->wx_term = dao_original->wx_term;
  dao_copy->wy_term = dao_original->wy_term;

  for(i=0;i<(2*MARGIN+1);i++){
    dao_copy->xt[i] = dao_original->xt[i];
    dao_copy->ext[i] = dao_original->ext[i];
    dao_copy->yt[i] = dao_original->yt[i];
    dao_copy->eyt[i] = dao_original->eyt[i];
  }
}


/*
 * daoFreePeaks()
 *
 * Frees a peakData array.
 *
 * peaks - Pointer to an array of peakData.
 * n_peaks - The size of the array.
 */
void daoFreePeaks(peakData *peaks, int n_peaks)
{
  int i;

  for(i=0;i<n_peaks;i++){
    free((daoPeak *)peaks[i].peak_model);
  }
  free(peaks);
}


/*
 * daoInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * clamp - The starting clamp values for each peak.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* daoInitialize(double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y)
{
  fitData* fit_data;

  fit_data = mFitInitialize(scmos_calibration, clamp, tol, im_size_x, im_size_y);

  fit_data->margin = MARGIN;
  fit_data->fit_model = (daoFit *)malloc(sizeof(daoFit));

  /* 
   * The default is not to do z fitting.
   *
   * We use this flag because it seems easier than having multiple 
   * versions of the daoNewPeaks() function.
   */
  ((daoFit *)fit_data->fit_model)->zfit = 0;
  
  /* Allocate storage for the working peak. */
  fit_data->working_peak->peak_model = (daoPeak *)malloc(sizeof(daoPeak));

  /* Set function pointers. */
  fit_data->fn_add_peak = &daoAddPeak;
  fit_data->fn_alloc_peaks = &daoAllocPeaks;
  fit_data->fn_calc_peak_shape = &daoCalcPeakShape;
  fit_data->fn_copy_peak = &daoCopyPeak;
  fit_data->fn_free_peaks = &daoFreePeaks;
  fit_data->fn_subtract_peak = &daoSubtractPeak;
  
  return fit_data;
}


/*
 * daoInitialize2DFixed()
 *
 * Initializes fitting for the 2DFixed model.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoInitialize2DFixed(fitData* fit_data)
{
  fit_data->jac_size = 4;
  
  fit_data->fn_calc_JH = &daoCalcJH2DFixed;
  fit_data->fn_check = &daoCheck2DFixed;
  fit_data->fn_update = &daoUpdate2DFixed;
}


/*
 * daoInitialize2D()
 *
 * Initializes fitting for the 2D model.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoInitialize2D(fitData* fit_data)
{
  fit_data->jac_size = 5;
    
  fit_data->fn_calc_JH = &daoCalcJH2D;
  fit_data->fn_check = &daoCheck2D;
  fit_data->fn_update = &daoUpdate2D;
}


/*
 * daoInitialize3D()
 *
 * Initializes fitting for the 3D model.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoInitialize3D(fitData* fit_data)
{
  fit_data->jac_size = 6;
    
  fit_data->fn_calc_JH = &daoCalcJH3D;
  fit_data->fn_check = &daoCheck3D;
  fit_data->fn_update = &daoUpdate3D;
}


/*
 * daoInitializeZ()
 *
 * Initializes fitting for the Z model with wx, wy dependence on z.
 *
 * fit_data - pointer to a fitData structure.
 * wx_vs_z - [wo, c, d, A, B] wx vs z coefficients.
 * wy_vs_z - [wo, c, d, A, B] wy vs z coefficients.
 * z_min - minimum allowed z value.
 * z_max - maximum allowed z value.
 */
void daoInitializeZ(fitData* fit_data, double *wx_vs_z, double *wy_vs_z, double z_min, double z_max)
{
  int i;
  daoFit *dao_fit;

  dao_fit = (daoFit *)fit_data->fit_model;

  fit_data->jac_size = 5;
    
  fit_data->fn_calc_JH = &daoCalcJHZ;
  fit_data->fn_check = &daoCheckZ;
  fit_data->fn_update = &daoUpdateZ;

  dao_fit->zfit = 1;
  for(i=0;i<5;i++){
    dao_fit->wx_z_params[i] = wx_vs_z[i];
    dao_fit->wy_z_params[i] = wy_vs_z[i];
  }
  dao_fit->wx_z_params[0] = dao_fit->wx_z_params[0] * dao_fit->wx_z_params[0];
  dao_fit->wy_z_params[0] = dao_fit->wy_z_params[0] * dao_fit->wy_z_params[0];

  dao_fit->min_z = z_min;
  dao_fit->max_z = z_max;
}


/*
 * daoNewPeaks
 *
 * fit_data - Pointer to a fitData structure.
 * peak_params - Input values for the peak parameters.
 * p_type - The type of the peak parameters.
 * n_peaks - The number of peaks.
 */
void daoNewPeaks(fitData *fit_data, double *peak_params, char *p_type, int n_peaks)
{
  int i,j,k,l,m,n;
  int start,stop;
  double width,sp,sx,t1;
  peakData *peak;
  daoPeak *dao_peak;

  if(VERBOSE){
    printf("dNP %d\n", n_peaks);
  }

  /* Generic initializations. */
  mFitNewPeaks(fit_data, n_peaks);

  /* 3D-DAOSTORM specific initializations. */
  start = fit_data->nfit;
  stop = fit_data->nfit + n_peaks;
  
  /*
   * 'finder' or 'testing' parameters, these are the peak x,y,z 
   * and sigma values as an n_peaks x 4 array.
   */
  if(!strcmp(p_type, "finder") || !strcmp(p_type, "testing")){
    for(i=start;i<stop;i++){
      j = 4*(i-start);
      peak = &fit_data->fit[i];
      dao_peak = (daoPeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;

      /* Initial width. */
      if(((daoFit *)fit_data->fit_model)->zfit){
	daoCalcWidthsFromZ(fit_data, peak);
      }
      else{
	width = 1.0/(2.0*peak_params[j+3]*peak_params[j+3]);
	peak->params[XWIDTH] = width;
	peak->params[YWIDTH] = width;
      }

      dao_peak->xc = (int)round(peak->params[XCENTER]);
      dao_peak->yc = (int)round(peak->params[YCENTER]);
      dao_peak->wx = daoCalcWidth(fit_data, peak->params[XWIDTH],-10.0);
      dao_peak->wy = daoCalcWidth(fit_data, peak->params[YWIDTH],-10.0);

      /* Calculate initial peak ROI. */
      daoCalcLocSize(peak);

      /* Estimate background. */
      peak->params[BACKGROUND] = fit_data->bg_estimate[dao_peak->yc * fit_data->image_size_x + dao_peak->xc];

      /* Arbitrary initial value for HEIGHT. */
      peak->params[HEIGHT] = 1.0;
      
      /* Copy into working peak. */
      daoCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	printf("Warning peak %d is bad!\n", (i-start));
	fit_data->working_peak->status = ERROR;
	daoCopyPeak(fit_data->working_peak, peak);
	continue;
      }

      if(!strcmp(p_type, "finder")){

	/* Calculate peak shape (of working peak). */
	daoCalcPeakShape(fit_data);

	/* 
	 * Estimate height. 
	 *
	 * Calculate the area under the peak of unit height and compare this to
	 * the area under (image - current fit - estimated background) x peak.
	 *
	 * We are minimizing : fi * (h*fi - xi)^2
	 * where fi is the peak shape at pixel i, h is the height and xi is
	 * is the data (minus the current fit & estimated background.
	 *
	 * Taking the derivative with respect to h gives fi*fi*xi/fi*fi*fi as the
	 * value for h that will minimize the above.
	 */
	dao_peak = (daoPeak *)fit_data->working_peak->peak_model;
	k = peak->yi * fit_data->image_size_x + peak->xi;
	sp = 0.0;  /* This is fi*fi*fi. */
	sx = 0.0;  /* This is fi*fi*xi. */
	for(l=0;l<peak->size_y;l++){
	  for(m=0;m<peak->size_x;m++){
	    n = l * fit_data->image_size_x + m + k;
	    t1 = dao_peak->eyt[l]*dao_peak->ext[m];
	    sp += t1*t1*t1;
	    sx += t1*t1*(fit_data->x_data[n] - fit_data->f_data[n] - fit_data->bg_estimate[n]);
	  }
	}
	fit_data->working_peak->params[HEIGHT] = sx/sp;
	
	/* Check that the initial height is positive, error it out if not. */
	if(fit_data->working_peak->params[HEIGHT] <= 0.0){
	  printf("Warning peak %d has negative estimated height!\n", (i-start));
	  fit_data->working_peak->status = ERROR;
	  daoCopyPeak(fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      daoAddPeak(fit_data);

      /* Copy values back from working peak. */
      daoCopyPeak(fit_data->working_peak, peak);
    }
  }
  /*
   * "pre-specified" parameters, these are the peak x,y,z and sigma 
   * values as an n_peaks x 7 array.
   */
  else if(!strcmp(p_type, "text") || !strcmp(p_type, "insight3")){
    for(i=start;i<stop;i++){
      j = 7*(i-start);
      peak = &fit_data->fit[i];
      dao_peak = (daoPeak *)peak->peak_model;

      /* Initial location. */
      peak->params[XCENTER] = peak_params[j] - fit_data->xoff;
      peak->params[YCENTER] = peak_params[j+1] - fit_data->yoff;
      peak->params[ZCENTER] = peak_params[j+2] - fit_data->zoff;
      peak->params[BACKGROUND] = peak_params[j+3];
      peak->params[HEIGHT] = peak_params[j+4];
      
      width = 1.0/(2.0*peak_params[j+5]*peak_params[j+5]);
      peak->params[XWIDTH] = width;

      width = 1.0/(2.0*peak_params[j+6]*peak_params[j+6]);
      peak->params[YWIDTH] = width;

      /* Correct initial width (z fitting model). */
      if(((daoFit *)fit_data->fit_model)->zfit){
	daoCalcWidthsFromZ(fit_data, peak);
      }

      dao_peak->xc = (int)round(peak->params[XCENTER]);
      dao_peak->yc = (int)round(peak->params[YCENTER]);
      dao_peak->wx = daoCalcWidth(fit_data, peak->params[XWIDTH],-10.0);
      dao_peak->wy = daoCalcWidth(fit_data, peak->params[YWIDTH],-10.0);

      /* Calculate initial peak ROI. */
      daoCalcLocSize(peak);
      
      /* Copy into working peak. */
      daoCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	printf("Warning peak %d is bad!\n", (i-start));
	fit_data->working_peak->status = ERROR;
	daoCopyPeak(fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      daoCalcPeakShape(fit_data);
      daoAddPeak(fit_data);

      /* Copy values back from working peak. */
      daoCopyPeak(fit_data->working_peak, peak);
    }
  }
  else{
    printf("Unrecognized peak type '%s'!\n", p_type);
  }

  /* Initial error calculation. */
  for(i=start;i<stop;i++){
    peak = &fit_data->fit[i];
    daoCopyPeak(peak, fit_data->working_peak);
    mFitCalcErr(fit_data);
    daoCopyPeak(fit_data->working_peak, peak);
  }

  fit_data->nfit = stop;

  /* Reset the clamp values on all the peaks. */
  if(USECLAMP){
    mFitResetClampValues(fit_data);
  }
}


/*
 * daoSubtractPeak()
 *
 * Subtract the peak out of the current fit, basically 
 * this just undoes addPeak().
 *
 * fit_data - pointer to a fitData structure.
 * peak - pointer to a peakData structure.
 */
void daoSubtractPeak(fitData *fit_data)
{
  int j,k,l,m;
  double bg,mag,tmp;
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  peak->added--;
  
  l = peak->yi * fit_data->image_size_x + peak->xi;
  bg = peak->params[BACKGROUND];
  mag = peak->params[HEIGHT];
  for(j=0;j<peak->size_y;j++){
    tmp = dao_peak->eyt[j];
    for(k=0;k<peak->size_x;k++){
      m = j * fit_data->image_size_x + k + l;
      fit_data->f_data[m] -= mag * tmp * dao_peak->ext[k];
      fit_data->bg_counts[m] -= 1;
      fit_data->bg_data[m] -= (bg + fit_data->scmos_term[m]);
    }
  }
}


/*
 * daoUpdate()
 *
 * Updates working_peak (integer) center.
 *
 * peak - pointer to a peakData structure.
 */
void daoUpdate(peakData *peak)
{
  daoPeak *dao_peak;

  dao_peak = (daoPeak *)peak->peak_model;

  /* Update peak (integer) center with hysteresis. */
  if(fabs(peak->params[XCENTER] - (double)dao_peak->xc) > HYSTERESIS){
    dao_peak->xc = (int)round(peak->params[XCENTER]);
  }
  if(fabs(peak->params[YCENTER] - (double)dao_peak->yc) > HYSTERESIS){
    dao_peak->yc = (int)round(peak->params[YCENTER]);
  }
}


/*
 * daoUpdate2DFixed()
 *
 * Update for the 2DFixed model.
 *
 * fit_data - pointer to a fitData structure.
 * delta - the deltas for different parameters.
 */
void daoUpdate2DFixed(fitData *fit_data, double *delta)
{
  peakData *peak;

  peak = fit_data->working_peak;

  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], YCENTER);
  mFitUpdateParam(peak, delta[3], BACKGROUND);

  /* Update center position. */
  daoUpdate(peak);

  /* Update ROI size and location. */
  daoCalcLocSize(peak);
}


/*
 * daoUpdate2D()
 *
 * Update for the 2D model.
 *
 * fit_data - pointer to a fitData structure.
 * delta - the deltas for different parameters.
 */
void daoUpdate2D(fitData *fit_data, double *delta)
{
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;
  
  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], YCENTER);
  mFitUpdateParam(peak, delta[3], XWIDTH);
  mFitUpdateParam(peak, delta[4], BACKGROUND);
  peak->params[YWIDTH] = peak->params[XWIDTH];

  /* Update center position. */
  daoUpdate(peak);

  /* Update ROI width. */
  dao_peak->wx = daoCalcWidth(fit_data, peak->params[XWIDTH], dao_peak->wx);
  dao_peak->wy = dao_peak->wx;

  /* Update ROI size and location. */
  daoCalcLocSize(peak);
}


/*
 * daoUpdate3D()
 *
 * Update for the 3D model.
 *
 * fit_data - pointer to a fitData structure.
 * delta - the deltas for different parameters.
 */
void daoUpdate3D(fitData *fit_data, double *delta)
{
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;
  
  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], XWIDTH);
  mFitUpdateParam(peak, delta[3], YCENTER);
  mFitUpdateParam(peak, delta[4], YWIDTH);
  mFitUpdateParam(peak, delta[5], BACKGROUND);

  /* Update center position. */
  daoUpdate(peak);

  /* Update ROI width. */
  dao_peak->wx = daoCalcWidth(fit_data, peak->params[XWIDTH], dao_peak->wx);
  dao_peak->wy = daoCalcWidth(fit_data, peak->params[YWIDTH], dao_peak->wy);

  /* Update ROI size and location. */
  daoCalcLocSize(peak);
}


/*
 * daoUpdateZ()
 *
 * Update for the Z model.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoUpdateZ(fitData *fit_data, double *delta)
{
  peakData *peak;
  daoPeak *dao_peak;
  daoFit *dao_fit;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;
  dao_fit = (daoFit *)fit_data->fit_model;

  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], YCENTER);
  mFitUpdateParam(peak, delta[3], ZCENTER);
  mFitUpdateParam(peak, delta[4], BACKGROUND);

  /* Clamp z value. */
  if(peak->params[ZCENTER] < dao_fit->min_z){
    peak->params[ZCENTER] = dao_fit->min_z;
  }

  if(peak->params[ZCENTER] > dao_fit->max_z){
    peak->params[ZCENTER] = dao_fit->max_z;
  }
    
  /* Update center position. */
  daoUpdate(peak);

  /* Update ROI width. */
  daoCalcWidthsFromZ(fit_data, peak);
  dao_peak->wx = daoCalcWidth(fit_data, peak->params[XWIDTH], dao_peak->wx);
  dao_peak->wy = daoCalcWidth(fit_data, peak->params[YWIDTH], dao_peak->wy);

  /* Update ROI size and location. */
  daoCalcLocSize(peak);
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
