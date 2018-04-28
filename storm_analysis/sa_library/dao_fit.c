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

#include "dao_fit.h"


/*
 * daoAllocPeaks()
 *
 * Allocate storage for daoPeaks.
 */
void daoAllocPeaks(peakData *new_peaks, int n_peaks)
{
  int i;
  daoPeak *dao_peak;

  for(i=0;i<n_peaks;i++){
    new_peaks[i].peak_model = (daoPeak *)malloc(sizeof(daoPeak));
    dao_peak = (daoPeak *)new_peaks[i].peak_model;
    dao_peak->xt = NULL;
    dao_peak->ext = NULL;
    dao_peak->yt = NULL;
    dao_peak->eyt = NULL;
  }
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
 * daoCalcJH2DLS()
 *
 * Calculate Jacobian and Hessian for the 2D model (least-squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int j,k,l,m;
  double fi,xi,xt,ext,yt,eyt,e_t,t1,a1,width;
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
      t1 = (fi - xi);
      jacobian[0] += t1*jt[0];
      jacobian[1] += t1*jt[1];
      jacobian[2] += t1*jt[2];
      jacobian[3] += t1*jt[3];
      jacobian[4] += t1*jt[4];
      
      // calculate hessian.
      hessian[0] += jt[0]*jt[0];
      hessian[1] += jt[0]*jt[1];
      hessian[2] += jt[0]*jt[2];
      hessian[3] += jt[0]*jt[3];
      hessian[4] += jt[0]*jt[4];
	  
      // hessian[5]
      hessian[6] += jt[1]*jt[1];
      hessian[7] += jt[1]*jt[2];
      hessian[8] += jt[1]*jt[3];
      hessian[9] += jt[1]*jt[4];
	    
      // hessian[10]
      // hessian[11]
      hessian[12] += jt[2]*jt[2];
      hessian[13] += jt[2]*jt[3];
      hessian[14] += jt[2]*jt[4];
	  
      // hessian[15]
      // hessian[16]
      // hessian[17]
      hessian[18] += jt[3]*jt[3];
      hessian[19] += jt[3]*jt[4];
	
      // hessian[20]
      // hessian[21]
      // hessian[22]
      // hessian[23]
      hessian[24] += jt[4]*jt[4];
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
 * daoCalcPeakShape()
 *
 * Calculate peak shape.
 *
 * fit_data - pointer to a fitData structure.
 */
void daoCalcPeakShape(fitData *fit_data)
{
  int j,k,l;
  double xo,xt,yo,yt;
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  xo = -(peak->params[XCENTER] - (double)peak->xi) - fit_data->xoff;
  yo = -(peak->params[YCENTER] - (double)peak->yi) - fit_data->xoff;

  /* Calculate peak shape. */
  for(j=0;j<peak->size_x;j++){
    xt = (double)j + xo;
    dao_peak->xt[j] = xt;
    dao_peak->ext[j] = exp(-xt*xt*peak->params[XWIDTH]);
  }
  for(j=0;j<peak->size_y;j++){
    yt = (double)j + yo;
    dao_peak->yt[j] = yt;
    dao_peak->eyt[j] = exp(-yt*yt*peak->params[YWIDTH]);
  }

  for(j=0;j<peak->size_y;j++){
    l = j*peak->size_x;
    for(k=0;k<peak->size_x;k++){
      peak->psf[l+k] = dao_peak->ext[k]*dao_peak->eyt[j];
    }
  }
}


/*
 * daoCalcWidthsFromZ()
 *
 * Updates peak widths given z.
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
 * daoCheck2DFixed()
 *
 * Check that the parameters of working_peak are still valid (2DFixed model).
 *
 * fit_data - pointer to a fitData structure.
 *
 * Return 0 if okay.
 */
int daoCheck2DFixed(fitData *fit_data){
  return mFitCheck(fit_data);
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
  
  if(mFitCheck(fit_data)){
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
  return mFitCheck(fit_data);
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
  daoPeak *dao_peak;

  /* Free fitting specific peaks data of each peak. */
  if(fit_data->fit != NULL){
    for(i=0;i<fit_data->max_nfit;i++){
      dao_peak = (daoPeak *)fit_data->fit[i].peak_model;
      free(dao_peak->xt);
      free(dao_peak->ext);
      free(dao_peak->yt);
      free(dao_peak->eyt);
      free(dao_peak);
    }
  }

  /* Free working peak. */
  dao_peak = (daoPeak *)fit_data->working_peak->peak_model;
  free(dao_peak->xt);
  free(dao_peak->ext);
  free(dao_peak->yt);
  free(dao_peak->eyt);
  free(dao_peak);
  
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

  /* Allocate storage, if necessary. */
  if(copy->psf == NULL){
    copy->psf = (double *)malloc(sizeof(double)*original->size_x*original->size_y);
    dao_copy->xt = (double *)malloc(sizeof(double)*original->size_x);
    dao_copy->ext = (double *)malloc(sizeof(double)*original->size_x);
    dao_copy->yt = (double *)malloc(sizeof(double)*original->size_y);
    dao_copy->eyt = (double *)malloc(sizeof(double)*original->size_y);
  }
  
  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(original, copy);

  /* Copy the parts that are specific to 3D-DAOSTORM / sCMOS. */
  dao_copy->wx_term = dao_original->wx_term;
  dao_copy->wy_term = dao_original->wy_term;
  
  for(i=0;i<original->size_x;i++){
    dao_copy->xt[i] = dao_original->xt[i];
    dao_copy->ext[i] = dao_original->ext[i];
  }

  for(i=0;i<original->size_y;i++){
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
  daoPeak *dao_peak;

  for(i=0;i<n_peaks;i++){
    dao_peak = (daoPeak *)peaks[i].peak_model;
    free(dao_peak->xt);
    free(dao_peak->ext);
    free(dao_peak->yt);
    free(dao_peak->eyt);
    free(dao_peak); 
  }
}


/*
 * daoInitialize()
 *
 * Initializes fitting things for fitting.
 *
 * rqe - Pixel relative quantum efficiency.
 * scmos_calibration - sCMOS calibration data, variance/gain^2 for each pixel in the image.
 * clamp - The starting clamp values for each peak.
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 * roi_size - size of the fitting roi.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* daoInitialize(double *rqe, double *scmos_calibration, double *clamp, double tol, int im_size_x, int im_size_y, int roi_size)
{
  daoPeak* dao_peak;
  fitData* fit_data;

  fit_data = mFitInitialize(rqe, scmos_calibration, clamp, tol, im_size_x, im_size_y);

  /*
   * I thought that even/odd might need different offsets, but based on some testing
   * it looks like the optimal offset is pretty similar / the same for both.
   */
  if((roi_size%2)==0){
    fit_data->xoff = 0.5*((double)roi_size) - 1.0;
    fit_data->yoff = 0.5*((double)roi_size) - 1.0;
  }
  else{
    fit_data->xoff = 0.5*((double)roi_size) - 1.0;
    fit_data->yoff = 0.5*((double)roi_size) - 1.0;
  }

  fit_data->fit_model = (daoFit *)malloc(sizeof(daoFit));

  ((daoFit *)fit_data->fit_model)->roi_size = roi_size;
  /*
   * The default is not to do z fitting.
   *
   * We use this flag because it seems easier than having multiple 
   * versions of the daoNewPeaks() function.
   */
  ((daoFit *)fit_data->fit_model)->zfit = 0;
  
  /* Allocate storage for the working peak. */
  fit_data->working_peak->psf = (double *)malloc(sizeof(double)*roi_size*roi_size);
    
  fit_data->working_peak->peak_model = (daoPeak *)malloc(sizeof(daoPeak));
  dao_peak = (daoPeak *)fit_data->working_peak->peak_model;
  dao_peak->xt = (double *)malloc(sizeof(double)*roi_size);
  dao_peak->ext = (double *)malloc(sizeof(double)*roi_size);
  dao_peak->yt = (double *)malloc(sizeof(double)*roi_size);
  dao_peak->eyt = (double *)malloc(sizeof(double)*roi_size);

  /* Set function pointers. */
  fit_data->fn_alloc_peaks = &daoAllocPeaks;
  fit_data->fn_calc_peak_shape = &daoCalcPeakShape;
  fit_data->fn_copy_peak = &daoCopyPeak;
  fit_data->fn_free_peaks = &daoFreePeaks;
  
  return fit_data;
}


/*
 * daoInitialize2DFixed()
 *
 * Initializes fitting for the 2DFixed model.
 *
 * fit_data - pointer to a fitData structure.
 * peak_size - the peak roi size in pixels.
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
void daoInitialize2D(fitData* fit_data, int ls_fit)
{
  fit_data->jac_size = 5;

  if(ls_fit == 0){
    fit_data->fn_calc_JH = &daoCalcJH2D;
    fit_data->fn_error_fn = &mFitCalcErr;
  }
  else{
    fit_data->fn_calc_JH = &daoCalcJH2DLS;
    fit_data->fn_error_fn = &mFitCalcErrLS;
  }
  
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
  int i,j,xc,yc;
  int start,stop;
  double width;
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

      /*
       * We used to allow variable peak sizes. Now this is mostly because
       * it is more convenient for the peak to know it's own size instead
       * of having to refer back to the fit_data structure.
       */
      peak->size_x = ((daoFit *)fit_data->fit_model)->roi_size;
      peak->size_y = ((daoFit *)fit_data->fit_model)->roi_size;
      
      /* Allocate space for peak psf and other arrays. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*peak->size_x*peak->size_y);
	dao_peak->xt = (double *)malloc(sizeof(double)*peak->size_x);
	dao_peak->ext = (double *)malloc(sizeof(double)*peak->size_x);
	dao_peak->yt = (double *)malloc(sizeof(double)*peak->size_y);
	dao_peak->eyt = (double *)malloc(sizeof(double)*peak->size_y);
      }
      
      /* Initial width. */
      if(((daoFit *)fit_data->fit_model)->zfit){
	daoCalcWidthsFromZ(fit_data, peak);
      }
      else{
	width = 1.0/(2.0*peak_params[j+3]*peak_params[j+3]);
	peak->params[XWIDTH] = width;
	peak->params[YWIDTH] = width;
      }

      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);

      /* Estimate background. */
      xc = (int)round(peak_params[j]);
      yc = (int)round(peak_params[j+1]);      
      peak->params[BACKGROUND] = fit_data->bg_estimate[yc * fit_data->image_size_x + xc];

      /* Arbitrary initial value for HEIGHT. */
      peak->params[HEIGHT] = 1.0;
      
      /* Copy into working peak. */
      daoCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	daoCopyPeak(fit_data->working_peak, peak);
	continue;
      }

      /* Calculate peak shape (of working peak). */
      daoCalcPeakShape(fit_data);

      if(!strcmp(p_type, "finder")){

	/* Estimate best starting height. */
	mFitEstimatePeakHeight(fit_data);

	/* 
	 * This is for the benefit of multi-plane fitting where some peaks might
	 * not be present in a particular plane, but we don't want to throw them
	 * out here as they presumably are present in other planes.
	 */
	if(fit_data->working_peak->params[HEIGHT] < fit_data->minimum_height){
	  fit_data->working_peak->params[HEIGHT] = fit_data->minimum_height;
	}
	
	/* Check that the initial height is positive, error it out if not. */
	if(fit_data->working_peak->params[HEIGHT] <= 0.0){
	  printf("Warning peak %d has negative estimated height!\n", (i-start));
	  fit_data->working_peak->status = ERROR;
	  daoCopyPeak(fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      daoCopyPeak(fit_data->working_peak, peak);
    }
  }
  /*
   * "pre-specified" parameters, these are the peak x,y,z and sigma 
   * values as an n_peaks x 7 array.
   */
  else if(!strcmp(p_type, "text") || !strcmp(p_type, "hdf5")){
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

      /* Set fitting size. */
      peak->size_x = ((daoFit *)fit_data->fit_model)->roi_size;
      peak->size_y = ((daoFit *)fit_data->fit_model)->roi_size;

      /* Allocate space for peak psf and other arrays. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*peak->size_x*peak->size_y);
	dao_peak->xt = (double *)malloc(sizeof(double)*peak->size_x);
	dao_peak->ext = (double *)malloc(sizeof(double)*peak->size_x);
	dao_peak->yt = (double *)malloc(sizeof(double)*peak->size_y);
	dao_peak->eyt = (double *)malloc(sizeof(double)*peak->size_y);
      }

      /* Initial width. */
      width = 1.0/(2.0*peak_params[j+5]*peak_params[j+5]);
      peak->params[XWIDTH] = width;

      width = 1.0/(2.0*peak_params[j+6]*peak_params[j+6]);
      peak->params[YWIDTH] = width;

      /* Correct initial width (z fitting model). */
      if(((daoFit *)fit_data->fit_model)->zfit){
	daoCalcWidthsFromZ(fit_data, peak);
      }

      peak->xi = (int)floor(peak->params[XCENTER]);
      peak->yi = (int)floor(peak->params[YCENTER]);
      
      /* Copy into working peak. */
      daoCopyPeak(peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	daoCopyPeak(fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      daoCalcPeakShape(fit_data);
      mFitAddPeak(fit_data);

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
    fit_data->fn_error_fn(fit_data);
    daoCopyPeak(fit_data->working_peak, peak);
  }

  fit_data->nfit = stop;

  /* Reset the clamp values on all the peaks. */
  if(USECLAMP){
    mFitResetClampValues(fit_data);
  }
}


/*
 * daoPeakSum()
 *
 * Return the integral of the PSF.
 */
double daoPeakSum(peakData *peak)
{
  int i,j;
  double sum,tmp;
  daoPeak *dao_peak;

  dao_peak = (daoPeak *)peak->peak_model;

  sum = 0.0;
  for(i=0;i<peak->size_y;i++){
    tmp = dao_peak->eyt[i];
    for(j=0;j<peak->size_x;j++){
      sum += tmp * dao_peak->ext[j];
    }
  }
  sum = sum*peak->params[HEIGHT];

  return sum;
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
  mFitUpdate(peak);
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

  peak = fit_data->working_peak;
  
  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], YCENTER);
  mFitUpdateParam(peak, delta[3], XWIDTH);
  mFitUpdateParam(peak, delta[4], BACKGROUND);
  peak->params[YWIDTH] = peak->params[XWIDTH];

  /* Update center position. */
  mFitUpdate(peak);
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

  peak = fit_data->working_peak;
  
  mFitUpdateParam(peak, delta[0], HEIGHT);
  mFitUpdateParam(peak, delta[1], XCENTER);
  mFitUpdateParam(peak, delta[2], XWIDTH);
  mFitUpdateParam(peak, delta[3], YCENTER);
  mFitUpdateParam(peak, delta[4], YWIDTH);
  mFitUpdateParam(peak, delta[5], BACKGROUND);

  /* Update center position. */
  mFitUpdate(peak);
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
  daoFit *dao_fit;

  peak = fit_data->working_peak;
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
  mFitUpdate(peak);

  /* Update ROI width. */
  daoCalcWidthsFromZ(fit_data, peak);
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
