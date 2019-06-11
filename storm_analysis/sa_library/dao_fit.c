/*
 * Routine(s) for attempting to fit multiple gaussians to an image.
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
 * Calculate 'b' and 'A' for the 2DFixed model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DFixed(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[4];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize b and A. */
  for(i=0;i<4;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<16;i++){
    hessian[i] = 0.0;
  }
  
  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;

    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei;

    /* Calculate b vector. */
    t1 = 2.0*(1.0 - xi/fi);
    for(l=0;l<4;l++){
      jacobian[l] += t1*jt[l];
    }

    /* Calculate A matrix. */
    t2 = 2.0*xi/(fi*fi);
    for(l=0;l<4;l++){
      for(m=l;m<4;m++){
	hessian[l*4+m] += t2*jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJH2DFixedALS()
 *
 * Calculate 'b' and 'A' for the 2DFixed model (Anscombe least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DFixedALS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[4];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize b and A. */
  for(i=0;i<4;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<16;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;
    
    /*
     * This is the LM algorithm according to wikipedia.
     *
     * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
     *
     * fi = 2.0 * (rqe_i * fit_fn_i(theta) + variance_i + 3.0/8.0) ^ 1/2
     * dfi/dtheta = (rqe_i * fit_fn_i(theta) + variance_i + 3.0/8.0) ^ -1/2 * rqe_i * dfit_fn_i(theta)/dtheta
     */
    t1 = fit_data->rqe[k]*2.0/fi;

    jt[0] = t1*e_t;
    jt[1] = t1*2.0*a1*width*xt*e_t;
    jt[2] = t1*2.0*a1*width*yt*e_t;
    jt[3] = t1;
	  
    /* Calculate b vector. */
    t2 = (fi - fit_data->as_xi[k]);
    for(l=0;l<4;l++){
      jacobian[l] += t2*jt[l];
    }

    /* Calculate A matrix. */
    for(l=0;l<4;l++){
      for(m=l;m<4;m++){
	hessian[l*4+m] += jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJH2DFixedLS()
 *
 * Calculate 'b' and 'A' for the 2DFixed model (Least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DFixedLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xt,ext,yt,eyt,e_t,t1,a1,width;
  double jt[4];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize b and A. */
  for(i=0;i<4;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<16;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;
    
    /*
     * This is the LM algorithm according to wikipedia.
     *
     * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
     *
     * fi = rqe_i * fit_fn_i(theta) + variance_i
     * dfi/dtheta = rqe_i * dfit_fn_i(theta)/dtheta
     */
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei;

    /* Calculate b vector. */
    t1 = (fi - fit_data->x_data[k]);
    for(l=0;l<4;l++){
      jacobian[l] += t1*jt[l];
    }

    /* Calculate A matrix. */
    for(l=0;l<4;l++){
      for(m=l;m<4;m++){
	hessian[l*4+m] += jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJH2DFixedDWLS()
 *
 * Calculate 'b' and 'A' for the 2DFixed model (data weighted least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DFixedDWLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,a1,width;
  double jt[4];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize b and A. */
  for(i=0;i<4;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<16;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;

    /*
     * fi = rqe_i * fit_fn_i(theta) + variance_i
     * dfi/dtheta = rqe_i * dfit_fn_i(theta)/dtheta
     */
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei;

    /* Calculate b vector. */
    t1 = (fi - xi)/xi;
    for(l=0;l<4;l++){
      jacobian[l] += t1*jt[l];
    }

    /* Calculate A matrix. */
    for(l=0;l<4;l++){
      for(m=l;m<4;m++){
	hessian[l*4+m] += jt[l]*jt[m]/xi;
      }
    }
  }
}


/* 
 * daoCalcJH2DFixedFWLS()
 *
 * Calculate 'b' and 'A' for the 2DFixed model (fit weighted least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DFixedFWLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[4];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize b and A. */
  for(i=0;i<4;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<16;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;
    
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei;

    /* Calculate b vector. */
    t1 = (fi - xi)/fi;
    for(l=0;l<4;l++){
      jacobian[l] += (2.0*t1 - t1*t1)*jt[l]; 
    }

    /* Calculate A matrix. */
    t2 = 1.0/fi;
    for(l=0;l<4;l++){
      for(m=l;m<4;m++){
	hessian[l*4+m] += 2.0*t2*(xi*t2 - t1 + t1*t1)*jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJH2D()
 *
 * Calculate 'b' and 'A' for the 2D model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2D(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize b and A. */
  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }
  
  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;
    
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei*(-a1*xt*xt*e_t-a1*yt*yt*e_t);
    jt[4] = rqei;
    
    /* Calculate b vector */
    t1 = 2.0*(1.0 - xi/fi);
    for(l=0;l<5;l++){
      jacobian[l] += t1*jt[l];
    }
    
    /* Calculate A matrix. */
    t2 = 2.0*xi/(fi*fi);
    for(l=0;l<5;l++){
      for(m=l;m<5;m++){
	hessian[l*5+m] += t2*jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJH2DALS()
 *
 * Calculate 'b' and 'A' for the 2D model (Anscombe least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DALS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;


  /* Initialize b and A. */
  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;

    /*
     * This is the LM algorithm according to wikipedia.
     *
     * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
     *
     * fi = 2.0 * (rqe_i * fit_fn_i(theta) + variance_i + 3.0/8.0) ^ 1/2
     * dfi/dtheta = (rqe_i * fit_fn_i(theta) + variance_i + 3.0/8.0) ^ -1/2 * rqe_i * dfit_fn_i(theta)/dtheta
     */
    t1 = fit_data->rqe[k]*2.0/fi;
    
    jt[0] = t1*e_t;
    jt[1] = t1*2.0*a1*width*xt*e_t;
    jt[2] = t1*2.0*a1*width*yt*e_t;
    jt[3] = t1*(-a1*xt*xt*e_t-a1*yt*yt*e_t);
    jt[4] = t1;
	  
    /* Calculate b vector. */
    t2 = (fi - fit_data->as_xi[k]);
    for(l=0;l<5;l++){
      jacobian[l] += t2*jt[l];
    }
    
    /* Calculate A matrix. */
    for(l=0;l<5;l++){
      for(m=l;m<5;m++){
	hessian[l*5+m] += jt[l]*jt[m];
      }
    }
      
    /*
     * This is I think more exact, but it takes almost exactly the
     * same number of iterations to converge.
     */
    /*
      jt[0] = e_t;
      jt[1] = 2.0*a1*width*xt*e_t;
      jt[2] = 2.0*a1*width*yt*e_t;
      jt[3] = (-a1*xt*xt*e_t-a1*yt*yt*e_t);
      jt[4] = 1.0;
    */
    /* Calculate b vector. */
    /*
      t1 = -(fit_data->as_xi[k] - fi)*fit_data->rqe[k]/fi;
      for(l=0;l<5;l++){
      jacobian[l] += t1*jt[l];
      }
    */
    /* Calculate A matrix. */
    /*
      t2 = fit_data->as_xi[k]*fit_data->rqe[k]/(2.0*fi*fi*fi);
      for(l=0;l<5;l++){
      for(m=l;m<5;m++){
      hessian[l*5+m] += t2*jt[l]*jt[m];
      }
      }
    */
  }
}


/* 
 * daoCalcJH2DLS()
 *
 * Calculate 'b' and 'A' for the 2D model (Least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xt,ext,yt,eyt,e_t,t1,a1,width;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;


  /* Initialize b and A. */
  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;

    /*
     * This is the LM algorithm according to wikipedia.
     *
     * https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
     *
     * fi = rqe_i * fit_fn_i(theta) + variance_i
     * dfi/dtheta = rqe_i * dfit_fn_i(theta)/dtheta
     */
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei*(-a1*xt*xt*e_t-a1*yt*yt*e_t);
    jt[4] = rqei;
	  
    /* Calculate b vector. */
    t1 = (fi - fit_data->x_data[k]);
    for(l=0;l<5;l++){
      jacobian[l] += t1*jt[l];
    }
    
    /* Calculate A matrix. */
    for(l=0;l<5;l++){
      for(m=l;m<5;m++){
	hessian[l*5+m] += jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJH2DDWLS()
 *
 * Calculate 'b' and 'A' for the 2D model (data weighted least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DDWLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,a1,width;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;


  /* Initialize b and A. */
  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;

    /*
     * fi = rqe_i * fit_fn_i(theta) + variance_i
     * dfi/dtheta = rqe_i * dfit_fn_i(theta)/dtheta
     */
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei*(-a1*xt*xt*e_t-a1*yt*yt*e_t);
    jt[4] = rqei;
	  
    /* Calculate b vector. */
    t1 = (fi - xi)/xi;
    for(l=0;l<5;l++){
      jacobian[l] += t1*jt[l];
    }
    
    /* Calculate A matrix. */
    for(l=0;l<5;l++){
      for(m=l;m<5;m++){
	hessian[l*5+m] += jt[l]*jt[m]/xi;
      }
    }
  }
}


/* 
 * daoCalcJH2DFWLS()
 *
 * Calculate 'b' and 'A' for the 2D model (fit weighted least squares).
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH2DFWLS(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,width;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;


  /* Initialize b and A. */
  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  width = peak->params[XWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;

    /*
     * fi = rqe_i * fit_fn_i(theta) + variance_i
     * dfi/dtheta = rqe_i * dfit_fn_i(theta)/dtheta
     */
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*width*xt*e_t;
    jt[2] = rqei*2.0*a1*width*yt*e_t;
    jt[3] = rqei*(-a1*xt*xt*e_t-a1*yt*yt*e_t);
    jt[4] = rqei;
	  
    /* Calculate b vector. */
    t1 = (fi - fit_data->x_data[k])/fi;
    for(l=0;l<5;l++){
      jacobian[l] += (2.0*t1 - t1*t1)*jt[l];
    }
    
    /* Calculate A matrix. */
    t2 = 1.0/fi;
    for(l=0;l<5;l++){
      for(m=l;m<5;m++){
	hessian[l*5+m] += 2.0*t2*(xi*t2 - t1 + t1*t1)*jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJH3D()
 *
 * Calculate 'b' and 'A' for the 3D model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJH3D(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,a3,a5;
  double jt[6];
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  /* Initialize b and A. */
  for(i=0;i<6;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<36;i++){
    hessian[i] = 0.0;
  }
  
  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  a3 = peak->params[XWIDTH];
  a5 = peak->params[YWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;
    
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*a3*xt*e_t;
    jt[2] = rqei*(-a1*xt*xt*e_t);
    jt[3] = rqei*2.0*a1*a5*yt*e_t;
    jt[4] = rqei*(-a1*yt*yt*e_t);
    jt[5] = rqei;

    /* Calculate b vector. */
    t1 = 2.0*(1.0 - xi/fi);
    for(l=0;l<6;l++){
      jacobian[l] += t1*jt[l];
    }

    /* Calculate A matrix. */
    t2 = 2.0*xi/(fi*fi);
    for(l=0;l<6;l++){
      for(m=l;m<6;m++){
	hessian[l*6+m] += t2*jt[l]*jt[m];
      }
    }
  }
}


/* 
 * daoCalcJHZ()
 *
 * Calculate 'b' and 'A' for the Z model.
 *
 * fit_data - pointer to a fitData structure.
 * jacobian - pointer to an array of double for Jacobian storage. 
 * hessian - pointer to an array of double for Hessian storage. 
 */
void daoCalcJHZ(fitData *fit_data, double *jacobian, double *hessian)
{
  int i,j,k,l,m;
  double fi,rqei,xi,xt,ext,yt,eyt,e_t,t1,t2,a1,a3,a5;
  double z0,z1,z2,zt,gx,gy;
  double jt[5];
  peakData *peak;
  daoPeak *dao_peak;
  daoFit *dao_fit;
  
  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;
  dao_fit = (daoFit *)fit_data->fit_model;

  for(i=0;i<5;i++){
    jacobian[i] = 0.0;
  }
  for(i=0;i<25;i++){
    hessian[i] = 0.0;
  }
    
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

  i = peak->yi * fit_data->image_size_x + peak->xi;
  a1 = peak->params[HEIGHT];
  a3 = peak->params[XWIDTH];
  a5 = peak->params[YWIDTH];
  for(j=0;j<fit_data->roi_n_index;j++){
    k = fit_data->roi_y_index[j]*fit_data->image_size_x + fit_data->roi_x_index[j] + i;
    fi = fit_data->t_fi[k];
    rqei = fit_data->rqe[k];
    xi = fit_data->x_data[k];
    
    l = fit_data->roi_y_index[j];
    yt = dao_peak->yt[l];
    eyt = dao_peak->eyt[l];

    l = fit_data->roi_x_index[j];
    xt = dao_peak->xt[l];
    ext = dao_peak->ext[l];
    
    e_t = ext*eyt;
    
    jt[0] = rqei*e_t;
    jt[1] = rqei*2.0*a1*a3*xt*e_t;
    jt[2] = rqei*2.0*a1*a5*yt*e_t;
    jt[3] = rqei*(-a1*xt*xt*gx*e_t-a1*yt*yt*gy*e_t);
    jt[4] = rqei;
    
    /* Calculate b vector. */
    t1 = 2.0*(1.0 - xi/fi);
    for(l=0;l<5;l++){
      jacobian[l] += t1*jt[l];
    }
    
    /* Calculate A matrix. */
    t2 = 2.0*xi/(fi*fi);
    for(l=0;l<5;l++){
      for(m=l;m<5;m++){
	hessian[l*5+m] += t2*jt[l]*jt[m];
      }
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
  int i,j,k;
  double xo,xt,yo,yt;
  peakData *peak;
  daoPeak *dao_peak;

  peak = fit_data->working_peak;
  dao_peak = (daoPeak *)peak->peak_model;

  xo = -(peak->params[XCENTER] - (double)peak->xi) - fit_data->xoff;
  yo = -(peak->params[YCENTER] - (double)peak->yi) - fit_data->xoff;

  /* Calculate peak shape. */
  for(i=0;i<fit_data->fit_size_x;i++){
    xt = (double)i + xo;
    dao_peak->xt[i] = xt;
    dao_peak->ext[i] = exp(-xt*xt*peak->params[XWIDTH]);
  }
  for(i=0;i<fit_data->fit_size_y;i++){
    yt = (double)i + yo;
    dao_peak->yt[i] = yt;
    dao_peak->eyt[i] = exp(-yt*yt*peak->params[YWIDTH]);
  }

  for(i=0;i<fit_data->roi_n_index;i++){
    j = fit_data->roi_y_index[i];
    k = fit_data->roi_x_index[i];
    peak->psf[i] = dao_peak->eyt[j]*dao_peak->ext[k];
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
  daoFit *dao_fit;

  peak = fit_data->working_peak;
  dao_fit = (daoFit *)fit_data->fit_model;
  
  if(mFitCheck(fit_data)){
    return 1;
  }

  /*
   * Width range clamp.
   *
   * In the previous versions we did not use a width clamp. We removed
   * peaks that had a negative width. However based on some experience
   * with multiplane DAO analysis it appears that using a range clamp
   * works better.
   */
  if(peak->params[XWIDTH] < dao_fit->width_min){
    peak->params[XWIDTH] = dao_fit->width_min;
  }
  else if(peak->params[XWIDTH] > dao_fit->width_max){
    peak->params[XWIDTH] = dao_fit->width_max;
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
  peakData *peak;
  daoFit *dao_fit;

  peak = fit_data->working_peak;
  dao_fit = (daoFit *)fit_data->fit_model;
  
  if(daoCheck2D(fit_data)){
    return 1;
  }

  if(peak->params[YWIDTH] < dao_fit->width_min){
    peak->params[YWIDTH] = dao_fit->width_min;
  }
  else if(peak->params[YWIDTH] > dao_fit->width_max){
    peak->params[YWIDTH] = dao_fit->width_max;
  }
  
  return 0;
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
 * fit_data - pointer to a fitData structure.
 * original - pointer to a peakData structure.
 * copy - pointer to a peakData structure.
 */
void daoCopyPeak(fitData *fit_data, peakData *original, peakData *copy)
{
  daoPeak *dao_copy, *dao_original;

  dao_copy = (daoPeak *)copy->peak_model;
  dao_original = (daoPeak *)original->peak_model;

  /* Allocate storage, if necessary. */
  if(copy->psf == NULL){
    copy->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
    dao_copy->xt = (double *)malloc(sizeof(double)*fit_data->fit_size_x);
    dao_copy->ext = (double *)malloc(sizeof(double)*fit_data->fit_size_x);
    dao_copy->yt = (double *)malloc(sizeof(double)*fit_data->fit_size_y);
    dao_copy->eyt = (double *)malloc(sizeof(double)*fit_data->fit_size_y);
  }
  
  /* This copies the 'core' properties of the structure. */
  mFitCopyPeak(fit_data, original, copy);

  /* Copy the parts that are specific to 3D-DAOSTORM / sCMOS. */
  dao_copy->wx_term = dao_original->wx_term;
  dao_copy->wy_term = dao_original->wy_term;

  memcpy(dao_copy->xt, dao_original->xt, sizeof(double)*fit_data->fit_size_x);
  memcpy(dao_copy->ext, dao_original->ext, sizeof(double)*fit_data->fit_size_x);

  memcpy(dao_copy->yt, dao_original->yt, sizeof(double)*fit_data->fit_size_y);
  memcpy(dao_copy->eyt, dao_original->eyt, sizeof(double)*fit_data->fit_size_y);
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
 * tol - The fitting tolerance.
 * im_size_x - size of the image in x.
 * im_size_y - size of the image in y.
 * roi_size - size of the fitting roi.
 *
 * Returns - Pointer to the fitdata structure.
 */
fitData* daoInitialize(double *rqe, double *scmos_calibration, double tol, int im_size_x, int im_size_y, int roi_size)
{
  daoFit *dao_fit;
  daoPeak* dao_peak;
  fitData* fit_data;

  fit_data = mFitInitialize(rqe, scmos_calibration, tol, im_size_x, im_size_y);
  mFitInitializeROIIndexing(fit_data, roi_size);

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
  dao_fit = (daoFit *)fit_data->fit_model;  

  dao_fit->roi_size = roi_size;
  /*
   * The default is not to do z fitting.
   *
   * We use this flag because it seems easier than having multiple 
   * versions of the daoNewPeaks() function.
   */
  dao_fit->zfit = 0;

  /*
   * Set allowed widths to zero, mostly so that it will be more obvious if we
   * mess up by not setting them to the correct values later.
   */
  dao_fit->width_min = 0.0;
  dao_fit->width_max = 0.0;
  
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
 * fit_data - Pointer to a fitData structure.
 */
void daoInitialize2DFixed(fitData* fit_data)
{
  fit_data->jac_size = 4;

  fit_data->fn_calc_JH = &daoCalcJH2DFixed;
  fit_data->fn_error_fn = &mFitCalcErr;
  fit_data->fn_check = &daoCheck2DFixed;
  fit_data->fn_update = &daoUpdate2DFixed;
}


/*
 * daoInitialize2DFixedALS()
 *
 * Initializes fitting for the 2DFixed model (Anscombe least squares).
 *
 * fit_data - Pointer to a fitData structure.
 */
void daoInitialize2DFixedALS(fitData* fit_data)
{
  fit_data->jac_size = 4;

  fit_data->fn_calc_JH = &daoCalcJH2DFixedALS;
  fit_data->fn_error_fn = &mFitCalcErrALS;  
  fit_data->fn_check = &daoCheck2DFixed;
  fit_data->fn_update = &daoUpdate2DFixed;
}


/*
 * daoInitialize2DFixedLS()
 *
 * Initializes fitting for the 2DFixed model (Least squares).
 *
 * fit_data - Pointer to a fitData structure.
 */
void daoInitialize2DFixedLS(fitData* fit_data)
{
  fit_data->jac_size = 4;

  fit_data->fn_calc_JH = &daoCalcJH2DFixedLS;
  fit_data->fn_error_fn = &mFitCalcErrLS;  
  fit_data->fn_check = &daoCheck2DFixed;
  fit_data->fn_update = &daoUpdate2DFixed;
}


/*
 * daoInitialize2DFixedDWLS()
 *
 * Initializes fitting for the 2DFixed model (Data weighted least squares).
 *
 * fit_data - Pointer to a fitData structure.
 */
void daoInitialize2DFixedDWLS(fitData* fit_data)
{
  fit_data->jac_size = 4;

  fit_data->fn_calc_JH = &daoCalcJH2DFixedDWLS;
  fit_data->fn_error_fn = &mFitCalcErrDWLS;
  fit_data->fn_check = &daoCheck2DFixed;
  fit_data->fn_update = &daoUpdate2DFixed;
}


/*
 * daoInitialize2DFixedFWLS()
 *
 * Initializes fitting for the 2DFixed model (Fit weighted least squares).
 *
 * fit_data - Pointer to a fitData structure.
 */
void daoInitialize2DFixedFWLS(fitData* fit_data)
{
  fit_data->jac_size = 4;

  fit_data->fn_calc_JH = &daoCalcJH2DFixedFWLS;
  fit_data->fn_error_fn = &mFitCalcErrFWLS;
  fit_data->fn_check = &daoCheck2DFixed;
  fit_data->fn_update = &daoUpdate2DFixed;
}


/*
 * daoInitialize2D()
 *
 * Initializes fitting for the 2D model.
 *
 * fit_data - Pointer to a fitData structure.
 * width_min - Minimum allowed peak width.
 * width_max - Maximum allowed peak width.
 */
void daoInitialize2D(fitData* fit_data, double width_min, double width_max)
{
  fit_data->jac_size = 5;

  ((daoFit *)fit_data->fit_model)->width_min = width_min;
  ((daoFit *)fit_data->fit_model)->width_max = width_max;

  fit_data->fn_calc_JH = &daoCalcJH2D;
  fit_data->fn_error_fn = &mFitCalcErr;
  fit_data->fn_check = &daoCheck2D;
  fit_data->fn_update = &daoUpdate2D;
}


/*
 * daoInitialize2DALS()
 *
 * Initializes fitting for the 2D model (Anscombe least squares).
 *
 * fit_data - Pointer to a fitData structure.
 * width_min - Minimum allowed peak width.
 * width_max - Maximum allowed peak width.
 */
void daoInitialize2DALS(fitData* fit_data, double width_min, double width_max)
{
  fit_data->jac_size = 5;

  ((daoFit *)fit_data->fit_model)->width_min = width_min;
  ((daoFit *)fit_data->fit_model)->width_max = width_max;

  fit_data->fn_calc_JH = &daoCalcJH2DALS;
  fit_data->fn_error_fn = &mFitCalcErrALS;
  fit_data->fn_check = &daoCheck2D;
  fit_data->fn_update = &daoUpdate2D;
}


/*
 * daoInitialize2DLS()
 *
 * Initializes fitting for the 2D model (Least squares).
 *
 * fit_data - Pointer to a fitData structure.
 * width_min - Minimum allowed peak width.
 * width_max - Maximum allowed peak width.
 */
void daoInitialize2DLS(fitData* fit_data, double width_min, double width_max)
{
  fit_data->jac_size = 5;

  ((daoFit *)fit_data->fit_model)->width_min = width_min;
  ((daoFit *)fit_data->fit_model)->width_max = width_max;

  fit_data->fn_calc_JH = &daoCalcJH2DLS;
  fit_data->fn_error_fn = &mFitCalcErrLS;
  fit_data->fn_check = &daoCheck2D;
  fit_data->fn_update = &daoUpdate2D;
}


/*
 * daoInitialize2DDWLS()
 *
 * Initializes fitting for the 2D model (Data weighted least squares).
 *
 * fit_data - Pointer to a fitData structure.
 * width_min - Minimum allowed peak width.
 * width_max - Maximum allowed peak width.
 */
void daoInitialize2DDWLS(fitData* fit_data, double width_min, double width_max)
{
  fit_data->jac_size = 5;

  ((daoFit *)fit_data->fit_model)->width_min = width_min;
  ((daoFit *)fit_data->fit_model)->width_max = width_max;

  fit_data->fn_calc_JH = &daoCalcJH2DDWLS;
  fit_data->fn_error_fn = &mFitCalcErrDWLS;
  fit_data->fn_check = &daoCheck2D;
  fit_data->fn_update = &daoUpdate2D;
}


/*
 * daoInitialize2DFWLS()
 *
 * Initializes fitting for the 2D model (Fit weighted least squares).
 *
 * fit_data - Pointer to a fitData structure.
 * width_min - Minimum allowed peak width.
 * width_max - Maximum allowed peak width.
 */
void daoInitialize2DFWLS(fitData* fit_data, double width_min, double width_max)
{
  fit_data->jac_size = 5;

  ((daoFit *)fit_data->fit_model)->width_min = width_min;
  ((daoFit *)fit_data->fit_model)->width_max = width_max;

  fit_data->fn_calc_JH = &daoCalcJH2DFWLS;
  fit_data->fn_error_fn = &mFitCalcErrFWLS;
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
void daoInitialize3D(fitData* fit_data, double width_min, double width_max)
{
  fit_data->jac_size = 6;

  ((daoFit *)fit_data->fit_model)->width_min = width_min;
  ((daoFit *)fit_data->fit_model)->width_max = width_max;
  
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
  int i,j;
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
      
      /* Allocate space for peak psf and other arrays. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
	dao_peak->xt = (double *)malloc(sizeof(double)*fit_data->fit_size_x);
	dao_peak->ext = (double *)malloc(sizeof(double)*fit_data->fit_size_x);
	dao_peak->yt = (double *)malloc(sizeof(double)*fit_data->fit_size_y);
	dao_peak->eyt = (double *)malloc(sizeof(double)*fit_data->fit_size_y);
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

      /* Arbitrary initial values for BACKGROUND, HEIGHT. */
      peak->params[BACKGROUND] = 1.0;
      peak->params[HEIGHT] = 1.0;
      
      /* Copy into working peak. */
      daoCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	daoCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }

      /* Calculate peak shape (of working peak). */
      daoCalcPeakShape(fit_data);

      /* Estimate best starting background. */
      mFitEstimatePeakBackground(fit_data);

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
	  daoCopyPeak(fit_data, fit_data->working_peak, peak);
	  continue;
	}
      }
      
      /* Add peak to the fit image. */
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      daoCopyPeak(fit_data, fit_data->working_peak, peak);
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

      /* Allocate space for peak psf and other arrays. */
      if(peak->psf == NULL){
	peak->psf = (double *)malloc(sizeof(double)*fit_data->fit_size_x*fit_data->fit_size_y);
	dao_peak->xt = (double *)malloc(sizeof(double)*fit_data->fit_size_x);
	dao_peak->ext = (double *)malloc(sizeof(double)*fit_data->fit_size_x);
	dao_peak->yt = (double *)malloc(sizeof(double)*fit_data->fit_size_y);
	dao_peak->eyt = (double *)malloc(sizeof(double)*fit_data->fit_size_y);
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
      daoCopyPeak(fit_data, peak, fit_data->working_peak);

      /* Check that the peak is okay. */
      if(fit_data->fn_check(fit_data)){
	if(TESTING){
	  printf("Warning peak %d is bad!\n", (i-start));
	}
	fit_data->working_peak->status = ERROR;
	daoCopyPeak(fit_data, fit_data->working_peak, peak);
	continue;
      }
      
      /* Add peak to the fit image. */
      daoCalcPeakShape(fit_data);
      mFitAddPeak(fit_data);

      /* Copy values back from working peak. */
      daoCopyPeak(fit_data, fit_data->working_peak, peak);
    }
  }
  else{
    printf("Unrecognized peak type '%s'!\n", p_type);
  }

  /* Initial error calculation. */
  for(i=start;i<stop;i++){
    peak = &fit_data->fit[i];
    daoCopyPeak(fit_data, peak, fit_data->working_peak);
    fit_data->fn_error_fn(fit_data);
    daoCopyPeak(fit_data, fit_data->working_peak, peak);
  }

  fit_data->nfit = stop;
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
