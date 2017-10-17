/*
 * C library for pupil function math.
 *
 * This is designed to work in the context of a fitting algorithm,
 * so basically the steps are:
 *
 * 1. pfInitialize() to set things up.
 * 2. pfSetPF() to set the pupil function.
 * 3. A. pfTranslate() to adjust position in x,y,z.
 *    B. pfGetPSF() to get the PSF at the adjusted position.
 *    C. pfGetPSFdx() to get the derivative of the PSF at the adjusted
 *       position.
 *    D. ...
 *
 * So the PSF will stay at the adjusted position until pfTranslate() 
 * is called again.
 *
 * Note: The boundary conditions are periodic, so the size of the 
 *       pupil function in pixels should be at least 1 pixel larger 
 *       than the PSF.
 *
 * Hazen 10/17.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <fftw3.h>

#include "pupil_function.h"


/*
 * pfnCleanup()
 */
void pfnCleanup(pupilData *pupil_data)
{
  free(pupil_data->kx);
  free(pupil_data->ky);
  free(pupil_data->kz);

  free(pupil_data->kx_c);
  free(pupil_data->kx_r);
  free(pupil_data->ky_c);
  free(pupil_data->ky_r);
  free(pupil_data->kz_c);
  free(pupil_data->kz_r);

  fftw_free(pupil_data->pf);
  fftw_free(pupil_data->ws);
  
  fftw_free(pupil_data->fftw_pf);
  fftw_free(pupil_data->fftw_psf);

  fftw_destroy_plan(pupil_data->fft_backward);
  
  free(pupil_data);
}

/*
 * pfnGetPSF()
 *
 * Get the PSF of the PF.
 */
void pfnGetPSF(pupilData *pupil_data, double *psf_r, double *psf_c)
{
  int i;

  /* Copy current PF into FFTW input. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->fftw_pf[i][0] = pupil_data->ws[i][0];
    pupil_data->fftw_pf[i][1] = pupil_data->ws[i][1];
  }

  /* Perform FFT inverse. */
  fftw_execute(pupil_data->fft_backward);

  /* Return PSF. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    psf_r[i] = pupil_data->fftw_psf[i][0];
    psf_c[i] = pupil_data->fftw_psf[i][1];
  }
}

/*
 * pfnGetPSFdx()
 *
 * Get the derivative of the PSF in x.
 */
void pfnGetPSFdx(pupilData *pupil_data, double *psf_dx_r, double *psf_dx_c)
{
  int i;

  /* Copy current PF multiplied by kx into FFTW input. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->fftw_pf[i][0] = pupil_data->ws[i][1]*pupil_data->kx[i];
    pupil_data->fftw_pf[i][1] = -1.0*pupil_data->ws[i][0]*pupil_data->kx[i];
  }

  /* Perform FFT inverse. */
  fftw_execute(pupil_data->fft_backward);

  /* Return magnitude. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    psf_dx_r[i] = pupil_data->fftw_psf[i][0];
    psf_dx_c[i] = pupil_data->fftw_psf[i][1];
  }
}

/*
 * pfnGetPSFdy()
 *
 * Get the derivative of the PSF in y.
 */
void pfnGetPSFdy(pupilData *pupil_data, double *psf_dy_r, double *psf_dy_c)
{
  int i;

  /* Copy current PF multiplied by ky into FFTW input. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->fftw_pf[i][0] = pupil_data->ws[i][1]*pupil_data->ky[i];
    pupil_data->fftw_pf[i][1] = -1.0*pupil_data->ws[i][0]*pupil_data->ky[i];
  }

  /* Perform FFT inverse. */
  fftw_execute(pupil_data->fft_backward);

  /* Return magnitude. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    psf_dy_r[i] = pupil_data->fftw_psf[i][0];
    psf_dy_c[i] = pupil_data->fftw_psf[i][1];
  }
}

/*
 * pfnGetPSFdz()
 *
 * Get the derivative of the PSF in z.
 */
void pfnGetPSFdz(pupilData *pupil_data, double *psf_dz_r, double *psf_dz_c)
{
  int i;

  /* Copy current PF multiplied by ky into FFTW input. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->fftw_pf[i][0] = pupil_data->ws[i][1]*pupil_data->kz[i];
    pupil_data->fftw_pf[i][1] = -1.0*pupil_data->ws[i][0]*pupil_data->kz[i];
  }

  /* Perform FFT inverse. */
  fftw_execute(pupil_data->fft_backward);

  /* Return magnitude. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    psf_dz_r[i] = pupil_data->fftw_psf[i][0];
    psf_dz_c[i] = pupil_data->fftw_psf[i][1];
  }
}

/*
 * pfnGetSize()
 *
 * Return the X/Y size (in pixels) of the pupil function.
 */
int pfnGetSize(pupilData *pupil_data)
{
  return pupil_data->size;
}

/*
 * pfnInitialize()
 *
 * Initialize pupilData structure. The expectation is that the Python side
 * will provide the values for kx, ky, kz so that the math that this
 * library has to do is minimal.
 */
pupilData *pfnInitialize(double *kx, double *ky, double *kz, int size)
{
  int i;
  pupilData *pupil_data;

  pupil_data = (pupilData *)malloc(sizeof(pupilData));
  pupil_data->size = size;

  pupil_data->kx = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->ky = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->kz = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);

  pupil_data->kx_c = (double *)malloc(sizeof(double)*pupil_data->size);
  pupil_data->kx_r = (double *)malloc(sizeof(double)*pupil_data->size);
  pupil_data->ky_c = (double *)malloc(sizeof(double)*pupil_data->size);
  pupil_data->ky_r = (double *)malloc(sizeof(double)*pupil_data->size);
  pupil_data->kz_c = (double *)malloc(sizeof(double)*(pupil_data->size/2+1)*(pupil_data->size/2+1));
  pupil_data->kz_r = (double *)malloc(sizeof(double)*(pupil_data->size/2+1)*(pupil_data->size/2+1));

  pupil_data->pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->ws = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);

  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->kx[i] = 2.0*M_PI*kx[i];
    pupil_data->ky[i] = 2.0*M_PI*ky[i];
    pupil_data->kz[i] = 2.0*M_PI*kz[i];
  }
  
  pupil_data->fftw_pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->fftw_psf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);

  pupil_data->fft_backward = fftw_plan_dft_2d(size, size, pupil_data->fftw_pf, pupil_data->fftw_psf, FFTW_BACKWARD, FFTW_MEASURE);

  return pupil_data;
}

/*
 * void pfnSetPF()
 *
 * Set/change the pupil function.
 */
void pfnSetPF(pupilData *pupil_data, double *r_pf, double *c_pf)
{
  int i,j,k,l;
  double norm;

  norm = 1.0/((double)pupil_data->size);

  /* 
   * Copy pupil function. 
   * 
   * Notes: 
   *   1. We are messing with the sign so that the results will be
   *      properly centered after the iFFT.
   * 
   *   2. We normalize the PF now so that we only have to do this
   *      once.
   */
  for(i=0;i<pupil_data->size;i++){
    j = i * pupil_data->size;
    for(k=0;k<pupil_data->size;k++){
      l = j+k;
      if(((i+k)%2)==0){
	pupil_data->pf[l][0] = r_pf[l] * norm;
  	pupil_data->pf[l][1] = c_pf[l] * norm;
      }
      else{
	pupil_data->pf[l][0] = -1.0*r_pf[l] * norm;
	pupil_data->pf[l][1] = -1.0*c_pf[l] * norm;
      }
    }
  }

  /*
   * Copy pupil function into working storage so that pfGetPSF() will 
   * work without having to call pfTranslate() first.
   */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->ws[i][0] = pupil_data->pf[i][0];
    pupil_data->ws[i][1] = pupil_data->pf[i][1];
  }
}

/*
 * void pfnTranslate()
 *
 * Translate the PF in x,y,z.
 *
 * X,Y are in units of pixels, Z is in units of microns. Note that dz
 * has the opposite sign from pupil_math.Geometry.changeFocus().
 */
void pfnTranslate(pupilData *pupil_data, double dx, double dy, double dz)
{
  int i,j,l,m,n,o;
  double dd,c1,c2,r1,r2;
  
  /* kx/ky calculations. */
  j = pupil_data->size;
  for(i=0;i<j;i++){
    dd = pupil_data->kx[i*j]*dx;
    pupil_data->kx_r[i] = cos(dd);
    pupil_data->kx_c[i] = -sin(dd);

    dd = pupil_data->kx[i*j]*dy;
    pupil_data->ky_r[i] = cos(dd);
    pupil_data->ky_c[i] = -sin(dd);
  }

  /* 
   * kz calculations. 
   *
   * This is a little complicated, basically kz has radial symmetry centered
   * on pupil_data->size/2 + 1. The idea then is that if we calculate 1/8th 
   * (basically a pie slice) of the values then we have calculated all of the 
   * unique values.
   */
  m = pupil_data->size/2;
  for(i=0;i<=m;i++){
    l = i*(m+1);
    for(j=i;j<=m;j++){
      n = (m-i)*pupil_data->size + (m-j);
      dd = pupil_data->kz[n]*dz;
      pupil_data->kz_r[l+j] = cos(dd);
      pupil_data->kz_c[l+j] = -sin(dd);
      pupil_data->kz_r[j*(m+1)+i] = pupil_data->kz_r[l+j];
      pupil_data->kz_c[j*(m+1)+i] = pupil_data->kz_c[l+j];
    }
  }

  for(i=0;i<pupil_data->size;i++){
    l = i*pupil_data->size;
    n = abs(i-m)*(m+1);
    for(j=0;j<pupil_data->size;j++){
      o = n + abs(j-m);
      r1 = pupil_data->kx_r[i];
      c1 = pupil_data->kx_c[i];

      r2 = r1*pupil_data->ky_r[j] - c1*pupil_data->ky_c[j];
      c2 = r1*pupil_data->ky_c[j] + c1*pupil_data->ky_r[j];

      r1 = r2*pupil_data->kz_r[o] - c2*pupil_data->kz_c[o];
      c1 = r2*pupil_data->kz_c[o] + c2*pupil_data->kz_r[o];
      
      pupil_data->ws[l+j][0] = r1*pupil_data->pf[l+j][0] - c1*pupil_data->pf[l+j][1];
      pupil_data->ws[l+j][1] = r1*pupil_data->pf[l+j][1] + c1*pupil_data->pf[l+j][0];
    }
  }
}
