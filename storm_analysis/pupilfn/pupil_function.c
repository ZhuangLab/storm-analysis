/*
 * C library for pupil function math.
 *
 * We think (and we claimed) that a cubic spline representation of 
 * the PSF is faster than a pupil function representation. Here we 
 * actually test this assumption..
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
 * Hazen 10/17.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <fftw3.h>

#include "pupil_function.h"


/*
 * pfCleanup()
 */
void pfCleanup(pupilData *pupil_data)
{
  free(pupil_data->kx);
  free(pupil_data->ky);
  free(pupil_data->kz);

  fftw_free(pupil_data->pf);
  fftw_free(pupil_data->ws);
  
  fftw_free(pupil_data->fftw_pf);
  fftw_free(pupil_data->fftw_psf);

  fftw_destroy_plan(pupil_data->fft_backward);
  
  free(pupil_data);
}

/*
 * pfGetPSF()
 *
 * Get the PSF of the PF.
 */
void pfGetPSF(pupilData *pupil_data, double *psf)
{
  int i;

  /* Copy current PF into FFTW input. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->fftw_pf[i][0] = pupil_data->ws[i][0];
    pupil_data->fftw_pf[i][1] = pupil_data->ws[i][1];
  }

  /* Perform FFT inverse. */
  fftw_execute(pupil_data->fft_backward);

  /* Return magnitude. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    psf[i] = pupil_data->fftw_psf[i][0]*pupil_data->fftw_psf[i][0] + pupil_data->fftw_psf[i][1]*pupil_data->fftw_psf[i][1];
  }
}

/*
 * pfInitialize()
 *
 * Initialize pupilData structure. The expectation is that the Python side
 * will provide the values for kx, ky, kz so that the math that this
 * library has to do is minimal.
 */
pupilData *pfInitialize(double *kx, double *ky, double *kz, int size)
{
  int i;
  pupilData *pupil_data;

  pupil_data = (pupilData *)malloc(sizeof(pupilData));
  pupil_data->size = size;

  pupil_data->kx = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->ky = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->kz = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);

  pupil_data->pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->ws = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);

  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->kx[i] = kx[i];
    pupil_data->ky[i] = ky[i];
    pupil_data->kz[i] = kz[i];
  }
  
  pupil_data->fftw_pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->fftw_psf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);

  pupil_data->fft_backward = fftw_plan_dft_2d(size, size, pupil_data->fftw_pf, pupil_data->fftw_psf, FFTW_BACKWARD, FFTW_MEASURE);

  return pupil_data;
}

/*
 * void pfSetPF()
 *
 * Set/change the pupil function.
 */
void pfSetPF(pupilData *pupil_data, double *r_pf, double *c_pf)
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
 * void pfTranslate()
 *
 * Translate the PF in x,y,z.
 *
 * X,Y are in units of pixels, Z is in units of microns. Note that dz
 * has the opposite sign from pupil_math.Geometry.changeFocus().
 */
void pfTranslate(pupilData *pupil_data, double dx, double dy, double dz)
{
  int i;
  double dd, dd_c, dd_r;

  for(i=0;i<(pupil_data->size*pupil_data->size);i++){

    dd = 2.0*M_PI*(pupil_data->kx[i]*dx + pupil_data->ky[i]*dy + pupil_data->kz[i]*dz);
    dd_r = cos(dd);
    dd_c = -sin(dd);

    pupil_data->ws[i][0] = dd_r*pupil_data->pf[i][0] - dd_c*pupil_data->pf[i][1];
    pupil_data->ws[i][1] = dd_r*pupil_data->pf[i][1] + dd_c*pupil_data->pf[i][0];
  }
}
