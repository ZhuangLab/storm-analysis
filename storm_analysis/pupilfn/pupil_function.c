/*
 * C library for pupil function math.
 *
 * We think (and we claimed) that a cubic spline representation of 
 * the PSF is faster than a pupil function representation. Here we 
 * actually test this assumption..
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
  fftw_free(pupil_data->pf);
  fftw_free(pupil_data->kx);
  fftw_free(pupil_data->ky);
  fftw_free(pupil_data->kz);
  
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
    pupil_data->fftw_pf[i][0] = pupil_data->pf[i][0];
    pupil_data->fftw_pf[i][1] = pupil_data->pf[i][1];
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
pupilData *pfInitialize(double *c_kx, double *r_kx, double *c_ky, double *r_ky, double *c_kz, double *r_kz, int size)
{
  int i;
  pupilData *pupil_data;

  pupil_data = (pupilData *)malloc(sizeof(pupilData));
  pupil_data->size = size;

  pupil_data->pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->kx = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->ky = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->kz = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);

  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->kx[i][0] = r_kx[i];
    pupil_data->kx[i][1] = c_kx[i];
    pupil_data->ky[i][0] = r_ky[i];
    pupil_data->ky[i][1] = c_ky[i];
    pupil_data->kz[i][0] = r_kz[i];
    pupil_data->kz[i][1] = c_kz[i];
  }
  
  pupil_data->fftw_pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);
  pupil_data->fftw_psf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pupil_data->size*pupil_data->size);

  pupil_data->fft_backward = fftw_plan_dft_2d(size, size, pupil_data->fftw_pf, pupil_data->fftw_psf, FFTW_BACKWARD, FFTW_MEASURE);

  return pupil_data;
}

/*
 * void pfSetPf()
 *
 * Set/change the pupil function.
 */
void pfSetPf(pupilData *pupil_data, double *r_pf, double *c_pf)
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
}
