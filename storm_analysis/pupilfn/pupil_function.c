/*
 * C library for pupil function math.
 *
 * This is primarily designed to work in the context of a fitting 
 * algorithm, so basically the steps are:
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
 * Note: 
 *  1. The boundary conditions are periodic, so the size of the 
 *     pupil function in pixels should be at least 1 pixel larger 
 *     than the PSF.
 *
 *  2. This library matches the conventions established by
 *     simulator/pupil_math.py, verified by test/test_pupilfn.py.
 *
 *  3. Some additional functionality has been added to provide
 *     faster versions of some calculations that simulator/pupil_math.py
 *     does.
 *
 * Hazen 01/19.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

  if (pupil_data->px_ex != NULL){
    free(pupil_data->px_ex);
    free(pupil_data->px_ey);
    free(pupil_data->py_ex);
    free(pupil_data->py_ey);
    free(pupil_data->pz_ex);
    free(pupil_data->pz_ey);
  }
  
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
 * Get the PSF of the PF. This is complex array.
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
 * pfnGetPSFIntensity()
 *
 * Get the intensity PSF of the PF. Basically this is the result from
 * pfnGetPSF() times it's complex conjugate.
 *
 * This function is not used by pupil_fit.c.
 */
void pfnGetPSFIntensity(pupilData *pupil_data, double *psf_r)
{
  int i;
  double c1,r1;

  /* Copy current PF into FFTW input. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->fftw_pf[i][0] = pupil_data->ws[i][0];
    pupil_data->fftw_pf[i][1] = pupil_data->ws[i][1];
  }

  /* Perform FFT inverse. */
  fftw_execute(pupil_data->fft_backward);

  /* Return PSF. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    r1 = pupil_data->fftw_psf[i][0];
    c1 = pupil_data->fftw_psf[i][1];
    psf_r[i] = r1*r1 + c1*c1;
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
    psf_dz_r[i] = -pupil_data->fftw_psf[i][0];
    psf_dz_c[i] = -pupil_data->fftw_psf[i][1];
  }
}

/*
 * pfnGetPNEN()
 *
 * Get the PSF of the PF multiplied by a dipole pattern. This is complex array.
 *
 * This is not used for PF fitting.
 */
void pfnGetPNEN(pupilData *pupil_data, double *psf_r, double *psf_c, double *dipole)
{
  int i;

  /* Copy current PF times dipole into FFTW input. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->fftw_pf[i][0] = pupil_data->ws[i][0]*dipole[i];
    pupil_data->fftw_pf[i][1] = pupil_data->ws[i][1]*dipole[i];
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
 * pfnGetPXEX()
 *
 * Get the PSF of the PF, X axis detection, X axis emitter.
 *
 * This is not used for PF fitting.
 */
void pfnGetPXEX(pupilData *pupil_data, double *psf_r, double *psf_c)
{
  pfnGetPNEN(pupil_data, psf_r, psf_c, pupil_data->px_ex);
}

/*
 * pfnGetPXEY()
 *
 * Get the PSF of the PF, Y axis detection, X axis emitter.
 *
 * This is not used for PF fitting.
 */
void pfnGetPXEY(pupilData *pupil_data, double *psf_r, double *psf_c)
{
  pfnGetPNEN(pupil_data, psf_r, psf_c, pupil_data->px_ey);
}

/*
 * pfnGetPYEX()
 *
 * Get the PSF of the PF, X axis detection, Y axis emitter.
 *
 * This is not used for PF fitting.
 */
void pfnGetPYEX(pupilData *pupil_data, double *psf_r, double *psf_c)
{
  pfnGetPNEN(pupil_data, psf_r, psf_c, pupil_data->py_ex);
}

/*
 * pfnGetPYEY()
 *
 * Get the PSF of the PF, Y axis detection, Y axis emitter.
 *
 * This is not used for PF fitting.
 */
void pfnGetPYEY(pupilData *pupil_data, double *psf_r, double *psf_c)
{
  pfnGetPNEN(pupil_data, psf_r, psf_c, pupil_data->py_ey);
}

/*
 * pfnGetPZEX()
 *
 * Get the PSF of the PF, X axis detection, Z axis emitter.
 *
 * This is not used for PF fitting.
 */
void pfnGetPZEX(pupilData *pupil_data, double *psf_r, double *psf_c)
{
  pfnGetPNEN(pupil_data, psf_r, psf_c, pupil_data->pz_ex);
}

/*
 * pfnGetPZEY()
 *
 * Get the PSF of the PF, Y axis detection, Z axis emitter.
 *
 * This is not used for PF fitting.
 */
void pfnGetPZEY(pupilData *pupil_data, double *psf_r, double *psf_c)
{
  pfnGetPNEN(pupil_data, psf_r, psf_c, pupil_data->pz_ey);
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

  /* Explicitly NULL these, they might not be used. */
  pupil_data->px_ex = NULL;
  pupil_data->px_ey = NULL;
  pupil_data->py_ex = NULL;
  pupil_data->py_ey = NULL;
  pupil_data->pz_ex = NULL;
  pupil_data->pz_ey = NULL;

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
 * void pfnSetPNEN()
 *
 * Set the dipole radiation patterns for a vectorial PSF.
 *
 * This is not used for PF fitting.
 */
void pfnSetPNEN(pupilData *pupil_data, double *px_ex, double *px_ey, double *py_ex, double *py_ey,double *pz_ex, double *pz_ey)
{
  int i;

  pupil_data->px_ex = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->px_ey = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->py_ex = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->py_ey = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->pz_ex = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);
  pupil_data->pz_ey = (double *)malloc(sizeof(double)*pupil_data->size*pupil_data->size);

  /* Copy dipole patterns. */
  for(i=0;i<(pupil_data->size*pupil_data->size);i++){
    pupil_data->px_ex[i] = px_ex[i];
    pupil_data->px_ey[i] = px_ey[i];
    pupil_data->py_ex[i] = py_ex[i];
    pupil_data->py_ey[i] = py_ey[i];
    pupil_data->pz_ex[i] = pz_ex[i];
    pupil_data->pz_ey[i] = pz_ey[i];
  }
}

/*
 * void pfnTranslate()
 *
 * Translate the PF in x,y,z.
 *
 * X,Y are in units of pixels, Z is in units of microns.
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
   *
   * Also, we use -dz to match the Z convention of simulator.pupil_math.
   */
  m = pupil_data->size/2;
  for(i=0;i<=m;i++){
    l = i*(m+1);
    for(j=i;j<=m;j++){
      n = (m-i)*pupil_data->size + (m-j);
      dd = -pupil_data->kz[n]*dz;
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

/*
 * void pfnTranslateZ()
 *
 * Translate the PF in z only.
 *
 * This function is not used by pupil_fit.c.
 */
void pfnTranslateZ(pupilData *pupil_data, double dz)
{
  int i,j,l,m,n,o;
  double dd,c1,r1;
  
  /* 
   * kz calculations. 
   *
   * This is a little complicated, basically kz has radial symmetry centered
   * on pupil_data->size/2 + 1. The idea then is that if we calculate 1/8th 
   * (basically a pie slice) of the values then we have calculated all of the 
   * unique values.
   *
   * Also, we use -dz to match the Z convention of simulator.pupil_math.
   */
  m = pupil_data->size/2;
  for(i=0;i<=m;i++){
    l = i*(m+1);
    for(j=i;j<=m;j++){
      n = (m-i)*pupil_data->size + (m-j);
      dd = -pupil_data->kz[n]*dz;
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

      r1 = pupil_data->kz_r[o];
      c1 = pupil_data->kz_c[o];
      
      pupil_data->ws[l+j][0] = r1*pupil_data->pf[l+j][0] - c1*pupil_data->pf[l+j][1];
      pupil_data->ws[l+j][1] = r1*pupil_data->pf[l+j][1] + c1*pupil_data->pf[l+j][0];
    }
  }
}
