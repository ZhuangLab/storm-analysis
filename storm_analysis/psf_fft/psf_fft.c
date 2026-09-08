/*
 * C library for PSF manipulation using the FFT.
 *
 * This is designed to work in the context of a fitting algorithm,
 * so basically the steps are:
 *
 * 1. pFTInitialize() to set things up.
 * 2. A. pFTTranslate() to adjust position in x,y,z.
 *    B. pFTGetPSF() to get the PSF at the adjusted position.
 *    C. pFTGetPSFdx() to get the derivative of the PSF at the adjusted
 *       position.
 *    D. ...
 *
 * So the PSF will stay at the adjusted position until pFTTranslate() 
 * is called again.
 *
 * Note: The boundary conditions are periodic, so the size of the 
 *       PSF should be large enough that it goes to zero at the edges.
 *
 * Hazen 10/17.
 */


/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "psf_fft.h"

/* Functions */

/*
 * pFTCalcShiftVector()
 *
 * Calculate the FFT shift vector.
 *
 * sr - The real part of the shift vector.
 * sc - The complex part of the shift vector.
 * dx - Shift delta (in pixels).
 * size - The size of the vector.
 */
void pFTCalcShiftVector(double *sr, double *sc, double dx, int size)
{
  int i;
  double t1, t2, tc, tr;

  sr[0] = 1.0;
  sc[0] = 0.0;
  t1 = -M_PI * 2.0 * dx/((double)size);
  for(i=1;i<(size/2+1);i++){
    t2 = t1 * (double)i;
    tr = cos(t2);
    tc = sin(t2);
    sr[size-i] = tr;
    sr[i] = tr;
    sc[size-i] = -tc;
    sc[i] = tc;
  }
}

/*
 * pFTCalcShiftVectorDerivative()
 *
 * Calculate the derivative of the FFT shift vector.
 *
 * sc - The complex part of the derivative of the shift vector.
 * size - The size of the vector.
 */
void pFTCalcShiftVectorDerivative(double *sc, int size)
{
  int i;
  double t1,t2;

  sc[0] = 0.0;
  t1 = -M_PI * 2.0/((double)size);
  for(i=1;i<(size/2+1);i++){
    t2 = t1 * (double)i;
    sc[size-i] = t2;
    sc[i] = -t2;
  }
}

/*
 * pFTCleanup()
 *
 * pfft - A pointer to a psfFFT structure.
 */
void pFTCleanup(psfFFT *pfft)
{
  free(pfft->kx_c);
  free(pfft->kx_r);
  free(pfft->ky_c);
  free(pfft->ky_r);
  free(pfft->kz_c);
  free(pfft->kz_r);
  
  fftw_free(pfft->fftw_real);
  fftw_free(pfft->fftw_fft);
  fftw_free(pfft->ws);
  fftw_free(pfft->psf);
  
  fftw_destroy_plan(pfft->fft_backward);

  free(pfft);
}

/*
 * pFTGetPSF()
 *
 * Return the current PSF.
 *
 * pfft - A pointer to a psfFFT structure.
 * psf - Pre-allocated storage for the result.
 */
void pFTGetPSF(psfFFT *pfft, double *psf)
{
  int i;
  int mid_z, size_xy;

  /* Copy current FFT of PSF into FFTW input, */
  for(i=0;i<pfft->fft_size;i++){
    pfft->fftw_fft[i][0] = pfft->ws[i][0];
    pfft->fftw_fft[i][1] = pfft->ws[i][1];
  }
  
  /* Do reverse transform. */
  fftw_execute(pfft->fft_backward);

  /* The 2D PSF is the middle plane of the 3D PSF. */
  size_xy = pfft->x_size*pfft->y_size;
  mid_z = size_xy * (pfft->z_size/2);
  
  for(i=0;i<size_xy;i++){
    psf[i] = pfft->fftw_real[mid_z+i];
  }
}

/*
 * pFTGetPSFdx()
 *
 * Return the derivative with respect to x of the current PSF.
 *
 * pfft - A pointer to a psfFFT structure.
 * dx - Pre-allocated storage for the result.
 */
void pFTGetPSFdx(psfFFT *pfft, double *dx)
{
  int i,j,k,t1,t2,t3;
  int mid_z,size_xy;
  double *dd;

  /* Calculate FFT derivative vector. */
  dd = pfft->kx_c;
  pFTCalcShiftVectorDerivative(dd, pfft->x_size);

  /* Copy current PF multiplied by kx into FFTW input. */
  for(i=0;i<pfft->z_size;i++){
    t1 = i * (pfft->y_size * pfft->fft_x_size);
    for(j=0;j<pfft->y_size;j++){
      t2 = j * pfft->fft_x_size;
      for(k=0;k<pfft->fft_x_size;k++){
	t3 = t1 + t2 + k;

	pfft->fftw_fft[t3][0] = pfft->ws[t3][1]*dd[k];
	pfft->fftw_fft[t3][1] = -1.0*pfft->ws[t3][0]*dd[k];
      }
    }
  }
  
  /* Do reverse transform. */
  fftw_execute(pfft->fft_backward);

  /* The 2D PSF is the middle plane of the 3D PSF. */
  size_xy = pfft->x_size*pfft->y_size;
  mid_z = size_xy * (pfft->z_size/2);
  
  for(i=0;i<size_xy;i++){
    dx[i] = pfft->fftw_real[mid_z+i];
  }
}

/*
 * pFTGetPSFdy()
 *
 * Return the derivative with respect to y of the current PSF.
 *
 * pfft - A pointer to a psfFFT structure.
 * dy - Pre-allocated storage for the result.
 */
void pFTGetPSFdy(psfFFT *pfft, double *dy)
{
  int i,j,k,t1,t2,t3;
  int mid_z,size_xy;
  double *dd;

  /* Calculate FFT derivative vector. */
  dd = pfft->ky_c;
  pFTCalcShiftVectorDerivative(dd, pfft->y_size);

  /* Copy current PF multiplied by ky into FFTW input. */
  for(i=0;i<pfft->z_size;i++){
    t1 = i * (pfft->y_size * pfft->fft_x_size);
    for(j=0;j<pfft->y_size;j++){
      t2 = j * pfft->fft_x_size;
      for(k=0;k<pfft->fft_x_size;k++){
	t3 = t1 + t2 + k;

	pfft->fftw_fft[t3][0] = pfft->ws[t3][1]*dd[j];
	pfft->fftw_fft[t3][1] = -1.0*pfft->ws[t3][0]*dd[j];
      }
    }
  }

  /* Do reverse transform. */
  fftw_execute(pfft->fft_backward);

  /* The 2D PSF is the middle plane of the 3D PSF. */
  size_xy = pfft->x_size*pfft->y_size;
  mid_z = size_xy * (pfft->z_size/2);
  
  for(i=0;i<size_xy;i++){
    dy[i] = pfft->fftw_real[mid_z+i];
  }
}

/*
 * pFTGetPSFdz()
 *
 * Return the derivative with respect to z of the current PSF.
 *
 * pfft - A pointer to a psfFFT structure.
 * dz - Pre-allocated storage for the result.
 */
void pFTGetPSFdz(psfFFT *pfft, double *dz)
{
  int i,j,k,t1,t2,t3;
  int mid_z,size_xy;
  double *dd;

  /* Calculate FFT derivative vector. */
  dd = pfft->kz_c;
  pFTCalcShiftVectorDerivative(dd, pfft->z_size);

  /* Copy current PF multiplied by ky into FFTW input. */
  for(i=0;i<pfft->z_size;i++){
    t1 = i * (pfft->y_size * pfft->fft_x_size);
    for(j=0;j<pfft->y_size;j++){
      t2 = j * pfft->fft_x_size;
      for(k=0;k<pfft->fft_x_size;k++){
	t3 = t1 + t2 + k;

	pfft->fftw_fft[t3][0] = pfft->ws[t3][1]*dd[i];
	pfft->fftw_fft[t3][1] = -1.0*pfft->ws[t3][0]*dd[i];
      }
    }
  }
  
  /* Do reverse transform. */
  fftw_execute(pfft->fft_backward);

  /* The 2D PSF is the middle plane of the 3D PSF. */
  size_xy = pfft->x_size*pfft->y_size;
  mid_z = size_xy * (pfft->z_size/2);
  
  for(i=0;i<size_xy;i++){
    dz[i] = -pfft->fftw_real[mid_z+i];
  }
}

/*
 * pFTGetXSize()
 */
int pFTGetXSize(psfFFT *pfft)
{
  return pfft->x_size;  
}

/*
 * pFTGetYSize()
 */
int pFTGetYSize(psfFFT *pfft)
{
  return pfft->y_size;  
}

/*
 * pFTGetZSize()
 */
int pFTGetZSize(psfFFT *pfft)
{
  return pfft->z_size;  
}

/*
 * pFTInitialize()
 *
 * Set things up for FFT based PSF calculations.
 *
 * Note: In the psf array, the z axis is the slowest followed
 *       by the y axis and then the x axis.
 *
 * psf - The psf (z_size, y_size, x_size).
 * z_size - The size of the psf in z (slowest dimension).
 * y_size - The size of the psf in y.
 * x_size - The size of the psf in x (fastest dimension). 
 *
 * return - A psf_fft structure.
 */
psfFFT *pFTInitialize(double *psf, int z_size, int y_size, int x_size)
{
  int i;
  double normalization;
  fftw_plan fft_forward;
  psfFFT *pfft;

  pfft = (psfFFT *)malloc(sizeof(psfFFT));
  
  /* Initialize some variables. */
  pfft->x_size = x_size;
  pfft->y_size = y_size;
  pfft->z_size = z_size;
  pfft->psf_size = z_size * y_size * x_size;

  pfft->fft_x_size = (x_size/2 + 1);
  pfft->fft_size = z_size * y_size * (x_size/2 + 1);

  normalization = 1.0/((double)(pfft->psf_size));

  /* Allocate storage. */
  pfft->kx_c = (double *)malloc(sizeof(double)*x_size);
  pfft->kx_r = (double *)malloc(sizeof(double)*x_size);
  pfft->ky_c = (double *)malloc(sizeof(double)*y_size);
  pfft->ky_r = (double *)malloc(sizeof(double)*y_size);
  pfft->kz_c = (double *)malloc(sizeof(double)*z_size);
  pfft->kz_r = (double *)malloc(sizeof(double)*z_size);

  pfft->fftw_real = (double *)fftw_malloc(sizeof(double)*pfft->psf_size);
  pfft->fftw_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pfft->fft_size);
  pfft->ws = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pfft->fft_size);
  pfft->psf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*pfft->fft_size);

  /* Create FFT plan. */
  pfft->fft_backward = fftw_plan_dft_c2r_3d(z_size, y_size, x_size, pfft->fftw_fft, pfft->fftw_real, FFTW_MEASURE);

  /* Compute FFT of psf and save. */
  fft_forward = fftw_plan_dft_r2c_3d(z_size, y_size, x_size, pfft->fftw_real, pfft->fftw_fft, FFTW_ESTIMATE);
  
  for(i=0;i<pfft->psf_size;i++){
    pfft->fftw_real[i] = psf[i] * normalization;
  }
  
  fftw_execute(fft_forward);
  
  for(i=0;i<pfft->fft_size;i++){
    pfft->psf[i][0] = pfft->fftw_fft[i][0];
    pfft->psf[i][1] = pfft->fftw_fft[i][1];
    pfft->ws[i][0] = pfft->fftw_fft[i][0];
    pfft->ws[i][1] = pfft->fftw_fft[i][1];
  }

  fftw_destroy_plan(fft_forward);

  return pfft;
}

/*
 * pFTTranslate()
 *
 * Translate the psf by dx, dy, dz.  
 */
void pFTTranslate(psfFFT *pfft, double dx, double dy, double dz)
{
  int i,j,k,t1,t2,t3;
  double c1,c2,r1,r2;
  double *sxc,*sxr,*syc,*syr,*szc,*szr;

  /* Calculate FFT translation vectors. */
  sxc = pfft->kx_c;
  sxr = pfft->kx_r;
  pFTCalcShiftVector(sxr, sxc, dx, pfft->x_size);
  
  syc = pfft->ky_c;
  syr = pfft->ky_r;
  pFTCalcShiftVector(syr, syc, dy, pfft->y_size);

  /* We use -dz to match the convention of Z convention of simulator.pupil_math. */
  szc = pfft->kz_c;
  szr = pfft->kz_r;
  pFTCalcShiftVector(szr, szc, -dz, pfft->z_size);

  /* Translate */
  for(i=0;i<pfft->z_size;i++){
    t1 = i * (pfft->y_size * pfft->fft_x_size);
    for(j=0;j<pfft->y_size;j++){
      t2 = j * pfft->fft_x_size;
      for(k=0;k<pfft->fft_x_size;k++){
	t3 = t1 + t2 + k;
	
	r1 = pfft->psf[t3][0];
	c1 = pfft->psf[t3][1];

	/* z shift. */
	r2 = r1 * szr[i] - c1 * szc[i];
	c2 = r1 * szc[i] + c1 * szr[i];

	/* y shift. */
	r1 = r2 * syr[j] - c2 * syc[j];
	c1 = r2 * syc[j] + c2 * syr[j];

	/* x shift. */
	r2 = r1 * sxr[k] - c1 * sxc[k];
	c2 = r1 * sxc[k] + c1 * sxr[k];		

	pfft->ws[t3][0] = r2;
	pfft->ws[t3][1] = c2;
      }
    }
  }  
}
