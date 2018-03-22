/*
 * C library for doing some of heavy lifting of pupil function 
 * calculations. This includes calculating zernike modes and
 * all the Fourier transforms for displacing the pupil function
 * in z, calculating the PSF and OTF rescaling.
 *
 * Hazen 02/18
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <fftw3.h>

typedef struct pfData
{
  int size;     /* The size of pupil function in X/Y in pixels. */
  int sf_set;   /* Scaling factor was set. */

  double *kz;   /* Z translation multiplier. */
  double *sf;   /* OTF scaling factor. */

  fftw_complex *pf;  
  fftw_complex *fftw_pf;  /* Used by the FFT. */
  fftw_complex *fftw_psf; /* Used by the FFT. */

  fftw_plan pf_to_psf;
  fftw_plan psf_to_pf;
} pfData;

void pfCleanup(pfData *);
void pfGetPSF(pfData *, double *, double);
int pfFactorial(int);
pfData *pfInitialize(double *, int);
double pfPreFac(int, int, int);
void pfResetScaling(pfData *);
void pfSetPF(pfData *, double *, double *);
void pfSetScaling(pfData *, double *);
double pfZernike(int, int, double, double);
void pfZernikeGrid(double *, int, double, double, double, int, int);
double pfZernikeRad(int, int, double);


/*
 * pfCleanup()
 *
 * Free the pfData structure.
 */
void pfCleanup(pfData *pf_data)
{
  free(pf_data->kz);
  free(pf_data->sf);

  fftw_free(pf_data->pf);
  fftw_free(pf_data->fftw_pf);
  fftw_free(pf_data->fftw_psf);

  fftw_destroy_plan(pf_data->pf_to_psf);
  fftw_destroy_plan(pf_data->psf_to_pf);

  free(pf_data);
}

/*
 * pfGetPSF()
 *
 * Returns the OTF scaled PSF at the requested z value.
 *
 * psf - Pre-allocated storage for the PSF.
 * dz - The z offset in microns.
 */
void pfGetPSF(pfData *pf_data, double *psf, double dz)
{
  int i;
  double dd,r1,c1;

  /* Translate PF in Z if necessary. */
  if (dz != 0.0){
    for(i=0;i<(pf_data->size*pf_data->size);i++){
      dd = -pf_data->kz[i]*dz;
      r1 = cos(dd);
      c1 = -sin(dd);
      pf_data->fftw_pf[i][0] = r1*pf_data->pf[i][0] - c1*pf_data->pf[i][1];
      pf_data->fftw_pf[i][1] = r1*pf_data->pf[i][1] + c1*pf_data->pf[i][0];
    }
  }
  else{
    for(i=0;i<(pf_data->size*pf_data->size);i++){
      pf_data->fftw_pf[i][0] = pf_data->pf[i][0];
      pf_data->fftw_pf[i][1] = pf_data->pf[i][1];
    }
  }

  /* FFT to real space. */
  fftw_execute(pf_data->pf_to_psf);

  /* Apply OTF scaling, if requested. */
  if (pf_data->sf_set){
    for(i=0;i<(pf_data->size*pf_data->size);i++){
      r1 = pf_data->fftw_psf[i][0];
      c1 = pf_data->fftw_psf[i][1];
      pf_data->fftw_psf[i][0] = r1*r1 + c1*c1;
      pf_data->fftw_psf[i][1] = 0.0;
    }

    /* Back to Fourier space. */
    fftw_execute(pf_data->psf_to_pf);

    /* Apply scaling. */
    for(i=0;i<(pf_data->size*pf_data->size);i++){
      pf_data->fftw_pf[i][0] = pf_data->fftw_pf[i][0] * pf_data->sf[i];
      pf_data->fftw_pf[i][1] = pf_data->fftw_pf[i][1] * pf_data->sf[i];
    }

    /* And back again. */
    fftw_execute(pf_data->pf_to_psf);

    /* Calculate PSF. */
    for(i=0;i<(pf_data->size*pf_data->size);i++){
      r1 = pf_data->fftw_psf[i][0];
      c1 = pf_data->fftw_psf[i][1];
      psf[i] = sqrt(r1*r1 + c1*c1);
    }
  }
  else{
    
    /* Calculate PSF. */
    for(i=0;i<(pf_data->size*pf_data->size);i++){
      r1 = pf_data->fftw_psf[i][0];
      c1 = pf_data->fftw_psf[i][1];
      psf[i] = r1*r1 + c1*c1;
    }
  }
}

/*
 * pfFactorial()
 *
 * Calculate factorial of a integer.
 *
 * n - The input number.
 *
 * Returns n!
 */
int pfFactorial(int n)
{
  int i, n_fac;

  n_fac = 1;
  for(i=1;i<=n;i++){
    n_fac = n_fac * i;
  }

  return n_fac;
}

/*
 * pfInitialize()
 *
 * Create a pfData structure.
 *
 * kz - Vector to use for z offsets.
 * size - The size of the pupil function.
 */
pfData *pfInitialize(double *kz, int size)
{
  int i;
  pfData *pf_data;

  pf_data = (pfData *)malloc(sizeof(pfData));
  pf_data->size = size;
  pf_data->sf_set = 0;

  pf_data->kz = (double *)malloc(sizeof(double)*size*size);
  pf_data->sf = (double *)malloc(sizeof(double)*size*size);

  pf_data->pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*size*size);
  pf_data->fftw_pf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*size*size);
  pf_data->fftw_psf = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*size*size);

  pf_data->pf_to_psf = fftw_plan_dft_2d(size, size, pf_data->fftw_pf, pf_data->fftw_psf, FFTW_BACKWARD, FFTW_MEASURE);
  pf_data->psf_to_pf = fftw_plan_dft_2d(size, size, pf_data->fftw_psf, pf_data->fftw_pf, FFTW_FORWARD, FFTW_MEASURE);

  for(i=0;i<(size*size);i++){
    pf_data->kz[i] = 2.0*M_PI*kz[i];
  }

  return pf_data;
}

/*
 * pfPreFac()
 *
 * Calculate factorial coefficient.
 *
 * m - zernike m value.
 * n - zernike n value.
 * k - index.
 *
 * Returns the factorial coefficient.
 */
double pfPreFac(int m, int n, int k)
{
  double d1, d2, d3, n1, sign;

  sign = pow(-1.0, k);
  n1 = (double)pfFactorial(n-k);
  d1 = (double)pfFactorial(k);
  d2 = (double)pfFactorial((n+m)/2 - k);
  d3 = (double)pfFactorial((n-m)/2 - k);

  return sign * n1/(d1 * d2 * d3);
}

/*
 * pfResetScaling()
 *
 * Turns off OTF scaling.
 */
void pfResetScaling(pfData *pf_data)
{
  pf_data->sf_set = 0;
}

/*
 * pfSetPF()
 *
 * Set the pupil function.
 *
 * r_pf - The real part of the pupil function.
 * c_pf - The complex part of the pupil function.
 */
void pfSetPF(pfData *pf_data, double *r_pf, double *c_pf)
{
  int i,j,k,l;
  double norm;

  norm = 1.0/((double)pf_data->size);

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
  for(i=0;i<pf_data->size;i++){
    j = i * pf_data->size;
    for(k=0;k<pf_data->size;k++){
      l = j+k;
      if(((i+k)%2)==0){
	pf_data->pf[l][0] = r_pf[l] * norm;
  	pf_data->pf[l][1] = c_pf[l] * norm;
      }
      else{
	pf_data->pf[l][0] = -1.0*r_pf[l] * norm;
	pf_data->pf[l][1] = -1.0*c_pf[l] * norm;
      }
    }
  }
}

/* 
 * pfSetScaling()
 *
 * Set OTF scaling term.
 *
 * scaling_factor - The OTF scaling factor.
 */
void pfSetScaling(pfData *pf_data, double *scaling_factor)
{
  int i;
  double norm;

  /* Normalize the OTF scaling now so that we only have to do this once. */
  norm = 1.0/(((double)pf_data->size)*((double)pf_data->size));
  
  pf_data->sf_set = 1;

  for(i=0;i<(pf_data->size*pf_data->size);i++){
    pf_data->sf[i] = scaling_factor[i]*norm;
  }
}

/*
 * pfZernike()
 *
 * Calculate the value of a zernike polynomial.
 *
 * m - zernike m value.
 * n - zernike n value.
 * rho - radius (0.0 - 1.0).
 * phi - angle (in radians).
 *
 * Returns the zernike polynomial value.
 */
double pfZernike(int m, int n, double rho, double phi)
{
  if (m > 0) return pfZernikeRad(m, n, rho) * cos(m * phi);
  if (m < 0) return pfZernikeRad(-m, n, rho) * sin(-m * phi);
  return pfZernikeRad(0, n, rho);
}

/*
 * pfZernikeGrid()
 *
 * Add zernike polynomial values to pre-defined grid.
 *
 * grid - pre-allocated & initialized double storage (square).
 * gsize - size of the grid in x / y.
 * gcenter - center point of the grid.
 * radius - radius on which to calculate the polynomial.
 * scale - scaling factor.
 * m - zernike m value.
 * n - zernike n value.
 */
void pfZernikeGrid(double *grid, int gsize, double gcenter, double radius, double scale, int m, int n)
{
  int i, j;
  double dd, dx, dy, phi, rr;

  rr = radius * radius;
  for(i=0;i<gsize;i++){
    dx = i - gcenter;
    for(j=0;j<gsize;j++){
      dy = j - gcenter;
      dd = dx * dx + dy * dy;
      if(dd <= rr){
	dd = sqrt(dd)/radius;
	phi = atan2(dy, dx);
	grid[i*gsize+j] += scale * pfZernike(m, n, dd, phi);
      }
    }
  }
}

/*
 * pfZernikeRad()
 *
 * Calculate the radial component of a Zernike polynomial.
 *
 * m - m coefficient.
 * n - n coefficient.
 * rho - radius (0.0 - 1.0).
 *
 * Returns the radial value.
 */
double pfZernikeRad(int m, int n, double rho)
{
  int k;
  double sum;

  if((n < 0) || (m < 0) || (abs(m) > n)) return 0.0;
  if(((n-m)%2) == 1) return 0.0;

  sum = 0.0;
  for(k=0;k<((n-m)/2+1);k++){
    sum += pfPreFac(m, n, k) * pow(rho, (n - 2*k));
  }
  
  return sum;
}
