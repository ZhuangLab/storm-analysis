/*
 * C library for OTF scaling of a pupil function PSF. These are
 * are always square.
 *
 * Hazen 3/16 
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fftw3.h>

/* Structures & Types */
typedef struct otfData {
  int fft_size;
  int size;

  double normalization;

  double *fft_vector;
  double *otf_vector;
  
  fftw_plan fft_backward;
  fftw_plan fft_forward;

  fftw_complex *fft_vector_fft;
} otfData;

/* Function Declarations */
void cleanup(otfData *);
otfData *initialize(int, int);
void scale(otfData *, double *, double *);
void setScale(otfData *, double *);


/*
 * cleanup()
 *
 * otf_data - A pointer to a otfData structure.
 */
void cleanup(otfData *otf_data)
{
  fftw_free(otf_data->fft_vector);
  fftw_free(otf_data->otf_vector);
  
  fftw_destroy_plan(otf_data->fft_backward);
  fftw_destroy_plan(otf_data->fft_forward);

  fftw_free(otf_data->fft_vector_fft);

  free(otf_data);
}


/*
 * initialize()
 *
 * Set things up for OTF scaling.
 *
 * size - the size of the OTF in x and y.
 * estimate - 0/1 to just use an estimated FFT plan. If you are only going to
 *            to use the FFT a few times this can be much faster.
 */
otfData *initialize(int size, int estimate)
{
  otfData *otf_data;

  otf_data = (otfData *)malloc(sizeof(otfData));
  
  /* Initialize some variables. */
  otf_data->fft_size = size * (size/2 + 1);
  otf_data->size = size;
  otf_data->normalization = 1.0/((double)(size * size));

  /* Allocate storage. */
  otf_data->fft_vector = (double *)fftw_malloc(sizeof(double)*size*size);
  otf_data->otf_vector = (double *)fftw_malloc(sizeof(double)*otf_data->fft_size);
  
  otf_data->fft_vector_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*size*size);

  /* Create FFT plans. */
  if (estimate){
    otf_data->fft_forward = fftw_plan_dft_r2c_2d(size, size, otf_data->fft_vector, otf_data->fft_vector_fft, FFTW_ESTIMATE);
    otf_data->fft_backward = fftw_plan_dft_c2r_2d(size, size, otf_data->fft_vector_fft, otf_data->fft_vector, FFTW_ESTIMATE);
  }
  else {
    otf_data->fft_forward = fftw_plan_dft_r2c_2d(size, size, otf_data->fft_vector, otf_data->fft_vector_fft, FFTW_MEASURE);
    otf_data->fft_backward = fftw_plan_dft_c2r_2d(size, size, otf_data->fft_vector_fft, otf_data->fft_vector, FFTW_MEASURE);    
  }

  return otf_data;
}


/*
 * scale()
 *
 * Apply OTF scaling to a psf.
 *
 * otf_data - A pointer to a otfData structure.
 * psf - The psf, must be the same size as the OTF.
 * result - Pre-allocated storage for the result of OTF scaling.
 */
void scale(otfData *otf_data, double *psf, double *result)
{
  int i,size;

  size = otf_data->size;
  
  /* Compute FFT of the PSF. */
  for(i=0;i<(size*size);i++){
    otf_data->fft_vector[i] = psf[i];
  }
  fftw_execute(otf_data->fft_forward);

  /* Multiple by the OTF and compute inverse FFT. */
  for(i=0;i<otf_data->fft_size;i++){
    otf_data->fft_vector_fft[i][0] = otf_data->fft_vector_fft[i][0] * otf_data->otf_vector[i];
    otf_data->fft_vector_fft[i][1] = otf_data->fft_vector_fft[i][1] * otf_data->otf_vector[i];
  }
  fftw_execute(otf_data->fft_backward);

  /* Normalize and copy into result. */
  for(i=0;i<(size*size);i++){
    result[i] = otf_data->fft_vector[i] * otf_data->normalization;
  }  
}


/*
 * setScale()
 *
 * Set array to use for OTF scaling.
 *
 * otf_data - A pointer to a otfData structure.
 * otf_scale - A vector containing the OTF scaling values.
 */
void setScale(otfData *otf_data, double *otf_scale)
{
  int i,j,k,size;

  size = otf_data->size;
  
  /* 
   * Copy OTF scaling array with some fiddling to get things in the correct order 
   * for post FFT multiplication.
   */
  i = 0;
  for(j=0;j<size;j++){
    for(k=0;k<((size/2)+1);k++){
      otf_data->otf_vector[i] = otf_scale[j*size+k];
      i++;
    }
  }
}
