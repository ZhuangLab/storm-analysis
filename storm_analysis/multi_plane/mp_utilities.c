/*
 * Utility functions for multi-plane analysis.
 *
 * FIXME: There is some redundancy here with sa_library/ia_utilities.c. At
 *        some point both libraries should be cleaned up and merged?
 *
 * Hazen 06/17
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "../sa_library/multi_fit.h"

typedef struct
{
  int im_size_x;        /* Image size in x (the fast axis). */
  int im_size_y;        /* Image size in y (the slow axis). */

  int n_channels;       /* The number of different channels / image planes. */
  int n_zplanes;        /* The number of different z planes. */

  double radius;        /* Distance for flagging neighbors as too close. */
  double neighborhood;  /* Distance for flagging neighbors as running. */

  double *xt_0toN;      /* Transform x coordinate from channel 0 to channel N. */
  double *yt_0toN;      /* Transform y coordinate from channel 0 to channel N. */  
} mpUtil;

void mpuBadPeakMask(mpUtil *, double *, uint8_t *, int);
void mpuCleanup(mpUtil *);
void mpuFilterPeaks(mpUtil *, double *, double *, uint8_t *, int, int);
mpUtil *mpuInitialize(double, double, int, int, int, int);
void mpuMarkClosePeaks(mpUtil *, double *, uint8_t *, int);
void mpuMergeNewPeaks(mpUtil *, double *, double *, double *, int, int);
void mpuSetTransforms(mpUtil *, double *, double *);
void mpuSplitPeaks(mpUtil *, double *, double *, int);


/*
 * mpuBadPeakMask()
 *
 * Change mask to zero where there was a bad peak. 
 *
 * Note this assumes that STATUS is identical across all channels so
 * only the first channel needs to be checked.
 */
void mpuBadPeakMask(mpUtil *mpu, double *peaks, uint8_t *mask, int n_peaks)
{
  int i;

  for(i=0;i<n_peaks;i++){
    if(peaks[i*NPEAKPAR+STATUS] == ERROR){
      mask[i] = 0;
    }
  }
}

/*
 * mpuCleanup()
 *
 * Free the mpUtil structure.
 */
void mpuCleanup(mpUtil *mpu)
{
  free(mpu->xt_0toN);
  free(mpu->yt_0toN);
  free(mpu);
}

/*
 * mpuFilterPeaks()
 *
 * Copies valid peaks, as specified by mask, into a new
 * (pre-allocated) peak array.
 *
 * Note: peaks with mask != 0 are copied.
 */
void mpuFilterPeaks(mpUtil *mpu, double *in_peaks, double *out_peaks, uint8_t *mask, int n_in, int n_out)
{
  int i,j,k,l,n_ch;

  n_ch = mpu->n_channels;
  j = 0;
  for(i=0;i<n_in;i++){
    if(mask[i]){
      for(k=0;k<n_ch;k++){
	for(l=0;l<NPEAKPAR;l++){
	  out_peaks[(k*n_out+j)*NPEAKPAR+l] = in_peaks[(k*n_in+i)*NPEAKPAR+l];
	}
      }
      j++;
    }
  }
}

/*
 * mpuInitialize()
 *
 * Capture some common parameters that most functions need in a 
 * single structure. This was mostly done to keep the function
 * arguments lists from being really long.
 */
mpUtil *mpuInitialize(double radius, double neighborhood, int im_size_x, int im_size_y, int n_channels, int n_zplanes)
{
  mpUtil *mpu;

  mpu = (mpUtil *)malloc(sizeof(mpUtil));

  mpu->im_size_x = im_size_x;
  mpu->im_size_y = im_size_y;
  mpu->n_channels = n_channels;
  mpu->n_zplanes = n_zplanes;
  mpu->radius = radius;
  mpu->neighborhood = neighborhood;

  mpu->xt_0toN = (double *)malloc(sizeof(double)*n_channels*3);
  mpu->yt_0toN = (double *)malloc(sizeof(double)*n_channels*3);

  return mpu;
}

/*
 * mpuMarkClosePeaks()
 *
 * Identify peaks that are too close to each other. Mark the dimmer
 * one for deletion and it's neighbors as RUNNING.
 *
 * Note: mask should be initialized to an array of 1s.
 *
 * FIXME: Use kdtree? The current approach scales poorly with the number of peaks.
 */
void mpuMarkClosePeaks(mpUtil *mpu, double *peaks, uint8_t *mask, int n_peaks)
{
  int bad,i,j,k;
  double dx,dy,h,neighborhood,radius,x,y;

  radius = mpu->radius * mpu->radius;
  neighborhood = mpu->neighborhood * mpu->neighborhood;

  /* Mark peaks that need to be removed. */
  for(i=0;i<n_peaks;i++){
    x = peaks[i*NPEAKPAR+XCENTER];
    y = peaks[i*NPEAKPAR+YCENTER];
    h = peaks[i*NPEAKPAR+HEIGHT];
    bad = 0;
    j = 0;
    while((j<n_peaks)&&(!bad)){
      if(j!=i){
	dx = x - peaks[j*NPEAKPAR+XCENTER];
	dy = y - peaks[j*NPEAKPAR+YCENTER];
	if(((dx*dx+dy*dy)<radius)&&(peaks[j*NPEAKPAR+HEIGHT]>h)){
	  bad = 1;
	}
      }
      j++;
    }
    if(bad){
      mask[i] = 0;
    }
  }

  /* Flag (non-bad) neighbors of bad peaks as running. */
  for(i=0;i<n_peaks;i++){
    if(mask[i]==0){
      x = peaks[i*NPEAKPAR+XCENTER];
      y = peaks[i*NPEAKPAR+YCENTER];
      for(j=0;j<n_peaks;j++){
	if(j!=i){
	  dx = x - peaks[j*NPEAKPAR+XCENTER];
	  dy = y - peaks[j*NPEAKPAR+YCENTER];
	  if(((dx*dx+dy*dy)<neighborhood)&&(mask[j])){
	    for(k=0;k<(mpu->n_channels);k++){
	      peaks[(k*n_peaks+j)*NPEAKPAR+STATUS] = RUNNING;
	    }
	  }
	}
      }
    }
  }
}

/*
 * mpuMergeNewPeaks()
 *
 * Merge the current peaks list with a new peak list. The number of peaks
 * arguments are the number of unique peaks, not the total number of peaks.
 *
 * After merging the list will need to be filtered to remove ERROR peaks.
 *
 * FIXME: If the peaks are far enough away in z then overlap should be allowed.
 *
 * FIXME: Use kdtree? The current approach scales poorly with the number of peaks.
 */
void mpuMergeNewPeaks(mpUtil *mpu, double *cur_peaks, double *new_peaks, double *merged_peaks, int n_cur, int n_new)
{
  int bad,i,j,k,l,m,n_ch,n_tot;
  double dd,dx,dy,radius,neighborhood,x,y;

  n_ch = mpu->n_channels;
  n_tot = n_cur + n_new;
  
  /* First copy cur_peaks into merged_peaks. */
  for(i=0;i<n_cur;i++){
    for(j=0;j<n_ch;j++){
      k = (j*n_tot+i)*NPEAKPAR;
      l = (j*n_cur+i)*NPEAKPAR;
      for(m=0;m<NPEAKPAR;m++){
	merged_peaks[k+m] = cur_peaks[l+m];
      }
    }
  }

  /* Then copy new_peaks into merged_peaks. */
  for(i=0;i<n_new;i++){
    for(j=0;j<n_ch;j++){
      k = (j*n_tot+i+n_cur)*NPEAKPAR;
      l = (j*n_new+i)*NPEAKPAR;
      for(m=0;m<NPEAKPAR;m++){
	merged_peaks[k+m] = new_peaks[l+m];
      }
    }
  }
  
  /* 
   * Then check new peaks and add if they are not too close to an existing 
   * peak. We do the check on channel 0, but mark the peaks in all the
   * channels.
   */
  radius = mpu->radius * mpu->radius;
  neighborhood = mpu->neighborhood * mpu->neighborhood;
  for(i=0;i<n_new;i++){
    x = new_peaks[i*NPEAKPAR+XCENTER];
    y = new_peaks[i*NPEAKPAR+YCENTER];
    bad = 0;
    j = 0;
    while((j<n_new)&&(!bad)){
      dx = x - cur_peaks[j*NPEAKPAR+XCENTER];
      dy = y - cur_peaks[j*NPEAKPAR+XCENTER];
      dd = dx*dx + dy*dy;
      if(dd < radius){
	bad = 1;
      }
      else if (dd < neighborhood){
	for(k=0;k<n_ch;k++){
	  merged_peaks[(k*n_tot+j)*NPEAKPAR+STATUS] = RUNNING;
	}
      }
      j++;
    }
    if(bad){
      for(k=0;k<n_ch;k++){
	merged_peaks[(k*n_tot+i)*NPEAKPAR+STATUS] = ERROR;
      }
    }
  }
}

/*
 * mpuSetTransforms()
 *
 * Set the channel0 to channelN transform arrays.
 */
void mpuSetTransforms(mpUtil *mpu, double *xt_0toN, double *yt_0toN)
{
  int i;

  for(i=0;i<(mpu->n_channels*3);i++){
    mpu->xt_0toN[i] = xt_0toN[i];
    mpu->yt_0toN[i] = yt_0toN[i];
  }
}
		      
/*
 * mpuSplitPeaks()
 *
 * Create peaks for all channels from channel 0 peaks.
 */
void mpuSplitPeaks(mpUtil *mpu, double *peaks, double *split_peaks, int n_peaks)
{
  int i,j,k;
  double xi,yi,t;

  /* First just duplicate everything. */
  for(i=0;i<(mpu->n_channels);i++){
    for(j=0;j<n_peaks;j++){
      for(k=0;k<NPEAKPAR;k++){
	split_peaks[(i*n_peaks+j)*NPEAKPAR+k] = peaks[j*NPEAKPAR+k];
      }
    }
  }

  /* Then fix x,y coordinates. */
  for(i=0;i<n_peaks;i++){
    xi = peaks[i*NPEAKPAR+XCENTER];
    yi = peaks[i*NPEAKPAR+YCENTER];
    for(j=1;j<(mpu->n_channels);j++){

      /* x transform. */
      t = mpu->xt_0toN[j*3] + mpu->xt_0toN[j*3+1]*xi + mpu->xt_0toN[j*3+2]*yi;
      split_peaks[(i+j*n_peaks)*NPEAKPAR+XCENTER] = t;

      /* y transform. */
      t = mpu->yt_0toN[j*3] + mpu->yt_0toN[j*3+1]*xi + mpu->yt_0toN[j*3+2]*yi;
      split_peaks[(i+j*n_peaks)*NPEAKPAR+YCENTER] = t;
    }
  }
}
