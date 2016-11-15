/*
 * Performs z fit based on wx, wy on a Insight3 file.
 *
 * Notes:
 *   (1) Works "in place" on the file.
 */


/* Include */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "insight.h"

typedef struct
{
  int size;

  double *z_values;
  double *wx_curve;
  double *wy_curve;
} zfitData;

zfitData *initialize(double *, double *, double, double, double);
float findBestZ(zfitData *, double, double, double);
int fitz(const char *, double *, double *, double, double, double, double);


/*
 * Initialize wx, wy pre-calculated array curves.
 */
zfitData *initWxWy(double *wx_params, double *wy_params, double z_min, double z_max, double z_step)
{
  int i,size;
  double z,zt,sum;
  zfitData *zfit_data;

  size = (int)((z_max - z_min)/z_step) + 1;
  
  zfit_data = (zfitData *)malloc(sizeof(zfitData));
  zfit_data->size = size;
  zfit_data->z_values = (double *)malloc(sizeof(double)*size);
  zfit_data->wx_curve = (double *)malloc(sizeof(double)*size);
  zfit_data->wy_curve = (double *)malloc(sizeof(double)*size);

  // wx
  i = 0;
  z = z_min;
  for(i=0;i<size;i++){
    zfit_data->z_values[i] = z;

    /* wx */
    zt = (z - wx_params[1])/wx_params[2];
    sum = 1.0 + zt*zt + wx_params[3]*zt*zt*zt + wx_params[4]*zt*zt*zt*zt;
    sum += wx_params[5]*zt*zt*zt*zt*zt + wx_params[6]*zt*zt*zt*zt*zt*zt;
    zfit_data->wx_curve[i] = sqrt(wx_params[0]*sqrt(sum));

    /* wy */
    zt = ((double)z - wy_params[1])/wy_params[2];
    sum = 1.0 + zt*zt + wy_params[3]*zt*zt*zt + wy_params[4]*zt*zt*zt*zt;
    sum += wy_params[5]*zt*zt*zt*zt*zt + wy_params[6]*zt*zt*zt*zt*zt*zt;
    zfit_data->wy_curve[i] = sqrt(wy_params[0]*sqrt(sum));
    
    z += z_step;
  }

  return zfit_data;
}

/*
 * Find the best fitting Z value.
 *
 * This just does a grid search.
 */
float findBestZ(zfitData *zfit_data, double wx, double wy, double cutoff)
{
  int i;
  double d,dwx,dwy,best_d,best_z;

  dwx = wx - zfit_data->wx_curve[0];
  dwy = wy - zfit_data->wy_curve[0];
  best_z = zfit_data->z_values[0];
  best_d = dwx*dwx+dwy*dwy;
  
  for(i=1;i<zfit_data->size;i++){
    dwx = wx - zfit_data->wx_curve[i];
    dwy = wy - zfit_data->wy_curve[i];
    d = dwx*dwx+dwy*dwy;
    if(d<best_d){
      best_z = zfit_data->z_values[i];
      best_d = d;
    }
  }

  // distance cut-off here
  if (best_d>cutoff){
    return (float)(zfit_data->z_values[0] - 1.0);
  }
  else {
    return (float)(best_z);
  }
}


/*
 * fitz
 *
 * i3_filename - Insight3 file on which to perform z calculations.
 * wx_params - 7 numbers (wx0, zc, d, A, B, C, D).
 * wy_params - 7 more numbers (wx0, zc, d, A, B, C, D).
 * cut_off - distance cutoff.
 * z_min - minimum z value.
 * z_max - maximum z value.
 * z_step - z step size.
 *
 * Expects calibration curves & molecule widths to be in nm.
 */

int fitz(const char *i3_filename, double *wx_params, double *wy_params, double cutoff, double z_min, double z_max, double z_step)
{
  int i,bad_cat,molecules,offset;
  float w,a,z;
  double minz,wx,wy;
  size_t n_read;
  zfitData *zfit_data;
  FILE *mlist;

  /* Setup */
  mlist = fopen(i3_filename, "rb+");
  cutoff = cutoff*cutoff;

  bad_cat = 9;

  fseek(mlist, MOLECULES, SEEK_SET);
  n_read = fread(&molecules, sizeof(int), 1, mlist);
  if(n_read != 1) return 1;  
  printf("Molecules: %d\n", molecules);

  zfit_data = initWxWy(wx_params, wy_params, z_min, z_max, z_step);
  minz = zfit_data->z_values[0] - 0.1;

  /* Analyze the file. */
  for(i=0;i<molecules;i++){
    if((i%50000)==0){
      printf("Processing molecule %d (fitz)\n", i);
      //printf(" (%f, %f, %f)\n", object_data[X], object_data[Y], object_data[Z]);
    }
    offset = DATA + i*OBJECT_DATA_SIZE*DATUM_SIZE;

    fseeko64(mlist, offset+WIDTH*DATUM_SIZE, SEEK_SET);
    n_read = fread(&w, sizeof(float), 1, mlist);
    if(n_read != 1) return 1;
  
    fseeko64(mlist, offset+ASPECT*DATUM_SIZE, SEEK_SET);
    n_read = fread(&a, sizeof(float), 1, mlist);
    if(n_read != 1) return 1;
    
    wx = sqrt(sqrt(w*w/a));
    wy = sqrt(sqrt(w*w*a));

    z = findBestZ(zfit_data, wx, wy, cutoff);
    // printf("%d %.3f %.3f %.3f\n", i, wx*wx, wy*wy, z);

    if(z<minz){
      fseek(mlist, offset+CAT*DATUM_SIZE, SEEK_SET);
      fwrite(&bad_cat, sizeof(int), 1, mlist);
    }
      
    fseeko64(mlist, offset+ZO*DATUM_SIZE, SEEK_SET);
    fwrite(&z, sizeof(float), 1, mlist);

    fseeko64(mlist, offset+Z*DATUM_SIZE, SEEK_SET);
    fwrite(&z, sizeof(float), 1, mlist);
  }

  /* Cleanup. */
  fclose(mlist);

  free(zfit_data->z_values);
  free(zfit_data->wx_curve);
  free(zfit_data->wy_curve);
  free(zfit_data);

  return 0;
}


/*
 * The MIT License
 *
 * Copyright (c) 2012 Zhuang Lab, Harvard University
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
