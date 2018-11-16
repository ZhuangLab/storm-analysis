/*
 * Finds the best z value based on the localizations width wx, wy and
 * calibration curves.
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

typedef struct
{
  int size;

  double cut_off;
  
  double *z_values;
  double *wx_curve;
  double *wy_curve;
} zfitData;

void cleanup(zfitData *);
zfitData *initialize(double *, double *, double, double, double, double);
double findBestZ(zfitData *, double, double);
double findMinimumDistance(zfitData *, double, double);


/*
 * cleanup()
 */
void cleanup(zfitData *zfit_data){
  free(zfit_data->z_values);
  free(zfit_data->wx_curve);
  free(zfit_data->wy_curve);
  free(zfit_data);
}

/*
 * initialize()
 *
 * Initialize a zfitData structure using z calibration data.
 */
zfitData *initialize(double *wx_params, double *wy_params, double z_min, double z_max, double z_step, double cut_off)
{
  int i,size;
  double z,zt,sum;
  zfitData *zfit_data;

  size = (int)((z_max - z_min)/z_step) + 1;
  
  zfit_data = (zfitData *)malloc(sizeof(zfitData));
  zfit_data->size = size;
  zfit_data->cut_off = cut_off;
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
 * findBestZ()
 *
 * Find the best fitting Z value, this just does a grid search.
 */
double findBestZ(zfitData *zfit_data, double wx, double wy)
{
  int i;
  double d,dwx,dwy,best_d,best_z,rval;

  wx = sqrt(wx);
  wy = sqrt(wy);

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
  if (best_d > zfit_data->cut_off){
    rval = zfit_data->z_values[0] - 1.0;
  }
  else {
    rval = best_z;
  }

  return rval;
}

/*
 * findMinimumDistance()
 *
 * Find the minimum distance to calibration curve 
 * for a given wx, wy.
 */
double findMinimumDistance(zfitData *zfit_data, double wx, double wy)
{
  int i;
  double d,dwx,dwy,best_d;

  wx = sqrt(wx);
  wy = sqrt(wy);

  dwx = wx - zfit_data->wx_curve[0];
  dwy = wy - zfit_data->wy_curve[0];
  best_d = dwx*dwx+dwy*dwy;
  
  for(i=1;i<zfit_data->size;i++){
    dwx = wx - zfit_data->wx_curve[i];
    dwy = wy - zfit_data->wy_curve[i];
    d = dwx*dwx+dwy*dwy;
    if(d<best_d){
      best_d = d;
    }
  }

  return best_d;
}

/*
 * The MIT License
 *
 * Copyright (c) 2018 Zhuang Lab, Harvard University
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
