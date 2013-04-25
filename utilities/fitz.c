/*
 * 07/11
 *
 * Performs z fit based on wx, wy on a Insight3 file.
 * Works "in place".
 *
 *
 * Hazen
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc fitz.c -o fitz -lm
 *
 */


/* Include */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "insight.h"

#define MINZ -500
#define MAXZ 500

double *wx_curve;
double *wy_curve;


/*
 * Initialize wx, wy pre-calculated array curves.
 */

void initWxWy(double *wx_params, double *wy_params)
{
  int i,size;
  double zt,sum;

  size = MAXZ-MINZ;

  // wx
  wx_curve = (double *)malloc(sizeof(double)*size);
  for(i=MINZ;i<MAXZ;i++){
    zt = ((double)i - wx_params[1])/wx_params[2];
    sum = 1.0 + zt*zt + wx_params[3]*zt*zt*zt + wx_params[4]*zt*zt*zt*zt;
    sum += wx_params[5]*zt*zt*zt*zt*zt + wx_params[6]*zt*zt*zt*zt*zt*zt;
    wx_curve[i-MINZ] = sqrt(wx_params[0]*sqrt(sum));
  }

  // wy
  wy_curve = (double *)malloc(sizeof(double)*size);
  for(i=MINZ;i<MAXZ;i++){
    zt = ((double)i - wy_params[1])/wy_params[2];
    sum = 1.0 + zt*zt + wy_params[3]*zt*zt*zt + wy_params[4]*zt*zt*zt*zt;
    sum += wy_params[5]*zt*zt*zt*zt*zt + wy_params[6]*zt*zt*zt*zt*zt*zt;
    wy_curve[i-MINZ] = sqrt(wy_params[0]*sqrt(sum));
  }
}


/*
 * Find the best fitting Z value.
 *
 * This just does a grid search.
 */
float findBestZ(double wx, double wy, double cutoff)
{
  int i,best_i;
  double d,dwx,dwy,best_d;

  best_i = 0;
  dwx = wx - wx_curve[0];
  dwy = wy - wy_curve[0];
  best_d = dwx*dwx+dwy*dwy;
  for(i=1;i<(MAXZ-MINZ);i++){
    dwx = wx - wx_curve[i];
    dwy = wy - wy_curve[i];
    d = dwx*dwx+dwy*dwy;
    if(d<best_d){
      best_i = i;
      best_d = d;
    }
  }

  // distance cut-off here
  if (best_d>cutoff){
    return (float)(MINZ-1.0);
  }
  else {
    return (float)(best_i + MINZ);
  }
}


/*
 * Main
 *
 * i3_file - insight3 file on which to perform z calculations.
 * cut_off - distance cutoff
 * wx_params - 7 numbers (wx0, zc, d, A, B, C, D).
 * wy_params - 7 more numbers (wx0, zc, d, A, B, C, D).
 *
 * Expects calibration curves & molecule widths to be in nm.
 */

int main(int argc, const char *argv[])
{
  int i,bad_cat,molecules,offset;
  float w,a,z;
  double cutoff;
  double wx,wy;
  double wx_params[7];
  double wy_params[7];
  FILE *mlist;

  if (argc == 1){
    printf("usage: best called using the appropriate python script.\n");
    exit(0);
  }

  // setup
  mlist = fopen(argv[1], "rb+");
  cutoff = atof(argv[2]);
  cutoff = cutoff*cutoff;

  bad_cat = 9;

  for(i=3;i<10;i++){
    wx_params[i-3] = atof(argv[i]);
  }

  for(i=10;i<17;i++){
    wy_params[i-10] = atof(argv[i]);
  }

  fseek(mlist, MOLECULES, SEEK_SET);
  fread(&molecules, sizeof(int), 1, mlist);
  printf("Molecules: %d\n", molecules);

  initWxWy(wx_params, wy_params);

  // analysis
  for(i=0;i<molecules;i++){
    if((i%50000)==0){
      printf("Processing molecule %d (fitz)\n", i);
      //printf(" (%f, %f, %f)\n", object_data[X], object_data[Y], object_data[Z]);
    }
    offset = DATA + i*OBJECT_DATA_SIZE*DATUM_SIZE;

    fseeko64(mlist, offset+WIDTH*DATUM_SIZE, SEEK_SET);
    fread(&w, sizeof(float), 1, mlist);

    fseeko64(mlist, offset+ASPECT*DATUM_SIZE, SEEK_SET);
    fread(&a, sizeof(float), 1, mlist);

    wx = sqrt(sqrt(w*w/a));
    wy = sqrt(sqrt(w*w*a));

    z = findBestZ(wx, wy, cutoff);

    // printf("%d %.3f %.3f %.3f\n", i, wx*wx, wy*wy, z);

    /*
    if(z<MINZ){
      fseek(mlist, offset+CAT*DATUM_SIZE, SEEK_SET);
      fwrite(&bad_cat, sizeof(int), 1, mlist);
    }
    */
      
    fseeko64(mlist, offset+ZO*DATUM_SIZE, SEEK_SET);
    fwrite(&z, sizeof(float), 1, mlist);

    fseeko64(mlist, offset+Z*DATUM_SIZE, SEEK_SET);
    fwrite(&z, sizeof(float), 1, mlist);
  }

  // cleanup
  free(wx_curve);
  free(wy_curve);
  fclose(mlist);
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
