/*
 * Applies drift correction to a molecule list,
 * works in place.
 *
 * Hazen
 * 12/11
 *
 * Compilation instructions:
 *
 * Linux:
 *  gcc apply-drift-correction.c -o apply-drift-correction
 *
 * Windows:
 *  gcc apply-drift-correction.c -o apply-drift-correction
 */


/* Include */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "insight.h"


/*
 * Main
 *
 * mlist - the molecule list file
 * drift - the drift correction file (in i3 standard format).
 *
 */

int main(int argc, const char *argv[])
{
  int i,cur_frame,frames,molecules,temp;
  int *object_data_int;
  char str[100];
  float *dx,*dy,*dz;
  float object_data[OBJECT_DATA_SIZE];
  FILE *mlist_fp,*drift_fp;

  if (argc != 3){
    printf("usage apply-drift-correction <mlist.bin file> <drift.txt file>\n");
    exit(0);
  }

  printf("Applying drift correction\n");

  /* 
   * Setup 
   */
  object_data_int = (int *)object_data;

  // Figure out how many molecules there are to process.
  mlist_fp = fopen(argv[1], "rb+");
  fseek(mlist_fp, MOLECULES, SEEK_SET);
  fread(&molecules, sizeof(int), 1, mlist_fp);
  printf(" Molecules: %d\n", molecules);

  // Determine size of drift correction file & load into memory.
  frames = 0;
  drift_fp = fopen(argv[2], "r");
  while(fgets(str,sizeof(str),drift_fp)!=NULL){
    frames++;
  }
  printf(" Frames: %d\n", frames);
  fseek(drift_fp, 0, SEEK_SET);

  dx = (float *)malloc(sizeof(float)*frames);
  dy = (float *)malloc(sizeof(float)*frames);
  dz = (float *)malloc(sizeof(float)*frames);
  for(i=0;i<frames;i++){
    fgets(str,sizeof(str),drift_fp);
    sscanf(str,"%d\t%f\t%f\t%f",&temp,&dx[i],&dy[i],&dz[i]);
  }
  fclose(drift_fp);

  // Print the last values as verification that 
  // the file was ready properly.
  for(i=(frames-6);i<frames;i++){
    printf("  %d %.3f %.3f %.3f\n",i,dx[i],dy[i],dz[i]);
  }

  /*
   * Go through all the molecules & apply drift correction.
   */
  cur_frame = 0;
  for(i=0;i<molecules;i++){
    //for(i=29800000;i<30000000;i++){
    if((i%50000)==0){
      printf(" Processing molecule %d in frame %d (apply-drift-correction)\n", i, cur_frame);
    }
    fseeko64(mlist_fp, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)i, SEEK_SET);
    fread(&object_data, sizeof(float), OBJECT_DATA_SIZE, mlist_fp);
    cur_frame = object_data_int[FRAME]-1;

    // range checking
    if(cur_frame < 0){
      cur_frame = 0;
    }
    if(cur_frame >= frames){
      cur_frame = frames;
    }
    
    // apply correction
    object_data[X] = object_data[XO]-dx[cur_frame];
    object_data[Y] = object_data[YO]-dy[cur_frame];
    object_data[Z] = object_data[ZO]-dz[cur_frame];

    fseeko64(mlist_fp, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)i, SEEK_SET);
    fwrite(object_data, sizeof(float), OBJECT_DATA_SIZE, mlist_fp);
  }
  printf("Finished applying drift correction\n");

  free(dx);
  free(dy);
  free(dz);
  fclose(mlist_fp);
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
