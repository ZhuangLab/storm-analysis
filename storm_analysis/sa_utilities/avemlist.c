/*
 * 07/11
 *
 * Compute average object molecule list using
 * a tracked molecule list.
 *
 *
 * Hazen
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc avemlist.c -o avemlist -lm
 *
 */


/* Include */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "insight.h"

/* Define */
#define NOAVERAGE 0
#define AVERAGE 1
#define TOTAL 2
#define TESTING 0

/* Functions */
void averageTrack(FILE *, FILE *, int, int);


/* These are as in the insight.h file */
static int average_flag[] = {AVERAGE,      /* XO */
			     AVERAGE,      /* YO */
			     AVERAGE,      /* X */
			     AVERAGE,      /* Y */
			     TOTAL,        /* HEIGHT */
			     TOTAL,        /* AREA */
			     AVERAGE,      /* WIDTH */
			     NOAVERAGE,    /* VISITED */
			     AVERAGE,      /* ASPECT */
			     AVERAGE,      /* BACKGROUND */
			     NOAVERAGE,    /* SUM */
			     NOAVERAGE,    /* CAT */
			     NOAVERAGE,    /* FITI */
			     NOAVERAGE,    /* FRAME */
			     NOAVERAGE,    /* TLEN */
			     NOAVERAGE,    /* LINK */
			     AVERAGE,      /* ZO */
			     AVERAGE};     /* Z */


/*
 * Follows links between molecules to generate the average track.
 * Values are weighted by the square root of the object fit height.
 */

void averageTrack(FILE *input_mlist, FILE *output_mlist, int molecule, int visited)
{
  int i,*object_data_int,elements,track_id;
  float weight, total_weight;
  float average_data[OBJECT_DATA_SIZE], object_data[OBJECT_DATA_SIZE];
  
  elements = 1;
  object_data_int = (int *)object_data;

  // load object data
  fseeko64(input_mlist, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)molecule, SEEK_SET);
  fread(&object_data, sizeof(float), OBJECT_DATA_SIZE, input_mlist);
  for(i=0;i<(OBJECT_DATA_SIZE);i++){
    average_data[i] = object_data[i];
  }
  track_id = object_data_int[FITI];

  // normalize relevant fields
  weight = sqrt(object_data[HEIGHT]);
  for(i=0;i<(OBJECT_DATA_SIZE);i++){
    if(average_flag[i] == AVERAGE){
      average_data[i] = average_data[i]*weight;
    }
  }
  total_weight = weight;
  
  // mark as visited
  object_data_int[VISITED] = visited;
  fseeko64(input_mlist, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)molecule, SEEK_SET);
  fwrite(&object_data, sizeof(float), OBJECT_DATA_SIZE, input_mlist);  

  // printf("\n");
  while(object_data_int[LINK]>0){
    // load object data
    molecule = object_data_int[LINK];
    // printf(" %d\n", molecule);
    fseeko64(input_mlist, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)molecule, SEEK_SET);
    fread(&object_data, sizeof(float), OBJECT_DATA_SIZE, input_mlist);  
    
    if (TESTING){
      if(track_id != object_data_int[FITI]){
	printf("Tracking error detected. %d %d\n", track_id, object_data_int[FITI]);
	printf("  %.3f %.3f %d %d\n", object_data[XO], object_data[YO], object_data_int[CAT], object_data_int[FRAME]);
      }
    }

    // average/sum relevant fields
    weight = sqrt(object_data[HEIGHT]);
    for(i=0;i<(OBJECT_DATA_SIZE);i++){
      if(average_flag[i] == AVERAGE){
	average_data[i] += object_data[i]*weight;
      }
      else if(average_flag[i] == TOTAL){
	average_data[i] += object_data[i];
      }
    }
    total_weight += weight;

    // mark as visited
    object_data_int[VISITED] = visited;
    fseeko64(input_mlist, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)molecule, SEEK_SET);
    fwrite(&object_data, sizeof(float), OBJECT_DATA_SIZE, input_mlist);  

    elements += 1;
  }

  // perform weighted averages
  for(i=0;i<(OBJECT_DATA_SIZE);i++){
    if(average_flag[i] == AVERAGE){
      average_data[i] = average_data[i]/total_weight;
    }
  }

  // save average object
  fwrite(&average_data, sizeof(float), OBJECT_DATA_SIZE, output_mlist);
}


/*
 * Main
 *
 *
 *
 */

int main(int argc, const char *argv[])
{
  int i, last_frame, molecules, tracks, unvisited;
  char header[DATA];
  int object_data[OBJECT_DATA_SIZE];
  FILE *input_mlist, *output_mlist;

  if (argc != 3){
    printf("usage avemlist <input file> <output file>\n");
    exit(0);
  }


  /* 
   * Setup 
   */

  input_mlist = fopen(argv[1], "rb+");
  if (!input_mlist){
    printf("avemlist: Could not open localization file %s\n", argv[1]);
    exit(0);
  }

  output_mlist = fopen(argv[2], "wb");
  if (!output_mlist){
    printf("avemlist: Could not open localization file %s\n", argv[1]);
    exit(0);
  }

  fread(&header, sizeof(char), DATA, input_mlist);
  fwrite(&header, sizeof(char), DATA, output_mlist);

  fseek(input_mlist, MOLECULES, SEEK_SET);
  fread(&molecules, sizeof(int), 1, input_mlist);
  // printf("Molecules: %d\n", molecules);

  fseek(input_mlist, DATA, SEEK_SET);
  fread(&object_data, sizeof(int), OBJECT_DATA_SIZE, input_mlist);
  unvisited = object_data[VISITED];
  // printf("Unvisited: %d\n", unvisited);

  /*
   * Go through all the molecules & generate averages.
   */
  last_frame = 0;
  tracks = 0;
  for(i=0;i<molecules;i++){
    if((i%50000)==0){
      printf("Processing molecule %d (avemlist)\n", i);
    }
    fseeko64(input_mlist, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)i, SEEK_SET);
    fread(&object_data, sizeof(int), OBJECT_DATA_SIZE, input_mlist);
    if (last_frame != object_data[FRAME]){
      fflush(input_mlist);
      last_frame = object_data[FRAME];
    }
    if (object_data[VISITED] == unvisited){
      if (object_data[CAT] >= 0){
	averageTrack(input_mlist, output_mlist, i, unvisited+1);
	tracks++;
      }
      else {
	object_data[VISITED] = unvisited+1;
	fseeko64(input_mlist, DATA + OBJECT_DATA_SIZE*DATUM_SIZE*(long long)i, SEEK_SET);
	fwrite(&object_data, sizeof(int), OBJECT_DATA_SIZE, input_mlist);
      }
    }
  }
  printf("Processed %d tracks\n", tracks);

  fseek(output_mlist, MOLECULES, SEEK_SET);
  fwrite(&tracks, sizeof(int), 1, output_mlist);

  fclose(input_mlist);
  fclose(output_mlist);
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
