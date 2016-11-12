/*
 * C library for rolling ball background removal.
 *
 * 02/16
 *
 * Hazen 
 *
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall rolling_ball_lib.c
 *  gcc -shared -Wl,-soname,rolling_ball_lib.so.1 -o rolling_ball_lib.so.1.0.1 rolling_ball_lib.o -lc
 *
 * Windows (mingw):
 *  gcc -c rolling_ball_lib.c
 *  gcc -shared -o rolling_ball_lib.dll rolling_ball_lib.o
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>


/* Structures */
typedef struct{
  int ball_size;
  double *ball;
} ballData;

/* Function Declarations */
void cleanup(ballData *);
void estimateBg(ballData *, double *, double *, int, int);
ballData* init(double *, int);


/*
 * cleanup()
 *
 * ball_data - Pointer to a ballData structure.
 */
void cleanup(ballData *ball_data)
{
  free(ball_data->ball);
  free(ball_data);
}

/*
 * estimateBg()
 *
 * ball_data - Pointer to a ballData structure.
 * image - The image to estimate the background of.
 * background - Pre-initialized storage for the background estimate.
 * image_x - The size of the image in x (slow dimension).
 * image_y - The size of the image in y (fast dimension).
 */
void estimateBg(ballData *ball_data, double *image, double *background, int image_x, int image_y)
{
  int bb,cx,cy,i,j,k,l;
  double min,cur;

  bb = (ball_data->ball_size - 1)/2;
  
  for(i=0;i<image_x;i++){
    for(j=0;j<image_y;j++){
      min = image[i*image_y+j];
      for(k=0;k<ball_data->ball_size;k++){
	cx = i + k - bb;
	if (cx < 0) continue;
	if (cx >= image_x) continue;
	for(l=0;l<ball_data->ball_size;l++){
	  cy = j + l - bb;
	  if (cy < 0) continue;
	  if (cy >= image_y) continue;
	  cur = image[cx*image_y+cy] - ball_data->ball[k*ball_data->ball_size+l];
	  if (cur < min){
	    min = cur;
	  }
	}
      }
      background[i*image_y+j] = min;
    }
  }
}

/*
 * init()
 *
 * py_ball - The python ball array (square).
 * py_ball_size - The size of py_ball (should be an odd number).
 */
ballData *init(double *py_ball, int py_ball_size)
{
  int i;
  ballData *ball_data;

  ball_data = (ballData *)malloc(sizeof(ballData));
  
  ball_data->ball_size = py_ball_size;

  ball_data->ball = (double *)malloc(sizeof(double) * ball_data->ball_size * ball_data->ball_size);

  for (i=0;i<(ball_data->ball_size * ball_data->ball_size);i++){
    ball_data->ball[i] = py_ball[i];
  }

  return ball_data;
}
