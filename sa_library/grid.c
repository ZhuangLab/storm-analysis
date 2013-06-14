/*
 * Functions for gridding data.
 * Basically specializing numpy.histogramdd in the hopes of
 * speeding things up for "standard" cases.
 *
 * Hazen
 * 12/12
 * 
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall grid.c
 *  gcc -shared -Wl,-soname,grid.so.1 -o grid.so.1.0.1 grid.o -lc
 *
 * Windows:
 *  gcc -c grid.c
 *  gcc -shared -o grid.dll grid.o
 */

/* Include */
#include <stdlib.h>
#include <stdio.h>

void grid2D(int *, int *, int *, int, int, int);
void grid3D(int *, int *, int *, int *, int, int, int, int);
void grid3DZInclusive(int *, int *, int *, int *, int, int, int, int);

/*
 * grid2D()
 *
 * place x,y data into a 2D histogram.
 *
 * grid - pre-allocated storage for the histogram.
 * i_x - array of data x locations.
 * i_y - array of data y locations.
 * x_size - grid size in x.
 * y_size - grid size in y.
 * n_x - number of x (and y) locations.
 */
void grid2D(int *grid, int *i_x, int *i_y, int x_size, int y_size, int n_x)
{
  int i,x,y;

  for(i=0;i<n_x;i++){
    x = i_x[i];
    y = i_y[i];
    if((x>=0)&&(x<x_size)&&(y>=0)&&(y<y_size)){
      grid[x*y_size+y]++;
    }
  }
}

/*
 * grid3D()
 *
 * place x,y,z data into a 3D histogram.
 *
 * grid - pre-allocated storage for the histogram.
 * i_x - array of data x locations.
 * i_y - array of data y locations.
 * i_z - array of data z locations.
 * x_size - grid size in x.
 * y_size - grid size in y.
 * z_size - grid size in z.
 * n_x - number of x (and y,z) locations.
 */
void grid3D(int *grid, int *i_x, int *i_y, int *i_z, int x_size, int y_size, int z_size, int n_x)
{
  int i,x,y,z;

  for(i=0;i<n_x;i++){
    x = i_x[i];
    y = i_y[i];
    z = i_z[i];
    if((x>=0)&&(x<x_size)&&(y>=0)&&(y<y_size)&&(z>=0)&&(z<z_size)){
      grid[x*y_size*z_size+y*z_size+z]++;
    }
  }
}

/*
 * grid3DZInclusive()
 *
 * Place x,y,z data into a 3D histogram. Z values that are outside
 * of the bin range are placed in the first (or last) bin.
 *
 * grid - pre-allocated storage for the histogram.
 * i_x - array of data x locations.
 * i_y - array of data y locations.
 * i_z - array of data z locations.
 * x_size - grid size in x.
 * y_size - grid size in y.
 * z_size - grid size in z.
 * n_x - number of x (and y,z) locations.
 */
void grid3DInclusize(int *grid, int *i_x, int *i_y, int *i_z, int x_size, int y_size, int z_size, int n_x)
{
  int i,x,y,z;

  for(i=0;i<n_x;i++){
    x = i_x[i];
    y = i_y[i];
    z = i_z[i];
    if(z<0){
      z = 0;
    }
    else if(z>=z_size){
      z = z_size - 1;
    }
    if((x>=0)&&(x<x_size)&&(y>=0)&&(y<y_size)){
      grid[x*y_size*z_size+y*z_size+z]++;
    }
  }
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
