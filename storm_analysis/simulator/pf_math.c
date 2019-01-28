/*
 * C library for doing some of heavy lifting of pupil function 
 * calculations.
 *
 * Hazen 01/19
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int pfFactorial(int);
double pfPreFac(int, int, int);
double pfZernike(int, int, double, double);
void pfZernikeGrid(double *, int, double, double, double, int, int);
double pfZernikeRad(int, int, double);


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
