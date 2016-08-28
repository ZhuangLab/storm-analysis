/*
 * C library for calculating Zernike polynomials.
 * Has issues for polynomials where n >= 13?
 *
 * Hazen 10/14
 *
 * Compilation instructions:
 *
 * Linux:
 *  gcc -fPIC -g -c -Wall zernike.c
 *  gcc -shared -Wl,-soname,zernike.so.1 -o zernike.so.1.0.1 zernike.o -lc
 *  ln -s zernike.so.1.0.1 zernike.so
 *
 * Windows:
 *  gcc -c zernike.c
 *  gcc -shared -o zernike.dll zernike.o
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int factorial(int);
double pre_fac(int, int, int);
double zernike(int, int, double, double);
void zernike_grid(double *, int, double, double, double, int, int);
double zernike_rad(int, int, double);

/*
 * factorial()
 *
 * Calculate factorial of a integer.
 *
 * n - The input number.
 *
 * Returns n!
 */
int factorial(int n)
{
  int i, n_fac;

  n_fac = 1;
  for(i=1;i<=n;i++){
    n_fac = n_fac * i;
  }

  return n_fac;
}

/*
 * pre_fac()
 *
 * Calculate factorial coefficient.
 *
 * m - zernike m value.
 * n - zernike n value.
 * k - index.
 *
 * Returns the factorial coefficient.
 */
double pre_fac(int m, int n, int k)
{
  double d1, d2, d3, n1, sign;

  sign = pow(-1.0, k);
  n1 = (double)factorial(n-k);
  d1 = (double)factorial(k);
  d2 = (double)factorial((n+m)/2 - k);
  d3 = (double)factorial((n-m)/2 - k);

  return sign * n1/(d1 * d2 * d3);
}

/*
 * zernike()
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
double zernike(int m, int n, double rho, double phi)
{
  if (m > 0) return zernike_rad(m, n, rho) * cos(m * phi);
  if (m < 0) return zernike_rad(-m, n, rho) * sin(-m * phi);
  return zernike_rad(0, n, rho);
}

/*
 * zernike_grid()
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
void zernike_grid(double *grid, int gsize, double gcenter, double radius, double scale, int m, int n)
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
	grid[i*gsize+j] += scale * zernike(m, n, dd, phi);
      }
    }
  }
}

/*
 * zernike_rad()
 *
 * Calculate the radial component of a Zernike polynomial.
 *
 * m - m coefficient.
 * n - n coefficient.
 * rho - radius (0.0 - 1.0).
 *
 * Returns the radial value.
 */
double zernike_rad(int m, int n, double rho)
{
  int k;
  double sum;

  if((n < 0) || (m < 0) || (abs(m) > n)) return 0.0;
  if(((n-m)%2) == 1) return 0.0;

  sum = 0.0;
  for(k=0;k<((n-m)/2+1);k++){
    sum += pre_fac(m, n, k) * pow(rho, (n - 2*k));
  }
  
  return sum;
}
