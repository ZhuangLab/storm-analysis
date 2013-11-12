/*
 * Header file for that which is common to all the
 * homotopy based image analyzers.
 *
 * Hazen 04/13
 *
 */

/* Function Declarations */
void analyzeImage(double *, double *);
void closeFile(void);
void finishUp(void);
int getPeaks(double *, double *, double *, double *, int *, int);
int openFile(char *);
int saveHighRes(double *, int);
void setImageParameters(double *, int, int, int, double, int, int, int, int, int);

/* Global Variables */

int hres_x;
int hres_y;
int scale;

double epsilon;

FILE *hres_fp;

double *xvec;
double *yvec;

