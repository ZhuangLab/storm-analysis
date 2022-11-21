/*
 * Header file for that which is common to all the
 * homotopy based image analyzers.
 *
 * Hazen 04/13
 *
 */

#ifndef HOMOTOPY_IMAGEA_COMMON_H
#define HOMOTOPY_IMAGEA_COMMON_H

/* Function Declarations */
void analyzeImage(double *, double *);
void closeFile(void);
void finishUp(void);
int getPeaks(double *, double *, double *, double *, int *, int);
int openFile(char *);
int saveHighRes(double *, int);
void setImageParameters(double *, int, int, int, double, int, int, int, int, int);

/* Global Variables */

static int hres_x;
static int hres_y;
static int scale;

static double epsilon;

static FILE *hres_fp;

static double *xvec;
static double *yvec;

#endif
