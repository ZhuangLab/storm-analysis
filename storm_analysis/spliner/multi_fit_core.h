/* 
 * Header file for the multiple peak fitter.
 *
 * Hazen 12/13
 */

/*
 * begin FIXME: 
 *
 * This part is just a copy of multi_fit.h.
 */

/* number of peak fitting parameters. */
#define NPEAKPAR 7

/* indexs for peak parameters. */
#define HEIGHT 0      
#define XCENTER 1
#define XWIDTH 2
#define YCENTER 3
#define YWIDTH 4
#define BACKGROUND 5
#define ZCENTER 6

/* number of parameters in the results array. */
#define NRESULTSPAR 9

/* additional indexs for results. */
#define STATUS 7
#define IERROR 8

/* peak status flags */
#define RUNNING 0
#define CONVERGED 1
#define CHOLERROR 2
#define BADPEAK 3

/*
 * end FIXME
 */

#define MAXPARAMS 6

/* Structures */
typedef struct
{
  int index;         /* Fit index or id number. */
  int n_params;      /* Number of parameters that define the fit. */
  int size_x;        /* Size of the fitting area in x. */
  int size_y;        /* Size of the fitting area in y. */
  int size_z;        /* Z range. */
  int status;        /* Status of the fit (running, converged, etc.). */
  int type;          /* Type of peak. */
  int xi;            /* X index of the upper left corner of the fit area. */
  int yi;            /* Y index of the upper left corner of the fit area. */
  int zi;            /* Z index. */ 
  int *sign;         /* Sign of the (previous) update vector. */
  double error;      /* Current error. */
  double lambda;     /* Lambda term in the LM fit. */
  double *clamp;     /* Clamp term to supress fit oscillations. */
  double *delta;     /* Parameter update vector. */
  double *params;    /* Current fit parameters. */
  double *values;    /* Current fit values (and all their derivatives) in the fit area. */
} fitData;

/* (API) Functions */
void addPeak(fitData *);
double calcError(fitData *, double);
void copyFitData(fitData *, fitData *);
void getResidual(double *);
void freeFitData(fitData *);
void initializeMultiFit(double *, double, int, int);
void mallocFitData(fitData *, int, int);
void multiFitCleanup(void);
void newImage(double *);
void resetFit(void);
void subtractPeak(fitData *);
void updateFit(fitData *, double);

/* Global Variables */
extern int image_size_x;    /* Image size in x. */
extern int image_size_y;    /* Image size in y. */
extern double tolerance;    /* Fit tolerance for convergence. */
