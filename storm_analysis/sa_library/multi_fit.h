/*
 * Common constants for multiple peak fitting.
 *
 * Hazen 10/16
 *
 */

/* debugging */
#define TESTING 0
#define VERBOSE 0

/* number of peak and results parameters. */
#define NFITTING 7
#define NPEAKPAR 9

/* indexes for peak fitting parameters. */
#define HEIGHT 0      
#define XCENTER 1
#define XWIDTH 2
#define YCENTER 3
#define YWIDTH 4
#define BACKGROUND 5
#define ZCENTER 6

/* additional indexes for results. */
#define STATUS 7
#define IERROR 8

/* peak status flags */
#define RUNNING 0
#define CONVERGED 1
#define ERROR 2
#define BADPEAK 3


/*
 * There is one of these for each peak to be fit.
 */
typedef struct
{
  int index;                /* peak id */
  int size_x;               /* size of the fitting area in x. */
  int size_y;               /* size of the fitting area in y. */
  int status;               /* status of the fit (running, converged, etc.). */
  int xi;                   /* (integer) location of the fitting area in x. */
  int yi;                   /* (integer) location of the fitting area in y. */
  
  double error;             /* current error. */
  double error_old;         /* error during previous fitting cycle. */

  int sign[NFITTING];       /* sign of the (previous) update vector. */
  double clamp[NFITTING];   /* clamp term to suppress fit oscillations. */
  double params[NFITTING];  /* [height x-center x-width y-center y-width background] */  

  void *peakModel;          /* Pointer to peak model specific data (i.e. spline data, etc.) */
} peakData;


/*
 * This structure contains everything necessary to fit
 * an array of peaks on an image.
 */
typedef struct
{
  /* These are for diagnostics. */
  int n_dposv;                  /* number lost to an error trying to solve Ax = b. */
  int n_margin;                 /* number lost because they were too close to the edge of the image. */
  int n_neg_fi;                 /* number lost to a negative fi. */
  int n_neg_height;             /* number lost to negative height. */
  int n_neg_width;              /* number lost to negative width. */


  int margin;                   /* size of the band around the edge of the image to avoid. */
  int nfit;                     /* number of peaks to fit. */
  int image_size_x;             /* size in x (fast axis). */
  int image_size_y;             /* size in y (slow axis). */

  double xoff;                  /* offset between the peak center parameter in x and the actual center. */
  double yoff;                  /* offset between the peak center parameter in y and the actual center. */
  double zoff;                  /* offset between the peak center parameter in z and the actual center. */

  double tolerance;             /* fit tolerance. */
  double min_z;                 /* minimum z value. */
  double max_z;                 /* maximum z value. */

  int *bg_counts;               /* number of peaks covering a particular pixel. */
  
  double *bg_data;              /* background data. */
  double *f_data;               /* fit (foreground) data. */
  double *scmos_term;           /* sCMOS calibration term for each pixel (var/gain^2). */
  double *x_data;               /* image data. */

  double clamp_start[NFITTING]; /* starting values for the peak clamp values. */

  peakData *fit;                /* The peaks to be fit to the image. */
  void *fitModel;               /* Other data/structures necessary to do the fitting, such as a cubic spline structure. */
  
} fitData;


/*
 * Functions.
 */
void calcErr(fitData *, peakData *);
void getResidual(fitData *, double *);
void getResults(fitData *, double *);
int getUnconverged(fitData *);
void newImage(fitData *, double *);
void updateParams(peakData *, double *);
