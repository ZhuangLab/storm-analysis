/*
 * Common constants for multiple peak fitting.
 *
 * Hazen 10/17
 *
 */

/* debugging */
#define TESTING 0
#define VERBOSE 0

/* number of peak and results parameters. */
#define NFITTING 7
#define NPEAKPAR 9

/* indexes for peak fitting parameters. */
#define HEIGHT 0        /* Height */     
#define XCENTER 1       /* X center */
#define XWIDTH 2        /* Width in x, only relevant for gaussians */
#define YCENTER 3       /* Y center */
#define YWIDTH 4        /* Width in y, only relevant for gaussians */
#define BACKGROUND 5    /* Background level under the peak */
#define ZCENTER 6       /* Z center */

/* additional indexes for results. */
#define STATUS 7        /* Status flag, see below */
#define IERROR 8        /* Error in the fit (integrated over the AOI) */

/* peak status flags */
#define RUNNING 0
#define CONVERGED 1
#define ERROR 2

#define HYSTERESIS 0.6 /* In order to move the AOI or change it's size,
			  the new value must differ from the old value
			  by at least this much (<= 0.5 is no hysteresis). */

#define MAXCYCLES 10 /* The maximum number of times to increase lambda to try
                        and get a fit that reduces the peak error. */

/*
 * There is one of these for each peak to be fit.
 */
typedef struct peakData
{
  int index;                /* peak id */
  int status;               /* status of the fit (running, converged, etc.). */
  int xi;                   /* location of the fitting area in x. */
  int yi;                   /* location of the fitting area in y. */

  int size_x;               /* size of the fitting area in x. */
  int size_y;               /* size of the fitting area in y. */

  double error;             /* current error. */
  double error_old;         /* error during previous fitting cycle. */

  double lambda;            /* Levenberg-Marquadt lambda term. */

  int sign[NFITTING];       /* sign of the (previous) update vector. */
  double clamp[NFITTING];   /* clamp term to suppress fit oscillations. */
  double params[NFITTING];  /* [height x-center x-width y-center y-width background] */  

  void *peak_model;         /* Pointer to peak model specific data (i.e. spline data, etc.) */
} peakData;


/*
 * This structure contains everything necessary to fit
 * an array of peaks on an image.
 */
typedef struct fitData
{
  /* These are for diagnostics. */
  int n_dposv;                  /* number lost to an error trying to solve Ax = b. */
  int n_iterations;             /* number of iterations of fitting. */
  int n_margin;                 /* number lost because they were too close to the edge of the image. */
  int n_neg_fi;                 /* number lost to a negative fi. */
  int n_neg_height;             /* number lost to negative height. */
  int n_neg_width;              /* number lost to negative width. */

  int jac_size;                 /* The number of terms in the Jacobian. */
  int margin;                   /* size of the band around the edge of the image to avoid. */
  int nfit;                     /* number of peaks to fit. */
  int image_size_x;             /* size in x (fast axis). */
  int image_size_y;             /* size in y (slow axis). */

  double xoff;                  /* offset between the peak center parameter in x and the actual center. */
  double yoff;                  /* offset between the peak center parameter in y and the actual center. */
  double zoff;                  /* offset between the peak center parameter in z and the actual center. */

  double tolerance;             /* fit tolerance. */

  int *bg_counts;               /* number of peaks covering a particular pixel. */
  
  double *bg_data;              /* background data. */
  double *f_data;               /* fit (foreground) data. */
  double *scmos_term;           /* sCMOS calibration term for each pixel (var/gain^2). */
  double *x_data;               /* image data. */

  double clamp_start[NFITTING]; /* starting values for the peak clamp values. */

  peakData *working_peak;       /* Working copy of the peak that we are trying to improve the fit of. */
  peakData *fit;                /* The peaks to be fit to the image. */

  void *fit_model;              /* Other data/structures necessary to do the fitting, such as a cubic spline structure. */

  /*
   * Specific fitter versions must provide these functions.
   */
  void (*fn_add_peak)(struct fitData *);                      /* Function for adding working peak to the fit image. */
  void (*fn_calc_JH)(struct fitData *, double *, double *);   /* Function for calculating the Jacobian and the Hessian. */
  int (*fn_check)(struct fitData *);                          /* Function for checking the validity of the working peak parameters. */
  void (*fn_copy_peak)(struct peakData *, struct peakData *); /* Function for copying peaks. */
  void (*fn_subtract_peak)(struct fitData *);                 /* Function for subtracting the working peak from the fit image. */
  void (*fn_update)(struct fitData *, double *);              /* Function for updating the working peak parameters. */
  
} fitData;


/*
 * Functions.
 *
 * These all start with mFit to make it easier to figure
 * out in libraries that use this code where they came from.
 */
int mFitCalcErr(fitData *);
void mFitCleanup(fitData *);
void mFitCopyPeak(peakData *, peakData *);
void mFitGetFitImage(fitData *, double *);
void mFitGetResidual(fitData *, double *);
void mFitGetResults(fitData *, double *);
int mFitGetUnconverged(fitData *);
fitData *mFitInitialize(double *, double *, double, int, int);
void mFitIterate(fitData *);
void mFitNewImage(fitData *, double *);
void mFitNewPeaks(fitData *, double *, int);
int mFitSolve(double *, double *, int);
void mFitUpdateParam(peakData *, double *, int);
