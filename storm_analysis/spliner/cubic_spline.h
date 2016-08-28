/*
 * API functions.
 *
 * Hazen 01/14
 *
 */

/* Function Declarations */
void computeDelta2D(double, double);
void computeDelta3D(double, double, double);

double dxfAt2D(int, int);
double dxfAt3D(int, int, int);
double dxfSpline2D(double, double);
double dxfSpline3D(double, double, double);

double dyfAt2D(int, int);
double dyfAt3D(int, int, int);
double dyfSpline2D(double, double);
double dyfSpline3D(double, double, double);

double dzfAt3D(int, int, int);
double dzfSpline3D(double, double, double);

double fAt2D(int, int);
double fAt3D(int, int, int);
double fSpline2D(double, double);
double fSpline3D(double, double, double);

int getXSize(void);
int getYSize(void);
int getZSize(void);

void initSpline2D(double *, int, int);
void initSpline3D(double *, int, int, int);

void splineCleanup(void);

