/*
 * 08/11
 * 
 * Common constants for multiple peak fitting.
 *
 *
 * Hazen
 *
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
#define ERROR 2
#define BADPEAK 3



