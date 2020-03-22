#include "common_structures_strat.h"

/* determine the version between mac and linux */

/*  mac  */
#define _malloc malloc
#define _realloc realloc

/*
#define linux 1
#define _malloc GC_MALLOC
#define _realloc GC_REALLOC
*/

/* Define the physical constants used in the new_tephra.c code
 * units are SI except WIND_INTERVAL (in km)
 * MAX_LINE, COL_STEPS and PART_STEPS dimensionless
*/
#define LOG_FILE "node.txt"
/* #define _PRINT 1 */
#define pi 3.141592654
#define DEG2RAD 0.017453293
/*#define DEBUG 0 */
#define MAX_LINE 200
#define MAX_LINE2 10000

/*mixed diffusion model*/
/*eddy diff for small particles in m2/s (400 cm2/s) */
extern double EDDY_CONST;

/* diffusion coeff for large particles (m2/s) */
extern double DIFFUSION_COEFFICIENT; 

/* diffusion time (s) */
extern double SIGMA_PLUME;   /* added by Kaz Dec-01-'10*/

/*threshold for change in diffusion (seconds fall time) */
extern double FALL_TIME_THRESHOLD; 


/* density of air */
/* air density in kg/m3 */
#define AIR_DENSITY 1.293
   
/* dynamic viscosity of air */
#define AIR_VISCOSITY 0.000018325
#define GRAVITY 9.81

/*density model for the pyroclasts */
/* These are now defined in config file */
extern double LITHIC_DENSITY;
extern double PUMICE_DENSITY;

#define LITHIC_DIAMETER_THRESHOLD 7.0
#define PUMICE_DIAMETER_THRESHOLD -1.0

/* #define BETA_LIMIT 100.0, now same as COL_STEPS */

enum{DISTR_1, DISTR_2};

/* These are now defined in config file. */
extern int COL_STEPS;
extern int GRAIN_DATA;
extern int PLUME_MODEL;
extern int PART_STEPS;

extern double PLUME_RATIO; /* Hb/Ht[area_of_release] of the laterally spreading cloud */
extern double SUZUKI_A;    /* Suzuki parameters added on 11-JAN-2011 */
extern double SUZUKI_LAMDA;
extern int WIND_DAYS;
extern int WIND_COLUMNS;

// Added due to annexation of windy.c (Feb, 2017)
extern double INITIAL_WATER_CONTENT;
extern double VENT_ELEVATION; // Original parameter in Tephra2 but not global
extern double MAGMA_DISCHARGE_RATE;
extern double MAGMA_TEMPERATURE;
extern double INITIAL_PLUME_VELOCITY;
extern double VENT_RADIUS;
extern double PLUME_HEIGHT;
extern double SLICE_HEIGHT;

extern double ds;
extern double S_MAX;
extern int S_STEPS;
extern double S_HT_MAX;
extern int S_HT_STEPS;

/* wind data read in at 0.5km intervals - this can be changed to suit available data
 * make sure input file intervals agree */


