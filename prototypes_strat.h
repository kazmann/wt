#ifdef linux
#include "parameters_strat.h"
#else
#include "parameters_strat_mac.h"
#endif

/* define the following functions */

int get_eruptions(void);
int get_points(FILE *in_points);
double get_grain(FILE *in, int *gl_ptr);    	 	/* added by K. Mannen (19-Jan-2011) modified 2014.10.05 */
double get_rho(double part_size, int *gl_ptr);		/* added by K. Mannen (19-Jan-2011) modified 2014.10.05 */
int get_wind(FILE *in_wind);
int get_wind2(FILE *in_wind);			/* added by K. Mannen (19-Jan-2011) */
int init_globals(char *in_config);
void set_global_values(FILE *log_file);
void set_eruption_values(ERUPTION *erupt, WIND *wind, int *gl_ptr); 	//2014.10.05

void free_memory_tephra_calc(void);	//2018.02.18
void free_memory_new_calc(void); 	//2018.03.18
	
double phi2m( double xx );
double particle_density (double phi_slice);
double part_diff_time( double col_ht );
double (*pdf)(double, int, double, double);
double plume_pdf0(double x, int slice, double none, double ht_section_width);
double plume_pdf1(double x, int slice, double column_height, double vent_elevation);
double part_fall_time(double col_ht, double layer, double ashdiam, double part_density);
double part_fall_time_vg(double particle_ht, double ground_alt, WIND *level, double *wind_element, double ashdiam, double part_density);/* added by K. Mannen (19-Jan-2011) */
double pdf_grainsize(double part_mean_size, double part_size_slice, double part_step_width, int *gl_ptr);  	//2014.10.05


//double strat_average(WIND *level, double col_ht, double xspace, double yspace, double fall_time, double sigma, double elev);
double strat_average( double average_wind_direction, 
                      double average_wind_speed,             
			                double xspace, double yspace, 
			                double total_fall_time,
			                double sigma);
void tephra_calc(ERUPTION *erupt, POINT *pt, WIND *level, STATS *stats);
double part_fall_time_vent2ground (double vent_height, double layer,  double ashdiam,double part_density);

double get_dir2(double level, double ht1, double ht0, double dir1, double dir0);	// bug fix related wind dir interpolation across 0 deg. (27-Feb-2017)

// Prototypes related to windy.c/////
// ## of func## means equation number appears in Woodhouse et al. (2013) JGR DOI 10.1029/2012JB009592
void windy(FILE *in_wind);
void rk(int);
int read_wind(FILE *in_wind);

double func12(double, double, double, double);
double func13(double, double, double, double, double);
double func14(double, double, double, double, double);
double func15(double, double, double, double, double, double, double);
double func16(double, double, double, double);
	
double func17(double, double, double, double);				// calc column density
double func18(double);				// particle content
double func19(double);				// Rg calc
double func20(double);				// Cp calc

double calcCp0(void);			// Cp0 calc (eq. 20.5; in text between eq20 and 21)
double func21(double);			// temperature profile
double calc_Tatm(double, int);
double calc_Patm(double, int);
double func22(double, double);	// pressure profile
double func23(double, double);	// atmospheric density
double func24(double);			// wind velocity

double get_V(double, int);
double get_dir(double, int);

void read_column_file(int);	// main outputs column.txt
								// and this function read it
//////////////////////////////////