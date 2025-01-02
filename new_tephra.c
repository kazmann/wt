#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef linux
#include <gc.h>
#endif

#include "prototypes_strat.h"

/*
Code: new_tephra.c
By: C.B. & L.J. Connor, T. Hincks, and C. Bonadonna
Copyright (C) 2003  C.B. Connor, L.J. Connor, C. Bonadonna, T. Hincks
See: http://www.cas.usf.edu/~cconnor/parallel/tephra/tephra.html

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/

/* The following Global Variables are assigned some default values */

double DIFFUSION_COEFFICIENT = 200.0;
double SIGMA_PLUME = 0;				/* added by Kaz Dec-02-'10*/
double FALL_TIME_THRESHOLD = 180.0;
double EDDY_CONST = .04;
double LITHIC_DENSITY = 2350.0;
double PUMICE_DENSITY = 1000.0;

int COL_STEPS = -9999;
int GRAIN_COLUMN = 2;	/* 18-Jan-2011 added by Kaz */
int GRAIN_DATA =0;		/* 18-Jan-2011 added by Kaz */
int PLUME_MODEL = DISTR_1;
int PART_STEPS = 100;
int WIND_DAYS = 1;
int WIND_COLUMNS = 3;

double PLUME_RATIO = 0.1;
double WIND_INTERVAL;
double PLUME_HEIGHT = -9999;
double PLUME_THICKNESS = -9999;        /* added by Kaz 09-Mar-2020 */
double SOURCE_DECAY_RATE = -9999;      /* added by Kaz 16-Sep-2021 */
double ERUPTION_MASS = 1e10;
double MAX_GRAINSIZE = -7.0;
double MIN_GRAINSIZE = 7.0;
double MEDIAN_GRAINSIZE = 1.5;
double STD_GRAINSIZE = 2.0;
double SUZUKI_A =1;
double SUZUKI_LAMDA =1;
double VENT_EASTING = 0.0;
double VENT_NORTHING = 0.0;
double VENT_ELEVATION = 0.0;

double INITIAL_WATER_CONTENT = -9999;

double MAGMA_DISCHARGE_RATE = -9999;
double MAGMA_TEMPERATURE = 1200;
double INITIAL_PLUME_VELOCITY = -9999;
double VENT_RADIUS = -9999;

double S_MAX= -9999;
int S_STEPS= -9999;
double S_HT_MAX = -9999;
int S_HT_STEPS =  -9999;

double KS_ENTRAIN_U = 0.09;	// Added on 27 Oct 2021
double KW_ENTRAIN_V = 0.9;	// Added on 27 Oct 2021


/*define the following data structures for the code
  full descriptions are found in common_structures.h */
static ERUPTION *erupt;
static GRAIN *g_ptr;	/* 18-Jan-2011 added by Kaz */
static WIND **W;
static WIND **W2;		/* 18-Jan-2011 added by Kaz; wind condition between vent height and sea level */
static POINT *pt;


/*define the following variables passed among functions in this file */
static int num_pts = 0; /*total number of points used in the analysis */
static int num_eruptions = 0; /*total number of eruptions used in the analysis */
/* static int num_wind_data = 0; number of wind data points in the analysis */

static int local_n;

//FILE *in_eruptions;
FILE *in_points;
FILE *in_wind;
FILE *log_file;
FILE *in_grain;

int cmp(double *x, double *y) {
  if (*x < *y) return -1;
  if (*x > *y) return 1;
  return 0;
}

void exit_now(int e) {
  
  //(void) fclose(in_eruptions);
  (void) fclose(in_points);
  (void) fclose(in_wind);
  
  (void) fclose(in_grain);	// 2018.03.18
#ifdef _PRINT
  (void) fclose(log_file);
#endif
#ifdef DEBUG
  (void) fclose(log_file);
#endif
  exit(e);
}

int main(int argc, char *argv[]) { /* MAIN CODE STARTS HERE */
  
  int i, j, bin, phi_bins;
  int g_ttl_num = 0;	// 2014.10.05
  int *gl_ptr;			// 2014.10.05
  
  gl_ptr = &g_ttl_num;	// 2014.10.05
  
  double val;
  STATS stats;
  char log_name[25];


  /* Check for correct number of comand line arguments */
  if (argc < 4) {   /* modified by K. Mannen (19-Jan-2011) */
      fprintf(stderr, 
	      "Missing comand line arguments,\nUSAGE: <program name> <config file> <points file> <wind file> (<grain file> option)\n\n");
    exit(1);
  }
  

  /* Each node opens a file for logging */
  sprintf(log_name, "%s", LOG_FILE);
  
  fprintf(stderr, "%s\n", log_name);
  
  log_file  = fopen(log_name, "w+");
  if (log_file == NULL) {
    fprintf(stderr, "Cannot open LOG file=[%s]:[%s]. Exiting.\n", 
	    log_name, strerror(errno));
    exit(1);
  }

  
  /* Initialize the global variables (see top of file) with inputs from the configuration file. */
  if ( init_globals(argv[1]) ) {
    exit(1);
  }
  

  
  /*make sure the eruptions file exists 
  in_eruptions = fopen(argv[2], "r");
  if (in_eruptions == NULL) {

    fprintf(stderr, "Cannot open eruptions  file=[%s]:[%s]. Exiting.\n", 
	    argv[2], strerror(errno));
    exit_now(1);
  }
  */
  
#ifdef _PRINT
  fflush(log_file); 
#endif

  /*make sure the points file exists*/
  in_points= fopen(argv[2], "r");
  if (in_points == NULL) {
    fprintf(stderr, "Cannot open points  file=[%s]:[%s]. Exiting.\n", 
	    argv[2], strerror(errno));
    exit_now(1);
  }
  
  /* Input the data points from a file using the
     get_points function. 
  */
  
  if (get_points(in_points) ) {
    exit_now(1);
  }

#ifdef _PRINT
  fflush(log_file); 
#endif
  
  /*make sure the wind file exists*/
  in_wind= fopen(argv[3], "r");
  if (in_wind == NULL) {
    fprintf(stderr, "Cannot open wind file=[%s]:[%s]. Exiting.\n", 
	    argv[3], strerror(errno));
    exit_now(1);
  }
  
  windy(in_wind);	// Feb. 6, 2017 by kaz
  
  if (get_wind2(in_wind) ) {
    exit_now(1);
  }
	
  if (get_wind(in_wind) ) {
    exit_now(1);
  }
	
	/* read grain size file by K. Mannen (19-Jan-2011) */
	if (argc == 5)
	{
		in_grain= fopen(argv[4], "r");
		GRAIN_DATA = 1;
		if (in_grain == NULL)
		{
			fprintf(stderr, "Cannot open grain file=[%s]:[%s]. Exiting.\n", argv[4], strerror(errno));
			exit_now(1);
		}
		if (get_grain(in_grain, gl_ptr))
		{
			exit_now(1);
		}
	}
	
	
  /* 
     
  Note: "local_n" is the number of points from the
  input file. Therefore local_n is
  declared a static variable and is assigned a value in the
  get_points function 
  
  Note: tephra_calc is the main computational part of the
  code.
  
  Input the data points from wind file using the get_wind function. 
  
 
  */
  
#ifdef _PRINT
  fflush(log_file); 
#endif
  
  set_global_values(log_file);
  
#ifdef _PRINT
  fflush(log_file); 
#endif
 
/* each node gets all of the eruption data */ 
  if (get_eruptions() ) {
    exit_now(1);
  }

/* Calculating an accumulation map */
    
  for (j=0;  j < num_eruptions; j++) { /* For each eruptive senario */
    set_eruption_values(erupt+j, W[j], gl_ptr);
    //fprintf(stderr, "[%d]PARTICLE STEPS=%d ", j, PART_STEPS);
      for (i = 0;i < local_n; i++) {  /* For each location */
	      /* Note: W[j]  means that if there are multiple eruptions, 
	      there should be multiple WIND_DAYS in the wind.in file, 
	      1 WIND_DAY for each eruption line */
			tephra_calc(erupt+j, pt+i, W2[j], &stats);  /* W[j] -> W2[j]; modified by K. Mannen (19-Jan-2011) */ 
	      (pt+i)->cum_mass += (pt+i)->mass;
      }
    } j--;
    fprintf(stderr,"\nMin Particle Fall Time = %gs\nMax Particle Fall Time = %gs\n",stats.min_falltime, stats.max_falltime);
    
	/* Here I am assuming that each eruptive senario has the same min and max particle size range
	   So, the grainsize distribution contains the same number of bins. The amount in each bin accumulates
	   between eruptive senarios.
	*/
	phi_bins = (int)((erupt+j)->max_phi - (erupt+j)->min_phi);
	printf("#LOC\tEAST\tNORTH\tELEV\tMASS");
	for (bin = 0; bin < phi_bins; bin++) 
	  printf("\t[%d->%d)", 
		(int)((erupt+j)->min_phi)+bin, 
		(int)((erupt+j)->min_phi)+bin+1);
	
	printf("\n");
	
	fprintf(stderr, "PART_STEPS=%d phi_bins=%d\n", PART_STEPS, phi_bins);
	for (i=0; i < num_pts; i++) {
	  printf("%s\t%.0f\t%.0f\t%.0f\t%.6e", 
		(pt+i)->pnt_name,
		(pt+i)->easting, 
		(pt+i)->northing,
		(pt+i)->elevation, 
		(pt+i)->cum_mass);
	  for (bin=0; bin < phi_bins; bin++) {
	  	val = (pt+i)->phi[bin];
	    printf("\t%.6e", val);
	  }
	  printf("\n");
	}



  fprintf(log_file, "Finished.\n");
  
  free(erupt);
  free(g_ptr);
  free(pt);
  free(W);
  free(W2);
  
  free_memory_tephra_calc();
  free_memory_new_calc();

  exit_now(0);
  return 1;
}

/**************************************************************
FUNCTION:  get_eruptions
DESCRIPTION:  This function reads eruption data into the
ERUPTION array. Each node stores 
all of the eruption  parameters which are varied and then 
used in calculating the mass loading value at each point. 
INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/
int get_eruptions(void) {
  
  //char line[MAX_LINE];
  int i;
  
#ifdef _PRINT
  fprintf(log_file,"ENTER[get_eruptions]\n");
#endif
  
  /* while (fgets(line, MAX_LINE, in) != NULL) {
    if (line[0] == '#' || line[0] == '\n') continue;
    num_eruptions++;
  } */
	num_eruptions = 1;

	erupt = (ERUPTION *)malloc((size_t)(num_eruptions +1) * sizeof(ERUPTION));


  
  if (erupt == NULL) {
    fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for eruptions:[%s]\n",
            0,0, strerror(errno));
    return -1;
  } 
  
#ifdef _PRINT 
  fprintf(log_file,"\t%d eruptions in file.\n", num_eruptions);
#endif  
  
  //rewind(in);

  i=0;
 /* while (fgets(line, MAX_LINE, in) != NULL) {
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      while (ret = sscanf(line,
			  "%lf %lf %lf %lf %lf %lf %lf %lf %lf",*/
			  (erupt+i)->volcano_easting = VENT_EASTING;
			  (erupt+i)->volcano_northing = VENT_NORTHING;
			  (erupt+i)->total_ash_mass = ERUPTION_MASS;
			  (erupt+i)->min_phi = MAX_GRAINSIZE;
			  (erupt+i)->max_phi= MIN_GRAINSIZE;
			  (erupt+i)->mean_phi = MEDIAN_GRAINSIZE;
			  (erupt+i)->sigma_phi= STD_GRAINSIZE;
			  (erupt+i)->vent_height = VENT_ELEVATION;	
			  (erupt+i)->max_plume_height = PLUME_HEIGHT;
	
//), ret != 9) { 
// &(erupt+i)->column_beta,
/*	
	if (ret == EOF && errno == EINTR) continue;
	fprintf(stderr, "[line=%d,ret=%d] Did not read in 9 parameters:[%s]\n", i+1,ret, strerror(errno));
	return -1;
      }
      i++;
    }
  }
 */ 
#ifdef _PRINT
  fprintf(log_file, "EXIT[get_eruptions].\n");	  
#endif
  return 0;
}

/**************************************************************
FUNCTION:  get_grain
DESCRIPTION:  This function reads grain data prepared by user
G, which is grain size distribution, is stored here

INPUTS: (IN) FILE *in  (file handle from which to read)

***************************************************************/

double get_grain(FILE *in, int *gl_ptr)
{
	int i=0, ret;
	int total_line_number=0;
	double grain_size, grain_cont;
	char line[MAX_LINE];
	

	
/* count total number of grain size data */	
	while (NULL != fgets(line, MAX_LINE, in))
	{
		if (line[0] == '#' || strlen(line) < GRAIN_COLUMN) continue;
		++total_line_number;
		/*printf("line = %s", line);*/
	}
	
	++total_line_number; //2014.10.04
	
	/*struct GRAIN G[total_line_number];*/
	/*struct GRAIN *g_ptr;*/
	
	g_ptr = (GRAIN *)malloc((size_t)total_line_number * sizeof(GRAIN));
	if (g_ptr == NULL)
	{
		fprintf(log_file, 
				"Cannot malloc memory for Integration Table:[%s]\n",
				strerror(errno));
		exit(1);
    }
	
	
	i=0;
	rewind(in);
	
	while (NULL != fgets(line, MAX_LINE, in))
	{
		if (line[0] == '#' || strlen(line) < GRAIN_COLUMN) continue;
		else{
			while (ret = sscanf(line,"%lf %lf",&grain_size,&grain_cont),ret != 2)
			{
			if (ret == EOF && errno == EINTR) continue;
	        
	        fprintf(stderr, "[line=%d,ret=%d] Did not read in 3 parameters:[%s]\n", 
	        i+1,ret, strerror(errno));
	        
	        return -1;
			}
		}
		g_ptr[i].sizeclass = grain_size;
		g_ptr[i].content   = grain_cont;

		/*printf("line 115, line = %s sizeclass = %1.1f sizecontent = %1.4f\n", line, g_ptr[i].sizeclass, g_ptr[i].content);
		*/
		++i;
	}
	/*printf("total line number = %d\n", total_line_number);*/
	*gl_ptr = i;
	//printf("total num of grain file in line 431 is %d\n\n", *gl_ptr);
	return 0;
}


/**************************************************************
FUNCTION:  get_points
DESCRIPTION:  This function reads eruption data into the
ERUPTION array. Each node stores 
all of the eruption  parameters which are varied and then 
used in calculating the mass loading value at each point. 
INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/

int get_points(FILE *in) {
  
  char line[MAX_LINE];
  int i, j, ret=0, my_start=0, pts_read=0;
  
#ifdef _PRINT
  fprintf(log_file,"ENTER[get_points]\n");
#endif
  while (fgets(line, MAX_LINE, in) != NULL)  {
    if (line[0] == '#' || line[0] == '\n') continue;
    num_pts++;
  }
  rewind(in);
  

  local_n = num_pts;
  

#ifdef _PRINT  
    fprintf(log_file, "Total locations: %d.\n", num_pts);
#endif    
	

	pt = (POINT *)malloc((size_t)local_n * sizeof(POINT));

	
  if (pt == NULL) {
    fprintf(stderr, "Cannot malloc memory for my points:[%s]\n", strerror(errno));
    return -1;
  } 
  
#ifdef _PRINT  
  fprintf(log_file, "\tReading in %d locations, starting at line %d.\n",
	  local_n, my_start);
#endif    
   
  i=0;
  while (i < num_pts) {
    fgets(line, MAX_LINE, in);
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
	  strcpy((pt+pts_read)->pnt_name, "9999");
      while (ret = sscanf(line,
			  "%lf %lf %lf %s",
			  &(pt+pts_read)->easting,
			  &(pt+pts_read)->northing,
			  &(pt+pts_read)->elevation,
			  (pt+pts_read)->pnt_name),
	      (ret < 3)) { 
	
	if (ret == EOF && errno == EINTR) continue;
	fprintf(stderr, "[line=%d,ret=%d] Did not read in 3 coordinates:[%s]\n", i+1,ret, strerror(errno));
	return -1;
      }
      /* Initialize some values in the point structure */
      (pt+pts_read)->cum_mass = 0.0;
      (pt+pts_read)->mass_pt = NULL;
	  
      for (j=0; j<20; j++)
	(pt+pts_read)->phi[j] = 0.0;
      if (i >= my_start) {
	pts_read++;
	if (pts_read == local_n) break;
      } 
    }
    i++;
  }
  
  
#ifdef _PRINT  
  fprintf(log_file, "EXIT[get_points].\n");
#endif		  
  return 0;
}

/**************************************************************
FUNCTION:  get_rho
DESCRIPTION:  This returns particle frequency of given size. 

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: particle_frequency
***************************************************************/

double get_rho(double part_size, int *gl_ptr)
{
	double func_rho=0.0;
	double size1, size2;
	//double cum0;
	double cum1, cum2;
	int i=1;
	//printf("total num of grain file in line 536 in new_tephra is %d\n\n", *gl_ptr);
	//printf("line 537 max grain size and its content are = %1.1f\t%1.1f\n", g_ptr[*gl_ptr-1].sizeclass, g_ptr[*gl_ptr-1].content);
	/*
	size1 = g_ptr[i-1].sizeclass;
	size2 = g_ptr[i].sizeclass;
	cum1  = g_ptr[i-1].content;
	cum2  = g_ptr[i].content;
	*/
	
	

	size1 = g_ptr[i-1].sizeclass;
	size2 = g_ptr[i].sizeclass;
	cum1  = g_ptr[i-1].content;
	cum2  = g_ptr[i].content;
	
	//printf("line 170 %1.1f, %1.3e; size = %1.1f; %1.1f, %1.3e\n", size1, cum1, part_size, size2, cum2);
	
	while(size1 <= part_size)
	{
		/*printf("here I come; i= %d\n", i);*/
		if(part_size < size2)
		{
			func_rho = (cum2 - cum1) / ((size2 - size1) / 0.1);
			//printf("line 176 (%1.1f, %1.3e); size = %1.1f rho = %1.1e; (%1.1f, %1.3e)\n", size1, cum1, part_size, func_rho, size2, cum2);
			
			break;
		}
		else
		{
			//printf("else 557\n");
			i++;
			if(i == *gl_ptr){func_rho = 0.0; break;}
			size1 = g_ptr[i-1].sizeclass;
			size2 = g_ptr[i].sizeclass;
			cum1  = g_ptr[i-1].content;
			cum2  = g_ptr[i].content;	
		}
	}
	/*printf("here I come; func_rho= %1.4f\n", func_rho);*/
	return func_rho;
}
/**************************************************************
FUNCTION:  get_wind
DESCRIPTION:  This function reads wind data into the
WIND array. Each node stores all of the wind data. 

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/

int get_wind(FILE *in) {
  
  int i=0, j=0, ret;
  double	z, s, x, 
		dir, northing, easting, 
		ta, p, dens_atm, 
		dens_column, gas_cont, Q, 
		Cp, Rg, V, 
		Ue, M, theta, 
		U, R, T;
  
#ifdef _PRINT
  fprintf(log_file,"ENTER[get_wind].\n");
#endif
  
  //WIND_INTERVAL = (PLUME_HEIGHT - VENT_ELEVATION)/COL_STEPS;
	

	W = (WIND**)malloc(WIND_DAYS * sizeof(WIND *));
	
  if (W == NULL) {
    fprintf(stderr, "Cannot malloc memory for wind columns:[%s]\n", strerror(errno));
    return -1;
  } else {
	for (i=0; i < WIND_DAYS; i++) {

		W[i] = (WIND *)malloc((S_STEPS+1) * sizeof(WIND));

      if (W[i] == NULL) {
	fprintf(stderr, "Cannot malloc memory for wind rows %d:[%s]\n", i, strerror(errno));
	return -1;
      }
    }
  }
  /* Assume one wind day */
  i=0;
    
  /* Do for each column step */
  /* j = 0 is for the interval between the vent and the ground.
   * Here we set the wind speed and direction to be at the level of the vent;
   * The values used in the calculations change for each location and are 
   * set in the tephra_calc routine when the point elevation is known. 
   * The last interval ends at the top of the column. 
   */
	FILE *op;
	op=fopen("column.txt", "r");
	char line[200];
	
	while(NULL != fgets(line, MAX_LINE2, op)){
		
		if(line[0] != '#'){
			sscanf(line,
						"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
						&z, &s, &x, 
						&dir, &northing, &easting, 
						&ta, &p, &dens_atm, 
						&dens_column, &gas_cont, &Q, 
						&Cp, &Rg, &V, 
						&Ue, &M, &theta, 
						&U, &R, &T
						);
						
						
						W[i][j].wind_height = z;
						W[i][j].wind_speed = V;
						W[i][j].wind_dir = dir;
						W[i][j].t_atm = ta;
						W[i][j].p_atm = p;
						
			  	      fprintf(log_file, 
			  	      "%f\t%f\t%f\n", 
			  				W[i][j].wind_height, W[i][j].wind_speed, W[i][j].wind_dir);
							
						W[i][j].wind_dir *= DEG2RAD; /* change to radians */
						j++;
		}
	}	
	fclose(op);

  
#ifdef _PRINT
  fprintf(log_file, "\tRead %d wind days with %d wind levels per day.\n", i, j);
  fprintf(log_file, "EXIT[get_wind].\n");
#endif	  	  
  return 0;
}

double get_dir2(double level, double ht1, double ht0, double dir1, double dir0){
	double wind_dir;
	
	if(dir1 - dir0 > 180){
		dir0 = dir0 + 360; 
	}else if(dir0 - dir1 > 180){
		dir1 = dir1 + 360; 
	}
	
	wind_dir = ((dir1 - dir0) * (level - ht0) / (ht1 - ht0)) + dir0;
	
	if(wind_dir>360){wind_dir=wind_dir-360;}
	if(wind_dir<0){wind_dir=wind_dir+360;}
	
	return wind_dir;
}
	        //((wind_dir - dir0) * (level - ht0) / (wind_height - ht0)) + dir0;

/**************************************************************
FUNCTION:  get_wind2
DESCRIPTION:  This function reads wind data into the
W2, which are the conditions between vent height and the sea level

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error

Written by K. Mannen (18-Jan-2011) modified from get_wind()
***************************************************************/

int get_wind2(FILE *in_wind) {
  
	int i=0, j=0, ret;
	int jmax =0;
	char line[MAX_LINE];
	double wind_height, wind_dir, windspeed, dir0, ht0, sp0;
	double level;
  
#ifdef _PRINT
	fprintf(log_file,"ENTER[get_wind].\n");
#endif
	
	COL_STEPS = 100;
	/* print wind header in _node */
	fprintf(log_file, "height\tspeed\tdirection\n");
  
	WIND_INTERVAL = ds;
	/*printf("wind interval = %1.1f\n", WIND_INTERVAL);*/
  
	jmax = floor(VENT_ELEVATION / WIND_INTERVAL);
	/*printf("jmax = %d\n", jmax);*/

	W2 = (WIND**)malloc(WIND_DAYS * sizeof(WIND *));
	
	if (W2 == NULL)
	{
		fprintf(stderr, "Cannot malloc memory for wind columns:[%s]\n", strerror(errno));
		return -1;
	}else{
	 for (i=0; i < WIND_DAYS; i++)
		{
			W2[i] = (WIND *)malloc((jmax+1) * sizeof(WIND));
			
		if (W2[i] == NULL) {
			fprintf(stderr, "Cannot malloc memory for wind rows %d:[%s]\n", i, strerror(errno));
			return -1;
		}
	 }
	}
 
 /* Assume one wind day */
 i=0;
 /* start at sea level */ 
 level = VENT_ELEVATION - jmax * WIND_INTERVAL;
 /* Do for each column step */
 /* j = 0 is for the interval between the vent and the ground.
  * Here we set the wind speed and direction to be at the level of the vent;
  * The values used in the calculations change for each location and are 
  * set in the tephra_calc routine when the point elevation is known. 
  * The last interval ends at the top of the column. 
  */
  for (j=0; j <= jmax; j++) { 
    W2[i][j].wind_height = 0.0;
	ht0 = 0;
    dir0 = 0.0;
    sp0 = 0.0;
	 
	/*--------------------*/
	 
    /* Find wind elevation just greater than current level */
    /* Start scanning the wind file for the best match.
     * Each new level starts scanning the file from the beginning.
     */
							/*---case of error---*/
	while (NULL != fgets(line, MAX_LINE, in_wind))      /*a*/
	{
		if (line[0] == '#' || strlen(line) < WIND_COLUMNS) continue;
		else {								/*b*/
	      while (ret = sscanf(line,		/*c*/
			      "%lf %lf %lf",
			      &wind_height,
			      &windspeed,
			      &wind_dir), ret != 3) { 
	    
	        if (ret == EOF && errno == EINTR) continue;
	        
	        fprintf(stderr, 
	        "[line=%d,ret=%d] Did not read in 3 parameters:[%s]\n", 
	        i+1,ret, strerror(errno));
	        return -1;
	      } /*end of while c*/
	  } /*end of else b*/
    /* This is the case where we find the first height that is equal to
     * or greater that the level that we are assigning.
	*/
	
							/*--------------------*/
	
	if (wind_height >= level) {
			  if(wind_height == level) {
				  W2[i][j].wind_dir = wind_dir;
				  W2[i][j].wind_speed = windspeed;
	        
				} else { /* interpolate */
				  W2[i][j].wind_dir = 
					((wind_dir - dir0) * (level - ht0) / (wind_height - ht0)) + dir0;
	      
				  W2[i][j].wind_speed = 
					((windspeed - sp0) * (level - ht0) / (wind_height - ht0)) + sp0;
			     }      /*end of interpolate */
			W2[i][j].wind_height = level;
	      fprintf(log_file, 
	      "%f\t%f\t%f\n", 
				W2[i][j].wind_height, W2[i][j].wind_speed, W2[i][j].wind_dir);
		
	      W2[i][j].wind_dir *= DEG2RAD; /* change to radians */
	      break; /* ready to rescan the file for a match for the next level */
	}  /* end if */
	
	
							/*--------------------*/
	
	
	    /* This is the case where the scanned height is less than the level
	     * we are assigning.
	     */
	    else{  /* else case of line 140  // in other words, wind_height < level */
	      /* Maintain the scanned values for possible interpolation 
	       * at the next level.
	       */
	      ht0 = wind_height;
	      dir0 = wind_dir;
	      sp0 = windspeed;
	    }  /* end of else */
	
	
    } /* end of while loop*/

	  /* If we finish scanning the file and all heights are below the level we are
	   * currently assigning, then just use the direction and speed
	   * at the upper-most height.
	   */
	  if (!W2[i][j].wind_height) 
	  {
	    W2[i][j].wind_height = level;
	    W2[i][j].wind_speed = sp0;
	    W2[i][j].wind_dir = dir0;
	  }
	  /* Go to the next column height */
	  		  /*printf("line 204 %d %f %f %f\n", 
				j, W2[i][j].wind_height, W2[i][j].wind_speed, W2[i][j].wind_dir); */
	  rewind(in_wind);
	  level += WIND_INTERVAL;  

  }/*end of for loop*/
		/*printf("line 210 j = %d\n", i);*/
		return 0;
}

/* end of the function get_wind2*/



/**************************************************************
FUNCTION:  get_config_data
DESCRIPTION:  This function reads the configuration file,
and sets some global variables.

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/

int init_globals(char *config_file) {

  FILE *in_config;
  char buf[1][30], **ptr1;
  char line[MAX_LINE];
  char space[4] = "\n\t ";
  char *token;
  
#ifdef _PRINT  
  fprintf(log_file, "ENTER[init_globals].\n");
#endif
  
  in_config = fopen(config_file, "r");
  if (in_config == NULL) {
    fprintf(stderr, 
	    "Cannot open configuration file=[%s]:[%s]. Exiting.\n", config_file, strerror(errno)); 
    return 1;
  }
  
  ptr1 = (char **)&buf[0];
  while (fgets(line, MAX_LINE, in_config) != NULL) {
    /*fprintf(stderr, "%s\n", line); */
    if (line[0] == '#' || line[0] == '\n') continue;
    
    token = strtok_r(line, space, ptr1);
    if (!strncmp(token, "DIFFUSION_COEFFICIENT", strlen("DIFFUSION_COEFFICIENT"))) {
      token = strtok_r(NULL,space,ptr1);
      DIFFUSION_COEFFICIENT = strtod(token, NULL);
      fprintf(stderr, "DIFFUSION_COEFFICIENT=%.1f\n", DIFFUSION_COEFFICIENT);
	}
	else if (!strncmp(token, "EDDY_CONST", strlen("EDDY_CONST"))) {
      token = strtok_r(NULL,space,ptr1);
      EDDY_CONST = strtod(token, NULL);
      fprintf(stderr, "EDDY_CONST=%g\n", EDDY_CONST);
    }
    else if (!strncmp(token, "FALL_TIME_THRESHOLD", strlen("FALL_TIME_THRESHOLD"))) {
      token = strtok_r(NULL,space,ptr1);
      FALL_TIME_THRESHOLD = strtod(token, NULL);
      fprintf(stderr, "FALL_TIME_THRESHOLD=%.1f\n", FALL_TIME_THRESHOLD);
    }
    else if (!strncmp(token, "LITHIC_DENSITY", strlen("LITHIC_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      LITHIC_DENSITY = strtod(token, NULL);
      fprintf(stderr, "LITHIC_DENSITY=%.1f\n", LITHIC_DENSITY);
    }
    else if (!strncmp(token, "PUMICE_DENSITY", strlen("PUMICE_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      PUMICE_DENSITY = strtod(token, NULL);
      fprintf(stderr, "PUMICE_DENSITY=%.1f\n", PUMICE_DENSITY);
    }
    else if (!strncmp(token, "PART_STEPS", strlen("PART_STEPS"))) {
      token = strtok_r(NULL, space, ptr1);
      PART_STEPS = (int)atoi(token);
      fprintf(stderr, "PART_STEPS = %d\n", PART_STEPS);
    }
    else if (!strncmp(token, "PLUME_MODEL", strlen("PLUME_MODEL"))) {
      token = strtok_r(NULL,space,ptr1);
      PLUME_MODEL = (int)atoi(token);
      
      if (PLUME_MODEL == 0) {
	      pdf = plume_pdf0;
	      fprintf(stderr, "PLUME_MODEL=[%d]%s\n", PLUME_MODEL, "Uniform Distribution with threshold");
      }
      else if (PLUME_MODEL == 1) {
	      pdf = plume_pdf1;
	      fprintf(stderr, "PLUME_MODEL=[%d]%s\n", PLUME_MODEL, "log-normal Distribution using beta");
      }
    }
    else if (!strncmp(token, "PLUME_RATIO", strlen("PLUME_RATIO"))) {
      token = strtok_r(NULL, space, ptr1);
      PLUME_RATIO = strtod(token, NULL);
      if (!PLUME_MODEL) fprintf(stderr, "PLUME_RATIO = %.2f\n", PLUME_RATIO);
    }
	else if (!strncmp(token, "SIGMA_PLUME", strlen("SIGMA_PLUME"))) {          /* added by Kaz 02-Dec-2010 */
      token = strtok_r(NULL,space,ptr1);
      SIGMA_PLUME = strtod(token, NULL);
      fprintf(stderr, "SIGMA_PLUME=%g\n", SIGMA_PLUME);
	}
	else if (!strncmp(token, "SUZUKI_A", strlen("SUZUKI_A"))) {                /* added by Kaz 11-Jan-2011 */
      token = strtok_r(NULL,space,ptr1);
      SUZUKI_A = strtod(token, NULL);
      fprintf(stderr, "SUZUKI_A=%g\n", SUZUKI_A);
	}
	else if (!strncmp(token, "SUZUKI_LAMDA", strlen("SUZUKI_LAMDA"))) {        /* added by Kaz 11-Jan-2011 */
      token = strtok_r(NULL,space,ptr1);
      SUZUKI_LAMDA = strtod(token, NULL);
      fprintf(stderr, "SUZUKI_LAMDA=%g\n", SUZUKI_LAMDA);
	}
	else if (!strncmp(token, "WIND_DAYS", strlen("WIND_DAYS"))) {
      token = strtok_r(NULL, space, ptr1);
      WIND_DAYS = (int)atoi(token);
      fprintf(stderr, "WIND_DAYS = %d\n", WIND_DAYS);
    }
    else if (!strncmp(token, "WIND_COLUMNS", strlen("WIND_COLUMNS"))) {
      token = strtok_r(NULL, space, ptr1);
      WIND_COLUMNS = (int)atoi(token);
      fprintf(stderr, "WIND_COLUMNS = %d\n", WIND_COLUMNS);
    }
    else if (!strncmp(token, "PLUME_THICKNESS", strlen("PLUME_THICKNESS"))) {        /* added by Kaz 09-Mar-2020 */
      token = strtok_r(NULL, space, ptr1);
      PLUME_THICKNESS = strtod(token, NULL);
      fprintf(stderr, "PLUME_THICKNESS = %.1f\n", PLUME_THICKNESS);
    }
    else if (!strncmp(token, "SOURCE_DECAY_RATE", strlen("SOURCE_DECAY_RATE"))) {        /* added by Kaz 16-SEP-2021 */
      token = strtok_r(NULL, space, ptr1);
      SOURCE_DECAY_RATE = strtod(token, NULL);
      fprintf(stderr, "SOURCE_DECAY_RATE = %g\n", SOURCE_DECAY_RATE);
    }
    else if (!strncmp(token, "KS_ENTRAIN_U", strlen("KS_ENTRAIN_U"))) {       /* added by Kaz 27-OCT-2021 */
      token = strtok_r(NULL, space, ptr1);
      KS_ENTRAIN_U = strtod(token, NULL);
      fprintf(stderr, "KS_ENTRAIN_U = %g\n", KS_ENTRAIN_U);
    }
    else if (!strncmp(token, "KW_ENTRAIN_V", strlen("KW_ENTRAIN_V"))) {        /* added by Kaz 27-OCT-2021 */
      token = strtok_r(NULL, space, ptr1);
      KW_ENTRAIN_V = strtod(token, NULL);
      fprintf(stderr, "KW_ENTRAIN_V = %g\n", KW_ENTRAIN_V);
    }
    else if (!strncmp(token, "ERUPTION_MASS", strlen("ERUPTION_MASS"))) {
      token = strtok_r(NULL, space, ptr1);
      ERUPTION_MASS = strtod(token, NULL);
      fprintf(stderr, "ERUPTION_MASS = %g\n", ERUPTION_MASS);
    }
    else if (!strncmp(token, "MAX_GRAINSIZE", strlen("MAX_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MAX_GRAINSIZE = strtod(token, NULL);
      fprintf(stderr, "MAX_GRAINSIZE = %.0f\n", MAX_GRAINSIZE);
    }
    else if (!strncmp(token, "MIN_GRAINSIZE", strlen("MIN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MIN_GRAINSIZE = strtod(token, NULL);
      fprintf(stderr, "MIN_GRAINSIZE = %.0f\n", MIN_GRAINSIZE);
    }
    else if (!strncmp(token, "MEDIAN_GRAINSIZE", strlen("MEDIAN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      MEDIAN_GRAINSIZE = strtod(token, NULL);
      fprintf(stderr, "MEDIAN_GRAINSIZE = %.2f\n", MEDIAN_GRAINSIZE);
    }
    else if (!strncmp(token, "STD_GRAINSIZE", strlen("STD_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      STD_GRAINSIZE = strtod(token, NULL);
      fprintf(stderr, "STD_GRAINSIZE = %.2f\n", STD_GRAINSIZE);
    }
    else if (!strncmp(token, "VENT_EASTING", strlen("VENT_EASTING"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_EASTING = strtod(token, NULL);
      fprintf(stderr, "VENT_EASTING = %.1f\n", VENT_EASTING);
    }
    else if (!strncmp(token, "VENT_NORTHING", strlen("VENT_NORTHING"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_NORTHING = strtod(token, NULL);
      fprintf(stderr, "VENT_NORTHING = %.1f\n", VENT_NORTHING);
    }
    else if (!strncmp(token, "VENT_ELEVATION", strlen("VENT_ELEVATION"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_ELEVATION = strtod(token, NULL);
      fprintf(stderr, "VENT_ELEVATION = %.1f\n", VENT_ELEVATION);
    }
    else if (!strncmp(token, "INITIAL_WATER_CONTENT", strlen("INITIAL_WATER_CONTENT"))) {
      token = strtok_r(NULL, space, ptr1);
      INITIAL_WATER_CONTENT = strtod(token, NULL);
      fprintf(stderr, "INITIAL_WATER_CONTENT = %.1f\n", INITIAL_WATER_CONTENT);
    }
    else if (!strncmp(token, "MAGMA_DISCHARGE_RATE", strlen("MAGMA_DISCHARGE_RATE"))) {
      token = strtok_r(NULL, space, ptr1);
      MAGMA_DISCHARGE_RATE = strtod(token, NULL);
      fprintf(stderr, "MAGMA_DISCHARGE_RATE = %.1f\n", MAGMA_DISCHARGE_RATE);
    }
    else if (!strncmp(token, "MAGMA_TEMPERATURE", strlen("MAGMA_TEMPERATURE"))) {
      token = strtok_r(NULL, space, ptr1);
      MAGMA_TEMPERATURE = strtod(token, NULL);
      fprintf(stderr, "MAGMA_TEMPERATURE = %.1f\n", MAGMA_TEMPERATURE);
    }
    else if (!strncmp(token, "INITIAL_PLUME_VELOCITY", strlen("INITIAL_PLUME_VELOCITY"))) {
      token = strtok_r(NULL, space, ptr1);
      INITIAL_PLUME_VELOCITY = strtod(token, NULL);
      fprintf(stderr, "INITIAL_PLUME_VELOCITY = %.1f\n", INITIAL_PLUME_VELOCITY);
    }
    else if (!strncmp(token, "VENT_RADIUS", strlen("VENT_RADIUS"))) {
      token = strtok_r(NULL, space, ptr1);
      VENT_RADIUS = strtod(token, NULL);
      fprintf(stderr, "VENT_RADIUS = %.1f\n", VENT_RADIUS);
    }
    else if (!strncmp(token, "S_MAX", strlen("S_MAX"))) {
      token = strtok_r(NULL, space, ptr1);
      S_MAX = strtod(token, NULL);
      fprintf(stderr, "S_MAX = %.1f\n", S_MAX);
    }
    else continue;
  }
  (void) fclose(in_config);

  
#ifdef _PRINT
  fprintf(log_file, "EXIT[init_globals].\n");
#endif

  return 0;
}

void free_memory_new_calc(void){

}
