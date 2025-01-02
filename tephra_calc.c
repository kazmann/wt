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

static double threshold;
static double AIR_VISCOSITY_x_18;
static double LITHIC_DENSITY_minus_PUMICE_DENSITY;
static double PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD;
static double ONE_THIRD;
static double AIR_VISCOSITY_x_225;
static double GRAV_SQRD_x_4;
static double BETA_x_SQRT_TWO_PI;
static double TWO_BETA_SQRD;
static double PDF_GRAINSIZE_DEMON1;
static double TWO_x_PART_SIGMA_SIZE;
static double EDDY_CONST_x_8_div_5;
static TABLE **T;
static COLUMN *C;
static FILE *log_file;

/*
Code: tephra_calc.c
By: C.B. and L.J. Connor & T. Hincks & C. Bonadonna
Copyright (C) 2003  C.B. Connor, L.J. Connor, C. Bonadonna, and T. Hincks
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
/* ---------------------------------------------------------------------
 * FUNCTION:tephra_calc.c 
 *
 * Purpose: This function calculates and returns the expected accumulation
 * of volcanic ash (kg/m2) at a specific geogrphic location (x,y)
 * due to an eruption with specific input parameters. 
 * These points may be random or on  a UTM grid (m)
 * 
 * This implementation accounts for variation in wind velocity with height.
 * The model is discretized w.r.t. height and particle size. 

 * This function is called for each point (x,y,) If more than one eruption is
 * involved, for example in a probabilistic analysis, the function is called for
 * each set of eruption parameters. 
 
 * INPUTS:
 * ERUPTION *erupt: pointer to array of eruption parameters
 * POINT *pt: pointer to an array of location specific parameters, 
 * WIND *level: pointer to a day of wind data :
              height asl in m; 
              wind speed in ms-1
              wind direction in degrees N
 *	     
 * OUTPUTs:
 *   the value of the mass accumulated at the input location (northing, easting) in kg/m2 
 *
 *   a distribution of particle sizes
 *   the exact number of binss (i.e. sizes) and phi size used per bin is an integer and 
 *   is determined by (erupt->max_phi - erupt->min_phi)
 *   each bin accumulates phi sizes up to its integer size
 *   ex. bin[0] holds grainsizes [min_phi to min_phi+1)  
 ***************************************************************************/

void tephra_calc(ERUPTION *erupt, POINT *pt, WIND *level2, STATS *stats) { /* tephra_calc starts ... */
 
   /**********************************************************************************
   * WIND structure:
   * level->windspeed: windspeed in m/s
   * level->wind_dir: wind direction in +/- degrees from north
   * level->wind_height: meters above sea level
   *
   * See common_structures_strat.h for structure definitions.
   **********************************************************************************/
  
  
	int i = 0, j = 0, bin = -1; 
	double new_xspace, new_yspace, ht_above_vent; /*cos_wind = 0.0, sin_wind = 0.0, windspeed = 0.0*/
	double sigma, demon2, demon3, ash_fall, layer, fall_time_adj = 0.0, total_fall_time=0.0;
	double average_windspeed_x, average_windspeed_y, average_wind_direction, average_windspeed =0.0;
	double wind_sum_x = 0.0, wind_sum_y = 0.0;
	double pnt_vent_x, pnt_vent_y; // Added by K. Mannen (Feb 11, 2017)
	
  /* these are added by K. Mannen (19-Jan-2011) */
	double particle_ht= erupt->vent_height;
	double part_density;
	double ground_alt=0;
	double ashdiam;
  /* ------------------------------------------ */
	
	static double min=10e6, max=0.0;
	
	/* array and its pointer to retrieve wind condition
	between vent height and the sea level
	*/
	double wind_element[2];
	double *ptr_wind;        /* arry used as data passing between part_fall_time_vg() */
	                         /* wind advection of x and y coordinate will be stored   */
	
	double temp0, temp1, temp2, xprime, yprime, demon1;
	
	/* they are added by K. Mannen (19-Jan-2011) */
	
	
  
#ifdef _PRINT
	fprintf(log_file, "IN tephra_calc ...");
#endif

  /* Initialize mass to zero */
	pt->mass = 0.0;
	wind_sum_x = 0.0;
	wind_sum_y = 0.0;
	
  /* Transform the volcano location coordinate to 0,0 
   */
	pnt_vent_x = pt->northing - erupt->volcano_northing;	// Edited by K. Mannen (Feb 11, 2017)
	pnt_vent_y = pt->easting - erupt->volcano_easting;		// Edited by K. Mannen (Feb 11, 2017)

  /* do the double integration over grainsize and column height 
   */
#ifdef _PRINT
	fprintf(log_file, "\nBeginning integration loops ... \n");
#endif
  
  /* Interpolate to fine the windspeed and direction below the height of the vent.
   * Find one average speed and direction between vent and grid elevation point. 
   * The first values in the wind array give the wind speed and direction at the
   * vent height.
   */
	layer = erupt->vent_height - pt->elevation;
	ground_alt = pt->elevation;
	/*windspeed = (level[0].wind_speed * pt->elevation) / erupt->vent_height;*/
	/*cos_wind = cos(level[0].wind_dir);*/ /* windspeed;*/
	/*sin_wind = sin(level[0].wind_dir);*/ /* windspeed;*/
	   
	for (i = 0; i < PART_STEPS; i++) { /* PART_STEPS_LOOP */
   
		fall_time_adj = 0.0;
		part_density = T[i][0].part_density;
		ashdiam = T[i][0].ashdiam;




		/*printf("153 dens = %1.4f diam = %1.4f\n", part_density, ashdiam);*/
		
		/* Accumulate the particle sizes into bins of whole numbered phi sizes */
		if (!(i % 10)) { 
			bin++;
#ifdef _PRINT
			fprintf(log_file, "PART_STEP=%d phi[%d] = %g\n", i, bin, pt->phi[bin]);
#endif
		}
  
		/* Adjust the total fall time of each particle size (i) 
		   by the time it takes to descend from vent height to the grid cell (pt) 
		*/           
		if (layer > 0) {
			fall_time_adj = 
			part_fall_time_vg(particle_ht, ground_alt, level2, ptr_wind, ashdiam, part_density);
			
			/*printf("particle cos pointer = %1.2f\n", wind_element[0]);
			  printf("particle sin time pointer = %1.2f\n", wind_element[1]);*/
			
			
#ifdef DEBUG 	
			fprintf(log_file, "%d %g %g\n",  i, layer, fall_time_adj) ;
#endif	 
		}


#ifdef _PRINT
	//fprintf(log_file, "LINE280\tashdiam\tparticle_ht\tv-pntx\tv-pnty\tcolsourcex\tcolsourcey\tcolheight\tnew_xspace\tnew_yspace");
	//fprintf(log_file, "\ttemp0\ttemp1\txprime\typrime\ttemp2\tdemon1\twind_sum_x\twind_sum_y");
	//fprintf(log_file, "\tave_wind_sp\tave_wind_dir\tSourceMag\tdemon2\tdemon3\tsigma\ttotal_fall_time\tmassloading\n");
#endif
	
		for (j = 0; j < S_STEPS; j++) { /* COL_STEPS_LOOP */
     
  		total_fall_time = T[i][j].total_fall_time + fall_time_adj;
		//printf("line 186 falltime plume-VE = %1.3f VE-GL = %1.3f plume-GL = %1.3f\n", T[i][j].total_fall_time, fall_time_adj, total_fall_time);
    	// fprintf(stderr, "%g %g %g ", T[i][j].total_fall_time, total_fall_time, fall_time_adj);

    	/* Sum the adjustments (windspeed and wind_direction) 
    	 * for each particle size  falling from each level.
		 */
			
   		wind_sum_x = ptr_wind[0]; /*modified by K. Mannen (19-Jan-2011) */
	  	wind_sum_y = ptr_wind[1]; /*modified by K. Mannen (19-Jan-2011) */
		//printf("line 196 wind_x = %1.3f wind_y = % 1.3f\n", wind_sum_x, wind_sum_y);
		
		
	    
	    /* Now add the summed adjustments to the already summed
	     * windspeeds and directions 
	     and
    	 Account for the wind:
    	 Find the average windspeed in the x and y directions 
    	 over the total fall time.
			*/
			average_windspeed_x = 
			(T[i][j].wind_sum_x + wind_sum_x)/total_fall_time;
			
			average_windspeed_y = 
			(T[i][j].wind_sum_y + wind_sum_y)/total_fall_time;
    	    	
			/* If zero, make windspeed a very small value (cannot divide by zero in next step) */
			if (!average_windspeed_x) average_windspeed_x = .001;
			if (!average_windspeed_y) average_windspeed_y = .001;
      
			/* Find the average wind direction (direction of the velocity vector) */
			if (average_windspeed_x < 0) {
				average_wind_direction = 
				atan(average_windspeed_y/average_windspeed_x ) + M_PI;
			} else 
    	  average_wind_direction = 
    	  atan(average_windspeed_y/average_windspeed_x);
    	
			/* Find the average windspeed ( magnitude of the velocity vector) */
			average_windspeed = 
			sqrt(average_windspeed_x*average_windspeed_x + average_windspeed_y*average_windspeed_y);
				
			if (total_fall_time > max) max = total_fall_time;
			if (total_fall_time < min) min = total_fall_time;
				
			/* calculate the value of sigma (dispersion) based on total_fall_time  
			 * to acct for the change in the shape of the column with ht - increasing radius 
			 */
			ht_above_vent = T[i][j].particle_ht - erupt->vent_height;
			
			/* falltime for fine particles */
			if (total_fall_time >= FALL_TIME_THRESHOLD) {
				sigma = 
				EDDY_CONST_x_8_div_5 * pow((total_fall_time + T[i][j].plume_diffusion_fine_particle), 2.5);
				//fprintf(stderr,"f");
			} else { /* falltime for coarse particles */
				sigma = 
				4.0 * DIFFUSION_COEFFICIENT * (total_fall_time + T[i][j].plume_diffusion_coarse_particle);
				//fprintf(stderr, "c");
			}

			demon2 =  M_PI * sigma;
			/*demon2 =  sqrt(demon2);*/
      
    	/* Modify fall time by the variation of wind velocity with height */
			new_xspace = pnt_vent_x - C[j+1].centre_x;
			new_yspace = pnt_vent_y - C[j+1].centre_y;
			//printf("line 254 average_wind direction = %1.3f speed = % 1.3f\n", average_wind_direction, average_windspeed);
			demon3 = 
			strat_average( average_wind_direction, 
      	            	 	 average_windspeed,             
				             new_xspace, new_yspace, 
				             total_fall_time,
				             sigma); 
			//printf("line261 demon1 = %1.4e, 2 = %1.4e, 3 = %1.4e\n", T[i][j].demon1, demon2, demon3);

			//if(i == 0 && j == 0){printf("line263\ti\tj\tdemon1\tdemon2\tdemon3\tAWD\tAWS\tX\tY\tTFT\tSig\n");}
			//printf("line264\t%d\t%d\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%1.4e\n", i, j, T[i][j].demon1, demon2, demon3, average_wind_direction, average_windspeed, new_xspace, new_yspace, total_fall_time, sigma);
/*
			if (!demon2 || isnan(demon2) || isinf(demon2) || isnan(demon3) || isinf(demon3)) {
 				fprintf(stderr, 
      	"[%d][%d] layer= %.1f totalfalltime=%g [falltimeadj=%g] demon1=%g demon2=%g demon3=%g sigma=%g\n",
      	i,j, layer,total_fall_time, fall_time_adj, T[i][j].demon1, demon2, demon3, sigma);
      	exit(-1);
			}
 */    
							 
			ash_fall = (T[i][j].demon1 / demon2) * demon3;
			pt->mass += ash_fall;
			pt->phi[bin] += ash_fall;
			
			//double temp0, temp1, temp2, xprime, yprime, demon1;
		    temp0 = cos(average_wind_direction);
		    temp1 = sin(average_wind_direction);
    
		    xprime = new_xspace * temp0 + new_yspace * temp1;
		    yprime = new_yspace * temp0 - new_xspace * temp1;
    
		    temp2 = xprime - average_windspeed * total_fall_time;
		    demon1 = temp2 * temp2 + yprime * yprime;
			
#ifdef _PRINT
	//fprintf(log_file, "LINE280\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", ashdiam, T[i][j].particle_ht, pnt_vent_x, pnt_vent_y, C[j+1].centre_x, C[j+1].centre_y, C[j+1].z, new_xspace, new_yspace);
	//fprintf(log_file, "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", temp0, temp1, xprime, yprime, temp2, demon1, T[i][j].wind_sum_x, T[i][j].wind_sum_y);
	//fprintf(log_file, "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", average_windspeed, average_wind_direction, T[i][j].demon1, demon2, demon3, sigma, total_fall_time, ash_fall);
#endif
	
		}   // end of column step (j)
		//if(i == 0){printf("line279\ti\tbin\tphi\n");}
		//printf("line280\t%d\t%d\t%g\n", i, bin, pt->phi[bin]);
	}       // end of particle step (i)
#ifdef _PRINT
  fprintf(log_file, "PART_STEP=%d phi[%d] = %g\n", i, bin, pt->phi[bin]);
  fprintf(log_file, "OUT\n");
#endif
  stats->min_falltime = min;
  stats->max_falltime = max;
}


/* ----------------- New Function Starts Here -------------------- */
/* Function strat_average accounts for the variation in wind velocity 
   with height by using the average velocity value
*/
double strat_average( double average_wind_direction, 
                      double average_windspeed,             
			                double xspace, double yspace, 
			                double total_fall_time,
			                double sigma) {
		
		double temp0, temp1, xprime, yprime, demon1, demon3;
			                  
    temp0 = cos(average_wind_direction);
    temp1 = sin(average_wind_direction);
    
    xprime = xspace * temp0 + yspace * temp1;
    yprime = yspace * temp0 - xspace * temp1;
    
    temp0 = xprime - average_windspeed * total_fall_time;
    demon1 = temp0 * temp0 + yprime * yprime;
    demon3 = exp(-demon1/sigma); /* where sigma is calculated for the total fall time */
    return demon3;
			          
}

/* ----------------- New Function Starts Here -------------------- */
/* function phi2m converts the ash diameter from 
   units of phi to m
*/

 double phi2m(double xx) {
   double cms;
   //printf("in phi2m");
   cms = 0.001 * pow(2, -xx);
   return cms;
 }
/* ----------------- New Function Starts Here -------------------- */
/* function particle_density calculates varying particle density based on their grain size diamete
   using a linear correlation between pumice_threshold (PHI) and lithic_threshold (PHI)
*/

 double particle_density (double phi_slice) {

  double mean_density = 0;

  if (phi_slice >= LITHIC_DIAMETER_THRESHOLD) mean_density = LITHIC_DENSITY;
  else if (phi_slice <= PUMICE_DIAMETER_THRESHOLD) mean_density = PUMICE_DENSITY;
  else if (phi_slice < LITHIC_DIAMETER_THRESHOLD && phi_slice > PUMICE_DIAMETER_THRESHOLD) 
    mean_density = 
    
      LITHIC_DENSITY - 
      LITHIC_DENSITY_minus_PUMICE_DENSITY * 
      (phi_slice - LITHIC_DIAMETER_THRESHOLD) / PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD;

   return mean_density;
 }



/* ----------------- New Function Starts Here -------------------- */
/* function part_fall_time determines the time of particle fall within each falling step
   falling steps are here:
   set = particle rising steps = ht_step_width
   
   returns the particle fall time within each falling step

This function follows the approach outlined in Bonadonna et al. (1998) 
Briefly, particle fall time is calculated based on terminal velocities in
layers that are 1000 m thick. The terminal velocity is a function of the
particle Reynolds number, which varies with grainsize, air properties.

The thickness of the first layer (closest to the ground) is equal to the vent 
height. The vent_height is in meters above sea level. The area the ash falls on
is considered to be at sea level. This leads to some assumptions (!) near the
volcano...
*/


double part_fall_time(double particle_ht, double layer, double ashdiam, double part_density) {
  
	double rho, hz, temp0, temp1;
   	double vtl, vti, vtt;
   	double reynolds_number;
   	double particle_term_vel;
   	double particle_fall_time;
 
  	particle_fall_time = 0.0;
  	hz = particle_ht;  /* height of the particle above sea level */
    
  	/*rho is the density of air (kg/m^3) at the elevation of the current particle*/
  	temp0 = -hz / 8200.0;
  	rho = AIR_DENSITY * exp(temp0);
  
	/*
   	(friction due to the air) :
    	vtl is terminal velocity (m/s) in laminar regime RE<6 
    	vti is terminal velocity (m/s) in intermediate regime 6<RE<500
    	vtt is terminal velocity (m/s) in turbulent regime RE>500
  	*/
  	vtl = (part_density - rho) * GRAVITY * ashdiam * ashdiam / AIR_VISCOSITY_x_18; /* 18.0 * AIR_VISCOSITY */
  
  	/*
    	vti = ashdiam * 
    	pow(((4.0*GRAVITY*GRAVITY*erupt->part_mean_density *erupt->part_mean_density )/		(225.0*AIR_VISCOSITY*rho)),(1.0/3.0));
    	vtt=sqrt(3.1*erupt->part_mean_density *GRAVITY*ashdiam/rho);
  	*/
  
  	/*
    	RE is calculated using vtl (RE is Reynolds Number)
  	*/
  	reynolds_number = ashdiam * rho * vtl / AIR_VISCOSITY;
  	particle_term_vel = vtl;
  	temp0 = ashdiam * rho;


  	/*
    	c...if laminar RE>6 (intermediate regime), RE is calculated again considering vti
  	*/
  
	if (reynolds_number >= 6.0) {

    		/*4.0 * GRAVITY * GRAVITY * part_density * part_density / AIR_VISCOSITY * 225.0 * rho */
    		temp1 = GRAV_SQRD_x_4 * (part_density - rho) * (part_density - rho) / (AIR_VISCOSITY_x_225 * rho); 
    		vti = ashdiam * pow(temp1, ONE_THIRD); /* ONE_THIRD = 1.0/3.0 */    
    		reynolds_number = temp0 * vti / AIR_VISCOSITY;
    		particle_term_vel = vti;
    		/*
    		c...if intermediate RE>500 (turbulent regime), RE is calculated again considering vtt 
  		*/ 
  		if (reynolds_number >= 500.0) {
    			vtt = sqrt( 3.1 * (part_density - rho) * GRAVITY * ashdiam / rho);
    			reynolds_number =  temp0 * vtt / AIR_VISCOSITY; 
    			particle_term_vel = vtt;
  		}  
  	}
/* Calculate the time it takes this particle to fall through this distance=layer */
  particle_fall_time = layer / particle_term_vel;
  
/* particle fall time is in sec   */
  
  //printf("i= %d, layer = %f, hz = %f, particle_term_vel = %f, diam=%f, reynolds = %f\n", i,layer_thickness, hz, particle_term_vel, a//shdiam, reynolds_number);
  
  return particle_fall_time;
}



/* ----------------- New Function Starts Here -------------------- */
/* function part_fall_time_vg determines the time of particle fall
  within each falling step between vent height and sea level


  this function utilizes part_fall_time()
  written by K. Mannen (18-Jan-2011)

  output is the fall time between vent elevaion and the sea level

  also average sin_wind and cos_wind are also dispatched using pointer *ptr_wind
  */

double part_fall_time_vg(double vent_elevation, double ground_elevation, WIND *level, double *ptr_wind, double ashdiam, double part_density)
{
	int j=0;
	double particle_ht, layer;			/* particle height, fall distance */
	double h0;							/* height particle fall to */
	double fall_time_ttl_vg, fall_time;	/* particle fall times between vent alt to ground, and through the slice */
	double cos_wind, sin_wind;                        /* wind advection in each slice*/
	double cos_sum, sin_sum;                 /* wind advection within vent-ground */

	
	fall_time = 0;
	fall_time_ttl_vg = 0;	// if ground_elevation > vent_elevation, fall time vg = 0
	h0 = level[j].wind_height;
	cos_wind=0; sin_wind=0; cos_sum=0; sin_sum=0;
	
	if(ground_elevation < vent_elevation) /* ground must be lower than vent */
	{
		while(ground_elevation >= level[j].wind_height)
		{
			h0 = level[j].wind_height;
			//printf("line 436 j= %d, h0 = %1.1f ground = %1.1f\n",j, h0, ground_elevation);
			j = j + 1;
		}
		
		particle_ht = level[j].wind_height;
		h0=ground_elevation;
		
	
		//printf("line 444 j= %d, particle_ht = %1.1f ground = %1.1f\n",j, particle_ht, h0);
	
		while(particle_ht <= vent_elevation)
		{
			layer = particle_ht - h0;
			if(layer<=0)break;
			
			
			//printf("452 particle_ht = %1.3f, wind_level = %1.3f, j = %d, diam = %1.4f dens = %1.4f\n", particle_ht, level[j].wind_height, j, ashdiam, part_density);
			fall_time = part_fall_time(particle_ht, layer, ashdiam, part_density);
			fall_time_ttl_vg += fall_time;
			
			//printf("456 j = %d, fall time = %1.4f ttl fall time vg = %1.4f\n", j, fall_time, fall_time_ttl_vg);
		
			/*  calc wind advection */
			cos_wind = cos(level[j].wind_dir) * level[j].wind_speed * fall_time;
			sin_wind = sin(level[j].wind_dir) * level[j].wind_speed * fall_time;
		
			cos_sum += cos_wind; 
			sin_sum += sin_wind;
		
			j += 1;
		
			h0 = particle_ht;
			particle_ht = level[j].wind_height;
		}
	}
	/*printf("ttl fall time = %1.1f, cos_wind = %1.1f, sin_wind = %1.1f \n", fall_time_ttl_vg, cos_sum, sin_sum);*/
	*ptr_wind =cos_sum; ++ptr_wind;
	*ptr_wind =sin_sum;
	return fall_time_ttl_vg;
}


/* ----------------- New Function Starts Here -------------------- */
/* this function calculates the expected fraction of particles
   in a given grainsize class (part_size_slice) assuming a normal 
   distribution in phi units about the mean, dmean,
   with standard deviation sigma.
   The probability that 

   modified by K. Mannen (19-Jan-2011) to use arbitrary size distribution 
*/

double pdf_grainsize(double part_mean_size, double part_size_slice, double part_step_width, int *gl_ptr) 
{
	double func_rho, temp;
	double demon3, demon2;
	
	func_rho=0;
	
	
	if(GRAIN_DATA==0)
	{
		
		/* PDF_sizeclass_DEMON1 = 1.0 / 2.506628 * erupt->part_sigma_size */
		//demon3   = part_size_slice - part_mean_size;
		demon3   = part_size_slice - part_step_width / 2 - part_mean_size; //20221226 modified
		temp = -demon3 * demon3 / TWO_x_PART_SIGMA_SIZE; /* 2.0 * erupt->part_sigma_size * erupt->part_sigma_size */
		demon2   = exp(temp);
		func_rho = PDF_GRAINSIZE_DEMON1 * demon2 * part_step_width; 
		if (func_rho < 0.0)
		{
			fprintf(log_file, "error in ash size distribution - method pdf_sizeclass");
		}
	}else if(GRAIN_DATA==1)
	{
		func_rho = get_rho(part_size_slice, gl_ptr);
	}
	else{/* GRAIN_DATA is not 0 nor 1 so needs to exit with error msg */}
	
	return func_rho;
}


/* 
   inputs:
   x: height of a particle within the plume, relative to vent height
   slice: integration step (index)
   ht_section_width: the width of an integration step
   none: not used
   
   output: the probability that a given grainsize will be released from a given height
*/

double plume_pdf0(double x, int slice, double none0, double none1) {
	/* In wt, 3rd and 4th arguments (none0 and none1) are not used */

	double probability;
	double fallout_threshold;
	static int num_slices_left = 0;
	static double plume_slice = 0.0;
  
  /* if (!slice) fprintf(stderr, "ENTER plume_pdf0 ....\n"); */
 
	probability = 0.0;
	fallout_threshold = threshold;
  if (x > fallout_threshold) {
     
    if (!num_slices_left) { // Since num_slices_left is "static int", it is not 0 after the first call
      num_slices_left = S_STEPS - slice;
      plume_slice = 1.0 / (double)num_slices_left;
      /* fprintf(stderr, "slices left = %d\n ", num_slices_left); */
    }

    probability = plume_slice;
  }
  
  //printf("x=%g threshold=%g plume_slice=%g prob=%g\n", x, threshold, plume_slice, probability); 
  //fprintf(stderr, "x=%g threshold=%g plume_slice=%g prob=%g\n", x, threshold, plume_slice, probability); 
  /*if (probability < 0.0) This only gets printed if an error occurs. 
    fprintf(stderr, "col_ht=%f prob=%f\n", x, probability); */
  return probability; 
}

double plume_pdf1(double z, int slice, double plume_height, double vent_elevation) {
	/*
		SUZUKI FUNCTION 
		based on Macedonio et al. (2005)
		this function is coded by K. Mannen (10-JAN-2011)
		z     particle release height
		slice integration step
		SUZUKI_A   suzuki parameter A
		SUZUKI_LAMDA   suzuki parameter lamda
	*/
	double probability = 0;
	static double normalization_constant = 0;
	
	int i = 0;
	double column_h = plume_height - vent_elevation;
	double z_bar_H = 0;
	double zz = z;
	double demon1, demon2, demon3, demon4;
	double ht_step_width;
  
  /* if (!slice) fprintf(stderr, "ENTER plume_pdf0 ....\n"); */
 
	probability = 0.0;
	z_bar_H = (z - vent_elevation)/column_h;
	
	/* calculation of normalization_constant */
	
	if (!normalization_constant) {
		ht_step_width = zz - vent_elevation;
		/*fprintf(log_file,"line 66, z, demon2, demon4, probability\n");*/
		for(i=0; i<S_STEPS; i++){
			z_bar_H = zz / plume_height;
			demon1 = 1 - z_bar_H;
			demon2 = SUZUKI_A * (-demon1);
			demon3 = exp(demon2);
			demon4 = demon1*demon3;
			probability = pow(demon4, SUZUKI_LAMDA);
			/*fprintf(log_file,"line 66 %1.4f, %1.4f, %1.4f, %1.4f\n",zz, demon2, demon4, probability);*/
			/*printf("h = %1.0f probability = %1.4f\n", h, probability);*/
			normalization_constant += probability;
			zz += ht_step_width;
		}
		/*printf("normalization_constant=%1.3e\n", normalization_constant);*/
    }
	
	z_bar_H = (z - vent_elevation)/column_h;
	demon1 = 1 - z_bar_H;
	demon2 = SUZUKI_A * (-demon1);
	demon3 = exp(demon2);
	demon4 = demon1*demon3;
	probability = pow(demon4, SUZUKI_LAMDA);


  /*fprintf(log_file, "line80, %1.0f, %1.4f, %1.4f, %1.4f, %1.4e\n",z, demon1, demon2, demon4, probability);*/
			
	probability = probability/normalization_constant;
	/*fprintf(log_file, "line80, %1.0f, %1.4f, %1.4f, %1.4f, %1.4e\n",z, demon1, demon2, demon4, probability);*/
	/*printf("%1.4e\n", probability);*/
  return probability; 
}

void set_global_values(FILE *log) {

  log_file = log;
#ifdef _PRINT  
  fprintf(log_file, "IN set_global_values ...");
#endif
  /* Set values for global static variables */
  
  AIR_VISCOSITY_x_18 = 18.0 * AIR_VISCOSITY;
  LITHIC_DENSITY_minus_PUMICE_DENSITY = LITHIC_DENSITY - PUMICE_DENSITY;
  PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD = PUMICE_DIAMETER_THRESHOLD - LITHIC_DIAMETER_THRESHOLD;
  ONE_THIRD = 1.0 / 3.0;
  AIR_VISCOSITY_x_225 = AIR_VISCOSITY * 225.0;
  GRAV_SQRD_x_4 = 4.0 * GRAVITY * GRAVITY;
  EDDY_CONST_x_8_div_5 = 8.0 * EDDY_CONST / 5.0;
  T = NULL;
#ifdef _PRINT  
  fprintf(log_file, "OUT");
#endif
}

void set_eruption_values(ERUPTION *erupt, WIND *wind, int *gl_ptr) { /* set_eruption_values */

  /* The following parameters are the properties of a eruption
   * each eruption must have all of these parameters defined:
   *
   * erupt->total_ash_mass is the total amount of ash erupted by
   * the volcano over the course of the entire eruption or calculation period
   * erupt->max_part_size is the maximum particle diameter considered
   * in the calculation. This is input in phi units (so it will likely be
   * a negative number like -5 and appear to be less than min_part_size)
   * erupt->min_part_size is the minimum particle diameter condsidered in the
   * calculation. This input is in phi units.
   *
   * Note: erupt->max/min_part_size are used to set the limits of integration
   * on the calculation. Particles outside this range are not considered at all.
   *
   * erupt->part_mean_size is the mean particle diameter erupted in phi units
   * erupt->part_sigma_size is the standard deviation in particle diameter in phi units
   * erupt-> vent_height is the elevation of the vent m.a.s.l. in meters
   * erupt->max_column_height is the eruption column height m.a.s.l. 
   * (not used) erupt->column_beta is the shape factor governing 
   * the particle size distribution in the eruption column. 
   */

  int i, j;

  double y, prob;
  double x, total_P_col, total_P_part, cum_prob_part, cum_prob_col, total_P;
  double particle_ht, cum_fall_time, wind_x, wind_y, ht_above_vent, temp;
  double col_prob, part_prob;
  double v_phi;

  double fall_v;

  double ht_section_width;
  double part_section_width;
  double ht_step_width;
  double part_step_width;
  double alpha = 0.0;
  double beta = 0.0;
	
  double daemon_sq_sigma, daemon_fine, daemon_coarse;  /* added by Kaz (Dec.02,2010) */

#ifdef _PRINT
  fprintf(log_file, "IN set_eruption_values ... ");
#endif

  PART_STEPS = (erupt->max_phi - erupt->min_phi) * 10;
#ifdef _PRINT
  fprintf(log_file, "PART_STEPS=%d\n", PART_STEPS);
#endif

  threshold = PLUME_RATIO * S_MAX; //PLUME_RATIO * (erupt->max_plume_height - erupt->vent_height);
  //printf("line 738 PR %1.4f Threshold %1.4f Smax %1.4f\n", PLUME_RATIO, threshold, S_MAX);
  BETA_x_SQRT_TWO_PI = erupt->column_beta * sqrt(2.0 * M_PI);
  TWO_BETA_SQRD = 2.0 * erupt->column_beta * erupt->column_beta;
  PDF_GRAINSIZE_DEMON1 = 1.0 / (2.506628 * erupt->sigma_phi);
  TWO_x_PART_SIGMA_SIZE = 2.0 * erupt->sigma_phi * erupt->sigma_phi;
  
	/*printf("PDF_GRAINSIZE_DEMON1 = %1.4f\n", PDF_GRAINSIZE_DEMON1);*/
	/*printf("TWO_x_PART_SIGMA_SIZE = %1.4f\n", TWO_x_PART_SIGMA_SIZE);*/
  
	
  /*define the limits of integration */ 
  ht_section_width = erupt->max_plume_height - erupt->vent_height; 
  part_section_width = erupt->max_phi - erupt->min_phi;
  //ht_step_width = ht_section_width / (double)COL_STEPS; 
  part_step_width = part_section_width / (double)PART_STEPS;

  /* steps for nomalization of probabilities */
  cum_prob_col = 0.0;
  x = erupt->vent_height;
  for (i=0; i < S_STEPS; i++) {
	x = C[i+1].s;
    
	//printf("%f\n", wind[i].wind_height);
	//printf("line 761 x= %f", x);
	prob = (*pdf)(x, (double)i, erupt->max_plume_height, erupt->vent_height); /*modified by Kaz; Jan-2011*/
	cum_prob_col += prob;
    //#ifdef _PRINT
    //fprintf(stderr, " slice_ht=%g, prob=%g, cum_prob=%g\n", x, prob, cum_prob_col); 
    //#endif
  }    
  total_P_col = cum_prob_col;
  //fprintf( stderr, "total_P_col=%g\n ", total_P_col);

  cum_prob_part = 0.0;
  y = (erupt)->min_phi;
  for (i=0; i < PART_STEPS; i++) {
    prob = pdf_grainsize(erupt->mean_phi, y, part_step_width, gl_ptr);
    cum_prob_part += prob;
	
    //fprintf(stderr, " grain_size=%.2f, prob=%g, cum_prob=%g\n", y, prob, cum_prob_part);
    y += part_step_width;
  }
  total_P_part = cum_prob_part;
  //fprintf( stderr, "total_P_part=%g \n", total_P_part);

  /* Normalization constant */
  total_P = (total_P_col * total_P_part);


  /* End of normalization steps */  

  /* Dynamically allocated table for storing integration data.
     Used in the double integration steps below for each point considered.
  */
    if (T == NULL) {
		
		T = (TABLE **)malloc((size_t)PART_STEPS * sizeof(TABLE *));
		
	if (T == NULL) {
        fprintf(log_file, 
        "Cannot malloc memory for Integration Table:[%s]\n", strerror(errno));
        exit(1);
      }
		for (i=0; i<PART_STEPS; i++) {
			T[i] = (TABLE *)malloc((size_t)S_STEPS * sizeof(TABLE));
			
        if (T[i] == NULL) {
          fprintf(log_file, 
          "Cannot malloc memory for Integration Table[%d]:[%s]\n", i, strerror(errno));
          exit(1);
        }
      }
    } else {

		T = (TABLE **)_realloc(T, (size_t)PART_STEPS * sizeof(TABLE *));

      if (T == NULL) {
        fprintf(log_file, 
        "Cannot malloc memory for Integration Table:[%s]\n", strerror(errno));
        exit(1);
      }
		for (i=0; i<PART_STEPS; i++) {
		T[i] = (TABLE *)_realloc(T[i], (size_t)S_STEPS * sizeof(TABLE));

        if (T[i] == NULL) {
          fprintf(log_file, 
          "Cannot malloc memory for Integration Table[%d]:[%s]\n", i, strerror(errno));
          exit(1);
        }
      }
    }
fprintf(log_file,
	      "\nPart_Ht\ts\tAsh_Diam\tPart-Den\tFalltime\tFallV\tDFP\tDCP\tTFalltime\tWsumX\tWsumY\tColumn_x\tColumn_y\treleased_mass\tv_phi\n");
	      
    /* Start with the maximum particle size */
	y = (erupt)->min_phi;
	
	/* set common values to calculate diffusion time; added by Kaz on Dec. 02, 2010 */
	daemon_sq_sigma = pow(SIGMA_PLUME, 2);
	daemon_coarse = daemon_sq_sigma / (4 * DIFFUSION_COEFFICIENT);
	daemon_fine = daemon_sq_sigma * 5 / (8 * EDDY_CONST);
	daemon_fine = pow(daemon_fine, -2.5);
	
    for (i = 0; i < PART_STEPS; i++) { /* PART_STEPS_LOOP */
		/*y += part_step_width;*/
      T[i][0].part_density  =  particle_density(y);    
      T[i][0].ashdiam = phi2m(y);
      part_prob = pdf_grainsize(erupt->mean_phi, y, part_step_width, gl_ptr);
      cum_fall_time = 0.0;
      wind_x = 0.0;
      wind_y = 0.0;
      
      /* Start at the height of the vent */
      particle_ht = erupt->vent_height;
	  
      for (j = 0; j < S_STEPS; j++) { /* COL_STEPS_LOOP */

	      /* define the small slice dz */
	      //particle_ht = ht_step_width;
		  T[i][j].particle_ht = C[j+1].z;
		  particle_ht = T[i][j].particle_ht;
		  
		  if (j > 0){
		  	ht_step_width = T[i][j].particle_ht - T[i][j-1].particle_ht;
		  } else {
		  	ht_step_width = T[i][j].particle_ht - erupt->vent_height;
		  }
	      
	      /* Calculate the time it takes a particle to fall from its release point
	         in the column to the next column release point.
	       */    
	      T[i][j].fall_time = 
	      part_fall_time(particle_ht, ht_step_width, T[i][0].ashdiam, T[i][0].part_density);
		  fall_v = ht_step_width / T[i][j].fall_time;
	      
    	 /* Particle diffusion time (seconds) */
    	 //ht_above_vent = particle_ht - erupt->vent_height; (commented out by Kaz Apr. 28, 2017)
	     
			//temp = 0.2 * ht_above_vent * ht_above_vent;
			//temp = 5 / (8 * EDDY_CONST) * pow(0.34/3, 2) * ht_above_vent * ht_above_vent; // fixed above bug by Kaz. (Feb 10, 2017)
		 	temp = 5 / (8 * EDDY_CONST) * C[j].radius * C[j].radius; // edited by Kaz (Feb 10, 2017)
			
			if(SIGMA_PLUME==0){								/* this if statement was added by Kaz Dec-02-'10*/
				/* if particle release occurs in a manner of Suzuki Function, plume_diffusion time would change
				depending on the height (uprising column or umbrella cloud)  */
				
			T[i][j].plume_diffusion_fine_particle = 		/* eq.9 of Bonadonna et al. (2005)*/
			pow(temp, 0.4); /* 0.4 = 2.0/5.0 */
       
			//T[i][j].plume_diffusion_coarse_particle = 		/* eq.7 of Bonadonna et al. (2005)*/
			//0.0032 * (ht_above_vent *  ht_above_vent) / DIFFUSION_COEFFICIENT;
			
			T[i][j].plume_diffusion_coarse_particle = 		/* eq.7 of Bonadonna et al. (2005); edited by Kaz (Feb 10, 2017)*/
			(C[j+1].radius * C[j+1].radius) / (DIFFUSION_COEFFICIENT * 4);
			
			}
			else{
				if(C[j+1].s <= threshold){
						T[i][j].plume_diffusion_fine_particle = /* eq.9 of Bonadonna et al. (2005)*/
						pow(temp, 0.4); /* 0.4 = 2.0/5.0 */
       
						//T[i][j].plume_diffusion_coarse_particle = 	/* eq.7 of Bonadonna et al. (2005)*/
						//	0.0032 * (ht_above_vent *  ht_above_vent) / DIFFUSION_COEFFICIENT;
						
						T[i][j].plume_diffusion_coarse_particle = 		/* eq.7 of Bonadonna et al. (2005)*/
						(C[j+1].radius * C[j+1].radius) / (DIFFUSION_COEFFICIENT * 4); //edited by Kaz (Feb 10, 2017)
				}
				else{
				/* if characteristic radius of column is given, plume diffusion time is set using C or K. */
				daemon_fine = (daemon_sq_sigma * 5  / (8 * EDDY_CONST));
				daemon_fine = pow(daemon_fine, 0.4);
				daemon_coarse = daemon_sq_sigma / (4 * DIFFUSION_COEFFICIENT);
				T[i][j].plume_diffusion_fine_particle = daemon_fine;
				T[i][j].plume_diffusion_coarse_particle = daemon_coarse;
				}
			}
			
#ifdef _PRINT
    //fprintf(log_file, " i=%g, j=%g, cum_prob=%g\n", y, prob, cum_prob_part);
#endif
	     
	      /* Sum the windspeed and wind_direction for each particle size 
	       * falling from each level. In the wind array, the first wind level
	       * gives wind speed and direction at the vent height. 
	       * Start with the next wind level, 
	       * so that we are using the wind speed and direction 
	       * starting from one step above the vent. 
	       */
     
	     wind_x += 
	     T[i][j].fall_time * wind[j+1].wind_speed * cos(wind[j+1].wind_dir);
	     
	     wind_y += 
	     T[i][j].fall_time * wind[j+1].wind_speed * sin(wind[j+1].wind_dir);
	      
	     T[i][j].wind_sum_x = wind_x;
         T[i][j].wind_sum_y = wind_y;
	     
		 // Meaning of i or j = 0 is following
		 // s = 0   at C[0], T[–][j], and wind[1]
		 // s = 100 at C[1], T[0][j], and wind[2]
		 // Column data inlude that of at the source vent, while TABLE (T) does not and starts from the next step of the source
		 
		 
	      /* Accumulate the time it takes each particle size to descend
	         from its release point down
	         to its final resting place.This part of the code just 
	         calculates the fall_time from the release point to the 
	         height of the vent.
	         The time it takes a particle to fall from the vent height 
	         to a grid cell will be calculated later. 
	       */
	    cum_fall_time += T[i][j].fall_time;
	    T[i][j].total_fall_time = cum_fall_time;
        col_prob = (*pdf)(C[j+1].s, j, erupt->max_plume_height, erupt->vent_height);
		/* Select Particle Segregation Pattern (or Source Magnitude Distribution; SMD)*/
		if(PLUME_THICKNESS == -9999 && SOURCE_DECAY_RATE == -9999){
			T[i][j].demon1 = (erupt->total_ash_mass * col_prob  * part_prob) / total_P;
		}else if(PLUME_THICKNESS > 0 && SOURCE_DECAY_RATE > 0){
	        fprintf(stderr, 
	    	      "ERROR\nYou cannot assign both of PLUME_THICKNESS and SOURCE_DECAY_RATE at once. Assign either or none of both.\nPROGRAM HAS BEEN HALTED tephra_calc.c L959\n\n");
	        exit(1);
		}else if(PLUME_THICKNESS > 0){
			v_phi = 1 / part_fall_time(PLUME_HEIGHT, 1, T[i][0].ashdiam, T[i][0].part_density); // terminal velocity of particle
			beta = v_phi / (PLUME_THICKNESS * WIND_HT);
			alpha = exp(-1 * beta * C[j].s) - exp(-1 * beta * C[j+1].s);
			T[i][j].demon1 = (erupt->total_ash_mass * part_prob) / total_P * alpha;
		}else if(SOURCE_DECAY_RATE > 0){
			beta = SOURCE_DECAY_RATE;
			alpha = exp(-1 * beta * C[j].s) - exp(-1 * beta * C[j+1].s);
			T[i][j].demon1 = (erupt->total_ash_mass * part_prob) / total_P * alpha;
		}else{
	        fprintf(stderr, 
	    	      "ERROR\nAssign appropriate value for PLUME_THICKNESS or SOURCE_DECAY_RATE or remain both -9999 for uniform particle release.\nPROGRAM HAS BEEN HALTED tephra_calc.c L959\n\n");
	        exit(1);
		}
        /* Normalization is now done here */
        
        
	      //T[i][j].particle_ht = C[j].z;
	      	
	      fprintf(log_file,
	      "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      T[i][j].particle_ht,
		  C[j+1].s,
	      T[i][j].ashdiam, 
	      T[i][j].part_density, 
	      T[i][j].fall_time,
		  fall_v,
	      T[i][j].plume_diffusion_fine_particle,
	      T[i][j].plume_diffusion_coarse_particle,
	      T[i][j].total_fall_time,
	      T[i][j].wind_sum_x,
	      T[i][j].wind_sum_y,
		  C[j+1].centre_x,
		  C[j+1].centre_y,
		  T[i][j].demon1,
		  v_phi); //T[i][j].demon1);
      } /* END COL_STEPS_LOOP */ 
      
      fprintf(log_file, "\n");
      y += part_step_width;     /*20101203*/
    } /* END PART_STEPS_LOOP */
  	fprintf(log_file, "OUT\n");
}

void read_column_file(int step_column){ // moved from windy.c on Feb. 10, 2017 by Kaz; called by windy.c
	int i=0;
	int base=1;
	int ret, ttlline;
	double	z, s, x, 
			dir, northing, easting, 
			ta, p, dens_atm, 
			dens_column, gas_cont, Q, 
			Cp, Rg, V, 
			Ue, M, theta, 
			U, R, T;
	
	C = (COLUMN *)malloc((step_column + 1) * sizeof(COLUMN)); // added 1; sea note of 2018.3.21
	
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
						
						if(i>=base){
							C[i-base].centre_x = northing;
							C[i-base].centre_y = easting;
							C[i-base].radius = R;
						
							C[i-base].s = s;
							C[i-base].z = z;
							C[i-base].wind_dir = dir;
							C[i-base].wind_v = V;
							//printf("L1055 i = %d\t%g\n", i, s);
						}
		}
		
		i++;
	}
	fclose(op);
	
}

void free_memory_tephra_calc(void){
	free(C);
	free(T);
}