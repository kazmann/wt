typedef struct {
  /* geographic location of a point where ash mass in estimated */
  double easting; /*utm coordinate in meters */
  double northing; /*utm coordinate in meters */
  double elevation; /*elevation at the point in masl */
  double mass;  /*mass of material accumulated at the point (gm) */
  //double md_phi; /*median grainsize at the point (phi) */
  //double sigma_phi; /*std deviation of grainsize (sorting) at the point (phi) */
  double phi[20]; /* the grainsize distribution at this location */
  double cum_mass; /* used to accumulate mass from multiple eruptions */
  double *mass_pt; /* pointer to an array of mass accumulations, one for each eruption */
  char pnt_name[8]; /* point name (added on Mar 20 2014 by Kaz) */
} POINT;

/* These values are calculated ahead of time and used in the double integration [PART_STEPS][COL_STEPS]*/
typedef struct {
  double particle_ht;
  double ashdiam;
  double part_density;
  double fall_time;							// total fall time from particle_ht to 0
  double plume_diffusion_fine_particle;
  double plume_diffusion_coarse_particle;
  double total_fall_time;
  double wind_sum_x;						// drift distance by wind
  double wind_sum_y;						// drift distance by wind
  double demon1;
} TABLE;

typedef struct {
  /* eruption parameters */

  double volcano_northing; /*volcano location in UTM north (meters) */
  double volcano_easting; /*volcano location in UTM east (meters) */
  double total_ash_mass; /* is the total amount of ash erupted (kg) */
  double min_phi;  /*the maximum particle diameter considered (phi)*/
  double max_phi; /* is the minimum particle diameter considered (phi) */
   
  /*Note: erupt->max/min_part_size are used to set the limits of integration
   * on the calculation. Particles outside this range are not considered at all */

  /*Note: phi uits are such that max will appear to be less than min, this is
    accounted for in the conversion to cm, which is internal */

  double mean_phi;  /*the mean particle diameter erupted in phi units */
  double sigma_phi; /*standard deviation in particle diameter in phi units */
  double vent_height;  /* elevation of the vent amsl in meters */
  double max_plume_height;  /*eruption column height amsl in meters */
 
  double column_beta;  /* parameter governing the particle size distribution in column */

  /* Note: A large value  of beta (1) places most of the particles
   * high in the eruption column, a low value of beta (0.01) spreads the particle density
   * lower in the column. Particle release models based on "corner" models etc strongly
   * suggest that a larger value for beta should be used. */

} ERUPTION;

typedef struct {
  int day;
  int hour;
  double wind_height; /* height a.s.l. in km */
  double wind_speed; 	/* the average windspeed in m/s */
  double wind_dir;  	/* average wind direction in +/- degrees from north */
  double t_atm;
  double p_atm;
} WIND;

typedef double (*PFI)(double, double, double);

typedef struct {
  double min_falltime;
  double max_falltime;
} STATS;

/* added by K. Mannen (19-Jan-2011) */
typedef struct {
  double sizeclass;
  double content;
} GRAIN;
/*     +++++     +++++     +++++    */

/* added by K. Mannen (26-Jan-2017) */
typedef struct {
  double s;
  double z;
  double wind_dir;
  double wind_v;
  double centre_x;	/* x-coordinate of centre position */
  double centre_y;	/* y-coordinate of centre position */
  double radius;
} COLUMN;
/*     +++++     +++++     +++++    */