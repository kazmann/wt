/*
	windy.c

	A code to simulate DRY volcanic eruption column in windy condition
	Based on Woodhouse et al. [2013] Interaction between volcanic plumes and wind during the 2010
									Eyjafjallajokull eruption Iceland, J. Geophys. Res. V118 92-109

	coded by Kaz Mannen (2017 Jan)

	func## in this code is correlated to eq## in Woodhouse et al. [2013]
	modification on Jan 19, 2022
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "prototypes_strat.h"


//#define _malloc malloc

// Global

double z, T;

// atmosphere structure
// int MAX_LINE=100;
int STEP=10;

//double SLICE_HEIGHT=100;

double lps_strt = 2/1000;
double lps_trop = 6.5/1000;
double H1 = 11000;
double H2 = 20000;

double V1 = 40;	//discarded
double g_dir;

double Ra =  285;
double Rg0 =  462;
//double g =  9.81;

double z0 =  0;
double t0 =  293;	//293;
double x, north, east;	// horizontal position; x means max length

double p = 100000;	//101325.0;

double ds = 100;
double dz;
double U = -9999;
double R = -9999;
double n0;
double T;
double WIND_HT = -9999;

double rho_s = 1200;
double rho_w = 1000;

double Ca = 998; //713;
double Cs = 1617; //1100;
double Cv = 1850; //1850

//double ks = 0.09;  // another k should be introduced for gas thurst region but uniform value in this code
//double kw = 0.9;

double theta = M_PI / 2; // PI / 2

// param in func11
double E, Cp, Cp0;

// param in func12 - 15
double M, Q, rho_c, rho_a, Ue, V;

// param in func16
int flag=0; //	0, gas-thrust;   1, buoyant;     2, umbrella
			//  rho_a > rho_c    rho_a < rho_c   rho_a > rho_c

// param in func17-19
double n, Q0, Rg;


double z, s;
double ta, dp_over_dz;
double max;
static WIND *W1;


//////////////////

void windy(FILE *in_wind){	// The main routine in this file
	int i = 0;
	int total = 0; // total line number of wind file
	
	double Hg, Hb, Ht;
	double U0, R0;
	double U_Ht, T_Ht, R_Ht, x_Ht;
	double Q_hb, rho_c_hb, t_hb;
	
	double g_dir_previous, north_previous, east_previous, ta_previous, p_previous, rho_a_previous, rho_c_previous, n_previous;
	double Q_previous, Cp_previous, Rg_previous, M_previous;
	
	int step_column;
	
	n0 = INITIAL_WATER_CONTENT;
	z0 = VENT_ELEVATION;
	
	Q_hb = 0, rho_c_hb = 0;
	
	total=read_wind(in_wind);	// read wind file and store in W1
	
	z = z0;
	s = 0;
	FILE *f, *f2;
	
	T = MAGMA_TEMPERATURE;
	
	f = fopen("column.txt", "w");
	
	// initialize
	
	//ta = func21(z);
	ta = calc_Tatm(z, total);
	p = calc_Patm(z, total);
	Cp0 = calcCp0();	// func 20.5
	
	rho_a = func23(p, ta);
	n=n0;
	
	Rg=Rg0;
	rho_c=func17(n0, p, Rg, T); // get rho_c
	rho_c_hb = 0;
	
	if(n0 < 0 || n0 > 1){
      fprintf(stderr, 
  	      "ERROR\nYou need proper INITIAL_WATER_CONTENT in config file\nPROGRAM HAS BEEN HALTED\n\n");
      exit(1);
	}
	
	if(MAGMA_DISCHARGE_RATE < 0 || INITIAL_PLUME_VELOCITY < 0 || VENT_RADIUS < 0){
		if(MAGMA_DISCHARGE_RATE < 0 && INITIAL_PLUME_VELOCITY > 0 && VENT_RADIUS > 0){
			Q = rho_c * INITIAL_PLUME_VELOCITY * VENT_RADIUS * VENT_RADIUS;
			U = INITIAL_PLUME_VELOCITY;
			R = VENT_RADIUS;
		}else if(MAGMA_DISCHARGE_RATE > 0 && INITIAL_PLUME_VELOCITY < 0 && VENT_RADIUS > 0){
			Q = MAGMA_DISCHARGE_RATE / M_PI;	// mass flux is defined as pi * Q in Woodhouse et al. (2013)
			U = Q / (rho_c * VENT_RADIUS * VENT_RADIUS);
			R = VENT_RADIUS;
		}else if(MAGMA_DISCHARGE_RATE > 0 && INITIAL_PLUME_VELOCITY > 0 && VENT_RADIUS < 0){
			Q = MAGMA_DISCHARGE_RATE / M_PI;	// mass flux is defined as pi * Q in Woodhouse et al. (2013)
			U = INITIAL_PLUME_VELOCITY;
			R = sqrt(Q / (rho_c * INITIAL_PLUME_VELOCITY));
		}else{
	        fprintf(stderr, 
	    	      "ERROR\nYou need to assign at least two parameters properly from U, R and Q in the config file\nPROGRAM HAS BEEN HALTED 179\n\n");
	        exit(1);
		}
	}else{
        fprintf(stderr, 
    	      "ERROR\nYou need to assign at least two parameters properly from U, R and Q in the config file\nPROGRAM HAS BEEN HALTED 184\n\n");
        exit(1);
	}

	
	// initialize (func 11)
	// Q = rho_c * U * R * R;
	M = rho_c * U * U * R * R;
	E = Q * Cp0 * T;
	
	Q0 = Q;
	U0 = U;
	R0 = R;
	
	Cp = Cp0;

	//V=func24(z);
	V = get_V(z, total);		// wind velocity
	g_dir = get_dir(z, total);	// wind direction
	x = 0;
	north = 0;
	east = 0;
	
	fprintf(f, "#z\ts\tx\tdir\tnorthing\teasting\tTa\tP\tatm_dens\tcol_dens\tn\t");
	fprintf(f, "Q\tCp\tRg\tV\tUe\tM\ttheta\t");
	fprintf(f, "U\tR\tT\n");
	
	if(S_MAX < 0){max = 99999;}else{max = S_MAX;}
	
	// Step until plume reaches to the Ht
	while(i < 100000 && M > 0.0 && theta > 0.0 && theta < M_PI && s <= max){
		Ht = z; // when M < 0 (static) or theta < 0 (windy), z just before the height is considered as Ht
		WIND_HT = V;
		U_Ht = U, T_Ht = T, R_Ht = R, x_Ht = x;
		// top of the gas thrust region
		if(flag==0 && rho_a - rho_c > 0){flag=1; Hg = z;}	// top of the gas-thrust region
		// top of the convective region
		if(flag==1 && rho_a - rho_c < 0){flag=2; Hb = z; Q_hb = Q; rho_c_hb = rho_c; t_hb = T;}	// top of the convective region
		fprintf(f, "%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", z, s, x, g_dir, north, east, ta, p, rho_a, rho_c, n);
		fprintf(f, "%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", Q, Cp, Rg, V, Ue, M, theta);
		fprintf(f, "%1.4f\t%1.4f\t%1.4f\n", U, R, T);
		
		g_dir_previous = g_dir; north_previous = north; east_previous = east; 
		ta_previous = ta; p_previous = p; rho_a_previous = rho_a; rho_c_previous = rho_c; n_previous = n;
		north_previous = north; east_previous = east;
		Q_previous = Q; Cp_previous = Cp; Rg_previous = Rg; M_previous = M;
		
		rk(total);
		i++;
	}
	U = WIND_HT, R = R_Ht, T = T_Ht, x = x_Ht;
	//fclose(f);
	
	f2 = fopen("column_parameters.txt", "w");
	
	//step_column = (int)((Ht - z0) / SLICE_HEIGHT);	// column step that is a parameter in tephra 2 is set to give each slice is about 100m
	fprintf(f2, "#Q0\tU0\tR0\tHg\tHb\tHt\tColumnT\tAtmT\tCOL_STEP\tn_final\tQ_final\tQ_hb\trho_c_hb\tt_hb\n");
	//fprintf(f2, "%1.4e\t%1.4e\t%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n", Q0 * M_PI, U0, R0, Hg, Hb, Ht, T_final, Ta_final, step_column, n_final, Q_final, Q_hb, rho_c_hb, g_dir_final);
	fprintf(f2, "%1.4e\t%1.4e\t%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n", Q0 * M_PI, U0, R0, Hg, Hb, Ht, T, ta, step_column, n, Q, Q_hb, rho_c_hb, t_hb);
	
	fclose(f2);

	COL_STEPS = step_column;
	PLUME_HEIGHT = Ht;
	
	s = s - ds;
	S_HT_MAX = s;
	S_HT_STEPS = i - 1;
	
	if(S_MAX < 0){
		S_MAX = S_HT_MAX;
		S_STEPS = S_HT_STEPS;
	}else{
		S_STEPS = S_MAX / ds;
	}
	
	g_dir = g_dir_previous; north = north_previous; east = east_previous; 
	ta = ta_previous; p = p_previous; rho_a = rho_a_previous; rho_c = rho_c_previous; n = n_previous;
	north = north_previous; east = east_previous;
	Q = Q_previous; Cp = Cp_previous; Rg = Rg_previous; M = M_previous;
	
	// After reaching Ht
	theta = 0; Ue = 0;
	for (i = i;  i < S_STEPS + 1; i++){
		s += ds;
		x += ds;
		north = north + ds * cos(g_dir / 360 * 2 * M_PI);
		east = east + ds * sin(g_dir / 360 * 2 * M_PI);
		fprintf(f, "%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", Ht, s, x, g_dir, north, east, ta, p, rho_a, rho_c, n);
		fprintf(f, "%1.4e\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t", Q, Cp, Rg, WIND_HT, Ue, M, theta); // V becomes WIND_HT
		fprintf(f, "%1.4f\t%1.4f\t%1.4f\n", U, R, T); // U becomes WIND_HT
	}
	
	fclose(f);
	//printf("kokodayo %1.4f\t%d\n", S_MAX, S_STEPS);
	read_column_file(S_STEPS);
	free(W1); //20180218
}

double thetacheck(double tmp){
	if(tmp < 0){tmp = 0;}
	
	return tmp;
}

void rk(int total){
	double dp_over_ds, dQ_over_ds, dM_over_ds, dtheta_over_ds, dE_over_ds;
	double dp_over_ds1, dQ_over_ds1, dM_over_ds1, dtheta_over_ds1, dE_over_ds1;
	double dp_over_ds2, dQ_over_ds2, dM_over_ds2, dtheta_over_ds2, dE_over_ds2;
	double dp_over_ds3, dQ_over_ds3, dM_over_ds3, dtheta_over_ds3, dE_over_ds3;
	double dp_over_ds4, dQ_over_ds4, dM_over_ds4, dtheta_over_ds4, dE_over_ds4;
	double dx;
	
	double E_tmp, M_tmp, n_tmp, rho_a_tmp, rho_c_tmp, p_tmp, theta_tmp, Q_tmp;
	
	E_tmp = E; M_tmp = M; n_tmp = n; rho_a_tmp=rho_a; rho_c_tmp=rho_c; theta_tmp = theta; Q_tmp = Q;


	/////////////////////////////////
	// k1     ///////////////////////
	dp_over_ds1 = func22(p, ta);
	dQ_over_ds1 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds1 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds1 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds1 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds1);

	/////////////////////////////////////////////////
	// generate next step parameters-----------------
	Q_tmp = Q + dQ_over_ds1 * ds * 0.5;	
	M_tmp = M + dM_over_ds1 * ds * 0.5;
	theta_tmp = theta + dtheta_over_ds1 * ds * 0.5;
	theta_tmp = thetacheck(theta_tmp);
	E_tmp = E + dE_over_ds1 * ds * 0.5;
	
	// centre position of the next step
	dz = ds * sin(theta_tmp) * 0.5;
	
	
	// atmosphreic content of the next step
	n = func18(Q_tmp);				// calc n
	p_tmp = calc_Patm(z + dz, total);
	ta = calc_Tatm(z + dz, total);  //func21(z + dz);
	rho_a_tmp = func23(p_tmp, ta);
	
	Cp = func20(n);					// calc Cp	
	//printf("k1\n"); 
	V  = get_V(z + dz, total);
	
	T = E_tmp / Q_tmp / Cp;
	
	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc column density (rho_c)
	
	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));
	
	/////////////////////////////////
	// k2     ///////////////////////
	dp_over_ds2 = func22(p_tmp, ta);
	dQ_over_ds2 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds2 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds2 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds2 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds1);

	/////////////////////////////////////////////////
	// generate next step parameters-----------------
	Q_tmp = Q + dQ_over_ds2 * ds * 0.5;	
	M_tmp = M + dM_over_ds2 * ds * 0.5;
	theta_tmp = theta + dtheta_over_ds2 * ds * 0.5;
	theta_tmp = thetacheck(theta_tmp);
	E_tmp = E + dE_over_ds2 * ds * 0.5;
	
	// centre position of the next step
	dz = ds * 0.5 * sin(theta_tmp);
	
	
	// atmosphreic content of the next step
	n = func18(Q_tmp);				// calc n
	p_tmp = calc_Patm(z + dz, total);
	ta = calc_Tatm(z + dz, total);  //func21(z + dz);
	rho_a_tmp = func23(p_tmp, ta);
	
	Cp = func20(n);					// calc Cp	
	//V=func24(z);
	//printf("k2\n"); 
	V = get_V(z + dz, total);
	
	T = E_tmp / Q_tmp / Cp;
	
	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc column density (rho_c)
	
	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));
	
	/////////////////////////////////
	// k3     ///////////////////////
	dp_over_ds3 = func22(p_tmp, ta);
	dQ_over_ds3 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds3 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds3 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds3 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds2);

	/////////////////////////////////////////////////
	// generate next step parameters-----------------
	Q_tmp = Q + dQ_over_ds3 * ds;	
	M_tmp = M + dM_over_ds3 * ds;
	theta_tmp = theta + dtheta_over_ds3 * ds;
	theta_tmp = thetacheck(theta_tmp);
	E_tmp = E + dE_over_ds3 * ds;
	
	// centre position of the next step
	dz = ds * sin(theta_tmp);
	
	
	// atmosphreic content of the next step
	n = func18(Q_tmp);				// calc n
	p_tmp = calc_Patm(z + dz, total);
	ta = calc_Tatm(z + dz, total);  //func21(z + dz);
	rho_a_tmp = func23(p_tmp, ta);
	
	Cp = func20(n);					// calc Cp	
	//V=func24(z);
	//printf("k3\n"); 
	V = get_V(z + dz, total);
	
	T = E_tmp / Q_tmp / Cp;
	
	Rg = func19(n);					// calc Rg
	rho_c_tmp = func17(n, p_tmp, Rg, T);	// calc column density (rho_c)
	
	Ue = func16(M_tmp, Q_tmp, theta_tmp, V);
	U = M_tmp / Q_tmp;
	R = sqrt(Q_tmp / (U * rho_c));
	
	/////////////////////////////////
	// k4     ///////////////////////
	dp_over_ds4 = func22(p_tmp, ta);
	dQ_over_ds4 = func12(M_tmp, rho_a_tmp, rho_c_tmp, Q_tmp);
	dM_over_ds4 = func13(rho_a_tmp, rho_c_tmp, M_tmp, theta_tmp, Q_tmp);
	dtheta_over_ds4 = func14(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, theta_tmp);
	dE_over_ds4 =     func15(M_tmp, Q_tmp, rho_a_tmp, rho_c_tmp, ta, theta, dQ_over_ds3);
	
	////////////////////////////////////
	////////////////////////////////////
	//// set new step value   //////////
	
	dp_over_ds = 	 (dp_over_ds1 + 2 * dp_over_ds2 + 2 * dp_over_ds3 + dp_over_ds4)/6;
	dQ_over_ds = 	 (dQ_over_ds1 + 2 * dQ_over_ds2 + 2 * dQ_over_ds3 + dQ_over_ds4)/6;
	dM_over_ds = 	 (dM_over_ds1 + 2 * dM_over_ds2 + 2 * dM_over_ds3 + dM_over_ds4)/6;
	dtheta_over_ds = (dtheta_over_ds1 + 2 * dtheta_over_ds2 + 2 * dtheta_over_ds3 + dtheta_over_ds4)/6;
	dE_over_ds =     (dE_over_ds1 + 2 * dE_over_ds2 + 2 * dE_over_ds3 + dE_over_ds4)/6;
	

	//printf("U = %1.1f\tdQ_over_ds=%1.4f\tdM_over_ds=%1.4f\n", U, dQ_over_ds, dM_over_ds);
	
	// generate next step parameters
	Q = Q + dQ_over_ds * ds;	
	M = M + dM_over_ds * ds;
	if(M<0){M=0;}
	theta = theta + dtheta_over_ds * ds;
	E = E + dE_over_ds * ds;
	
	// centre position of the next step
	dx = ds * cos(theta);
	dz = ds * sin(theta);
	north = north + dx * cos(g_dir / 360 * 2 * M_PI);
	east = east + dx * sin(g_dir / 360 * 2 * M_PI);
	x = x + dx;
	s = s + ds;
	z = z + dz;
	
	// atmosphreic content of the next step
	n = func18(Q);				// calc n
	p = calc_Patm(z, total);
	ta = calc_Tatm(z, total);  //func21(z + dz);
	rho_a = func23(p, ta);
	
	Cp = func20(n);					// calc Cp	
	//V=func24(z);
	//printf("k4\n"); 
	V = get_V(z, total);
	g_dir = get_dir(z, total);


	
	T = E / Q / Cp;
	
	Rg = func19(n);					// calc Rg
	rho_c = func17(n, p, Rg, T);	// calc column density (rho_c)
	
	Ue = func16(M, Q, theta_tmp, V);
	U = M / Q;
	R = sqrt(Q / (U * rho_c));
}

double func12(double M_tmp, double rho_a_tmp, double rho_c_tmp, double Q_tmp){			// column mass flux
	double dQ_over_ds;
	
	dQ_over_ds = 2 * rho_a_tmp * Ue * Q_tmp / sqrt(rho_c_tmp * M_tmp);
	
	//printf("%1.4f\n", dQ_over_ds);
	return dQ_over_ds;
}

double func13(double rho_a_tmp, double rho_c_tmp, double M_tmp, double theta_tmp, double Q_tmp){
	double dM_over_ds;
	
	dM_over_ds = GRAVITY * (rho_a_tmp - rho_c) * Q_tmp * Q_tmp / (rho_c_tmp * M_tmp) * sin(theta_tmp);
	dM_over_ds = dM_over_ds + 2 * rho_a_tmp * Q_tmp / sqrt(rho_c_tmp * M_tmp) * Ue * V * cos(theta_tmp);
	
	return dM_over_ds;
}

double func14(double M_tmp, double Q_tmp, double rho_a_tmp, double rho_c_tmp, double theta_tmp){
	double dtheta_over_ds;
	
	dtheta_over_ds = GRAVITY * (rho_a_tmp - rho_c_tmp) * Q_tmp * Q_tmp * cos(theta_tmp) / (rho_c_tmp * M_tmp * M_tmp);
	dtheta_over_ds = dtheta_over_ds - 2 * rho_a_tmp * Q_tmp * Ue * V * sin(theta_tmp) / (M * sqrt(rho_c_tmp * M_tmp));
	
	//printf("dtheta = %1.4f\n", dtheta_over_ds);
	return dtheta_over_ds;
}

double func15(double M_tmp, double Q_tmp, double rho_a_tmp, double rho_c_tmp, double Ta, double theta_tmp, double dQ_over_ds){
	double dE_over_ds;
	double term1, term2, term3, term4;
	
	term1 = (Ca * Ta + Ue * Ue / 2) * dQ_over_ds;
	term2 = M_tmp * M_tmp / (2 * Q_tmp * Q_tmp) * dQ_over_ds;
	term3 = rho_a_tmp / rho_c_tmp * Q_tmp * GRAVITY * sin(theta_tmp);
	term4 = 2 * rho_a_tmp * Ue * V * cos(theta_tmp) * sqrt(M_tmp / rho_c_tmp);
	dE_over_ds = term1 + term2 - term3 - term4;
	
	return dE_over_ds;
}

double func16(double M_tmp, double Q_tmp, double theta_tmp, double V_tmp){
	double ue_tmp;
	double ks_tmp;
		
	//if(flag==0){ks_tmp=sqrt(rho_a/rho_c)/16;}
	//else{ks_tmp=ks;}	// use these lines when you use ks for gas thrust region; include rho_a and rho_c as local
	
	ks_tmp=KS_ENTRAIN_U;	// gas thrust region also assumes 0.09
				// remove this line when you take
				// ks = f(rho_a. rho_c)
		
	//printf("flag=%d\tks=%1.4f\n", flag, ks_tmp);
	
	ue_tmp = ks_tmp * fabs(M_tmp/Q_tmp - V_tmp * cos(theta_tmp)) + KW_ENTRAIN_V * fabs(V_tmp * sin(theta_tmp));
	
	return ue_tmp;
}

double func17(double n_tmp, double p_tmp, double Rg_tmp, double T_tmp){			// column density
	double rho_tmp;
	
	rho_tmp = (1 - n_tmp) / rho_s + n_tmp * Rg_tmp * T_tmp / p_tmp;
	rho_c = 1 / rho_tmp;
	//printf("rho     = %1.4f\n", rho);
	return rho_c;
}

double func18(double Q_tmp){			// solid content in the column
	double n_tmp;
	n_tmp = 1 - (1 - n0) * Q0 / Q_tmp;
	
	return n_tmp;
}

double func19(double n_tmp){
	double Rg_tmp;
	Rg_tmp = Ra + (Rg0 - Ra) * n0 * (1 - n_tmp) / (n_tmp * (1 - n0));
	
	return Rg_tmp;
}

double func20(double n_tmp){
	double Cp_tmp;
	Cp_tmp = Ca + (Cp0 - Ca) * (1 - n_tmp) / (1 - n0);
	
	return Cp_tmp;
}

double calcCp0(){
	double Cp0_tmp;
	Cp0_tmp = n0 * Cv + (1 - n0) * Cs;
	return Cp0_tmp;
}

double func21(double z){	// atmospheric temperature
	double t;
	
	if(z < H1){
		t = t0 - lps_trop * z;
	} else if (z < H2){
		t = t0 - lps_trop * H1;
	} else {
		t = t0 - lps_trop * H1 + lps_trop * (z - H2);
	}
	return t;
}

double calc_Tatm(double z, int total){	//return wind velocity based on discrete wind data
	int i=1;
	double t_atm=9999.9;
	
	t_atm=W1[total-1].t_atm;
	
	while(i<total){
		if(z < W1[i].wind_height){
			t_atm = W1[i-1].t_atm + (W1[i].t_atm - W1[i-1].t_atm) * (z - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			break;
		}
		i++;
	}
	//printf("%d\t%1.4f\t%1.4f\n", i, z, v);
	return t_atm;
}

double calc_Patm(double z, int total){	//return wind velocity based on discrete wind data
	int i=1;
	double t_atm0=9999.9, t_atm=9999.9;
	double p_atm0=9999.9, p_atm=9999.9;
	double a;
	
	t_atm=W1[total-1].t_atm;
	t_atm0=W1[total-1].t_atm;
	p_atm0=W1[total-1].p_atm;
	
	while(i<total){
		if(z < W1[i].wind_height){
			t_atm = W1[i-1].t_atm + (W1[i].t_atm - W1[i-1].t_atm) * (z - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			t_atm0 = W1[i-1].t_atm;
			p_atm0 = W1[i-1].p_atm;
			break;
		}
		i++;
	}
	
	a = (z - W1[i-1].wind_height) * GRAVITY * 2 / (Ra * (t_atm + t_atm0));
	p_atm = p_atm0 / exp(a) * 100; // hPa -> Pa
	//printf("z = %1.4f\ta_0 = %1.4f\tta_0 = %1.4f\tta_1 = %1.4f\te = %1.4f\tp = %1.4f\tp0 = %1.4f\n", z, a, t_atm0, t_atm, exp(a), p_atm, p_atm0);
	return p_atm;
}

double func22(double p_tmp, double t_tmp){	// atmospheric pressure
	double dp_over_ds;
	
	dp_over_ds = -1 * (GRAVITY * p_tmp) / (Ra * t_tmp);
	return dp_over_ds;
}

double func23(double p_tmp, double t_tmp){	// atmospheric density
	double rho_tmp;
	
	rho_tmp = p_tmp / (Ra * t_tmp);
	return rho_tmp;
}

double func24(double zz){	// DISCARDED; wind velocity change
	double v_tmp;
	
	if (zz < H1) {
		v_tmp = V1 * zz / H1;
	} else {
		v_tmp = V1;
	}
	return v_tmp;
}

double get_V(double z, int total){	//return wind velocity based on discrete wind data
	int i=1;
	double v=9999.9;
	
	v=W1[total-1].wind_speed;
	
	while(i<total){
		if(z < W1[i].wind_height){
			v = W1[i-1].wind_speed + (W1[i].wind_speed - W1[i-1].wind_speed)*(z - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			break;
		}
		i++;
	}
	//printf("%d\t%1.4f\t%1.4f\n", i, z, v);
	return v;
}

double get_dir(double z_tmp, int total){	//return wind direction based on discrete wind data
	int i=1;
	double dir=9999.9;
	double ratio, wind1, wind2;
	
	dir=W1[total-1].wind_dir;
	
	while(i<total){
		wind1 = W1[i-1].wind_dir;
		wind2 = W1[i].wind_dir;
		
		if(z < W1[i].wind_height){
			ratio = (z_tmp - W1[i-1].wind_height) / (W1[i].wind_height - W1[i-1].wind_height);
			if(wind2 - wind1 > 180){
				wind1 = wind1 + 360; 
			}else if(wind1 - wind2 > 180){
				wind2 = wind2 + 360; 
			}
			
			dir = ratio * (wind2 - wind1) + wind1;
			
			if(dir>360){dir=dir-360;}
			if(dir<0){dir=dir+360;}
			
			break;
		}
		i++;
	}
	//printf("dir\t%d\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n", i, z, ratio, dir, wind1, wind2);
	return dir;
}

int read_wind(FILE *in_wind){
	int i=0, z, ret, i2=0, total=0, count=0;
	double wind_height, wind_speed, wind_dir;
	double wind_temp, wind_pres;
	double ratio;
	double delta_wind, wind_tmp, wind1, wind2;
	char line[100];
	
	// Modified following lines till rewind(in_wind) on Jan 17, 2022
	while(NULL != fgets(line, MAX_LINE, in_wind)){
		if(line[0] == '#')continue;
		i++;
	}
	
	double height[i+1], speed[i+1], dir[i+1];
	double atm_temp[i+1], atm_pres[i+1];
	int level;
	
	i=0;
	rewind(in_wind);
	level = z0;
	
	//FILE *in_wind;
	
	//in_wind=fopen("wind.wind", "r");
	
	// read wind file and store the lines to local memory
	while(NULL != fgets(line, MAX_LINE, in_wind)){
		if(line[0] == '#')continue;
		else{
			while(ret=sscanf(line,
			"%lf %lf %lf %lf %lf",
			&wind_height,
			&wind_speed,
			&wind_dir,
			&wind_temp,
			&wind_pres	), ret != 5){}
		}
		// when there is no wind data of 0m
		// create 0m data based on the model that assume
		// 0m/s speed and same direction to the lowest observation
		if(wind_height == 0 && i==0){
			height[i] = wind_height;
			speed[i] = wind_speed;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp;
			atm_pres[i] = wind_pres;
		}else if(i==0){
			height[i] = 0;
			speed[i] = 0;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp;
			atm_pres[i] = wind_pres;
			i++;
			height[i] = wind_height;
			speed[i] = wind_speed;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp;
			atm_pres[i] = wind_pres;
		}else{
			height[i] = wind_height;
			speed[i] = wind_speed;
			dir[i] = wind_dir;
			atm_temp[i] = wind_temp;
			atm_pres[i] = wind_pres;
		}
		//printf("wind_height\t%1.1f\n", height[i]);
		i++;
	}
	
	total=i;
	i=0;
	
	// transfer the wind data to the global memory
	W1 = (WIND *)malloc((total) * sizeof(WIND));
	
	while(i < total){
	    W1[i].day=0;
	    W1[i].hour=0;
	    W1[i].wind_height=height[i];
	    W1[i].wind_speed=speed[i];
	    W1[i].wind_dir=dir[i];
		W1[i].t_atm=atm_temp[i];
		W1[i].p_atm=atm_pres[i];
		i++;
	}
	
	i=0;
	
	/*
	while(i < total){
	    printf("%1.0f, %1.4f, %1.4f\n", W1[i].wind_height, W1[i].wind_speed, W1[i].wind_dir);
		i++;
	}	*/
	return total;
}