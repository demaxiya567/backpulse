#include "udf.h"
#define Ax 55e-13 /*permeability in x direction */
#define Ay 55e-13 /*permeability in y direction */
#define PAR_DEN 2700 /*particle real density*/
#define CELL_POROUSITY 0.524 /*porousity of dust cake*/
#define DEPOSIT_FACE_ID 17 /*dust deposit face thread id*/
#define AIR_VISOCOSITY 1.83e-5
#define FACE_VELOCITY 0.02
#define ADJUST_DYNAMIC_SHAPE_FACTOR 1.25
#define MOLE_FREE_PATH 6.9e-8
#define FLAG_USE_FACE_VELOCITY 1
/*for data output*/
real source_x = 0;
real source_y = 0;

/*count the cunningham correction factor in a certain cell*/
real CUNNINGHAN_CORRECTION(real particle_diam, real molecular_free_path);

/*count the adjust dynamic shape factor factor in a certain cell*/
real ADJUST_DYN_SHA_FAC(real vm_diam, real vm_cunn, real sm_diam, real sm_cunn);

/*count the total dust loading mass for a certain cell*/
real TOTAL_MASS_C(cell_t cp, Thread *tp, real timestepp);

/*count the total dust loading face area for a certain cell*/
int FACE_AREA_C(real area_f[], cell_t cp, Thread *tp, int deposit_face_ID);

/*function to count the dustcake in xyz direction at a certain cell, not used now*/
/*real DUSTCAKE_THICKNESS(real c_porousity, real p_real_density, int deposit_face_ID, cell_t cp, Thread *tp);*/

/*count the total dust loading mass for a certain cell*/
int FACE_DUST_LOAD_C(real face_dust_load[], cell_t cp, Thread *tp);

/*count the total porousity for a certain cell*/
real TOTALL_POROUSITY(real face_velocity, real face_dust_load);

/*count the total pressure drop on the dintust cake  for a certain cell
 theroy detail showed on paper "Compression properties of dust cake of fine fly ashes from a fluidized bed coal combustor on a ceramic filter" */
real PRESSURE_DROP_SOURCE(real c_porousity, real vcsity,
		real ad_dynamic_shape_fac, real particle_density,
		real geometric_mean_diameter, real geometric_standard_deviation,
		real face_velocity, real face_dust_load);

DEFINE_SOURCE(xmom_source, c, t, dS, eqn) {
	/*reserved for a classic use*/
	/*	 real mu_coeiff;
	 real thickness[ND_ND];
	 real x[ND_ND];
	 real con, source;
	 real effu;
	 thickness = dustcake_thickness(Cell_porosity, PAR_DEN, DE_FACE_ID, c, t);
	 effu = C_MU_EFF(c, t);
	 con = effu / Ax*thickness[0];*/

	real porous;
	real air_vc;/*Viscosity*/
	real adsf;/*adjust_dynamic_shape_factor*/
	real pd;/*particle density*/
	real gmd;/*geometric_mean_diameter*/
	real gsd;/*geometric_standard_deviation*/
	real fvx;/*face_velocity*/
	real flx;/*face_dust_load*/
	real source; /*pressure drop on x direction*/
	real face_load_xy[ND_ND];
	real x[ND_ND];

	if (FLAG_USE_FACE_VELOCITY)
		fvx = FACE_VELOCITY;
	else
		fvx = C_U(c, t);

	C_CENTROID(x, c, t);
	/*the area you want to apply the source*/
	if (sqrt(x[0] * x[0] + x[1] * x[1]) > 0.03) {
		FACE_DUST_LOAD_C(face_load_xy, c, t);
		flx = face_load_xy[0];
		if (face_load_xy[0] != 0) {
			Message("your face_load in x now is %g \n", face_load_xy[0]);
		}
		porous = TOTALL_POROUSITY(fvx, flx);
		air_vc = AIR_VISOCOSITY;
		adsf = ADJUST_DYNAMIC_SHAPE_FACTOR;
		gmd = 2.5e-6;
		gsd = 1.5;
		pd = PAR_DEN;
		Message(
				"your porous air_vc adsf pd gmd gsd fvx flx is %g %g %g %g %g %g %g %g \n",
				porous, air_vc, adsf, pd, gmd, gsd, fvx, flx);
		source = -PRESSURE_DROP_SOURCE(porous, air_vc, adsf, pd, gmd, gsd, fvx,
				flx);
		source_x = source;
		C_UDMI(c,t,0) = source;
		Message("your x pressure now is %g\n", source);
		/*
		 if (c>=100 && c<=120)
		 {
		 printf('the x face thickness is %g\n ', thickness[0]);
		 }*/

		dS[eqn] = source / fvx;
		return source;
	} else {
		dS[eqn] = 0;
		C_UDMI(c,t,0) = source;
		return source = 0;
	}
}

DEFINE_SOURCE(ymom_source, c, t, dS, eqn) {

	real porous;
	real air_vc;
	real adsf;
	real pd;
	real gmd;
	real gsd;
	real fvy;
	real fly;
	real source;
	real face_load_xyy[ND_ND];
	real y[ND_ND];
	if (FLAG_USE_FACE_VELOCITY)
		fvy = FACE_VELOCITY;
	else
		fvy = C_V(c, t);

	C_CENTROID(y, c, t);	C_CENTROID(y, c, t);
	if (sqrt(y[0] * y[0] + y[1] * y[1]) > 0.03) {

		FACE_DUST_LOAD_C(face_load_xyy, c, t);
		fly = face_load_xyy[1];
		if (face_load_xyy[1] != 0) {
			Message("your face load in y now is %g \n", face_load_xyy[1]);
		}
		porous = TOTALL_POROUSITY(fvy, fly);
		air_vc = AIR_VISOCOSITY;
		adsf = ADJUST_DYNAMIC_SHAPE_FACTOR;
		gmd = 2.5e-6;
		gsd = 1.5;
		pd = PAR_DEN;
		Message(
				"your porous air_vc adsf pd gmd gsd fvx flx is %g %g %g %g %g %g %g %g \n",
				porous, air_vc, adsf, pd, gmd, gsd, fvy, fly);
		source = -PRESSURE_DROP_SOURCE(porous, air_vc, adsf, pd, gmd, gsd, fvy,
				fly);
		source_y = source;
		Message("your y pressure now is %g\n", source);
		/*
		 if (c >= 100 && c <= 120)
		 {
		 printf('the y face thickness is %g\n ', thickneseal FACE_AREA_C(real *area_f, cell_t cp,Thread *tp,int deposit_face_ID)s[1]);
		 }
		 source = -con*C_V(c, t);*/
		C_UDMI(c,t,1) = source;
		dS[eqn] = source / fvy;
		return source;
	} else {
		C_UDMI(c,t,1) = source;
		dS[eqn] = 0;
		return source = 0;
	}
}

DEFINE_OUTPUT_PARAMETER(pressure_x,n,parlist) {
	real prex = source_x;
	return prex;
}

DEFINE_OUTPUT_PARAMETER(pressure_y,n,parlist) {
	real prey = source_y;
	return prey;
}

real CUNNINGHAN_CORRECTION(real particle_diam, real molecular_free_path) {
	if (particle_diam > 2 * molecular_free_path) {
		return 1 + 2.468 * molecular_free_path / particle_diam;
	} else {
		return 1 + 3.294 * molecular_free_path / particle_diam;
	}

}

real ADJUST_DYN_SHA_FAC(real vm_diam, real vm_cunn, real sm_diam, real sm_cunn) {

	return CUNNINGHAN_CORRECTION(vm_diam, MOLE_FREE_PATH)
			/ CUNNINGHAN_CORRECTION(sm_diam, MOLE_FREE_PATH)
			* (vm_diam / sm_diam) * (vm_diam / sm_diam);
}

/*real DUSTCAKE_THICKNESS(real c_porousity, real p_real_density,int deposit_face_ID,cell_t cp, Thread *tp)
 {int
 real total_mass_c;
 real real_v_c;
 real v_dustcake_c;
 real deposit_face_area_c;
 real dust_cake_height_c[ND_ND];	C_CENTROID(y, c, t);
 real A[ND_ND];
 int n;
 total_mass_c = C_DPMS_CONCENTRATION(cp, tp)*C_VOLUME(cp, tp)*CURRENT_TIMESTEP;
 real_v_c = total_mass_c / p_real_density;
 v_dustcake_c = real_v_c / c_porousity;
 c_face_loop(cp, tp, n)
 {
 if (THREAD_ID(C_FACE_THREAD(cp, tp, n)) == deposit_face_ID)
 {
 F_AREA(A, C_FACE(cp, tp, n), C_FACE_THREAD(cp, tp, n));
 }
 }

 dust_cake_height_c[0] = v_dustcake_c / A[0];
 dust_cake_height_c[1] = v_dustcake_c / A[1];
 dust_cake_height_c[2] = 0;
 return dust_cake_height_c;
 }*/

real TOTAL_MASS_C(cell_t cp, Thread *tp, real timestepp) {
	if (C_DPMS_CONCENTRATION(cp,tp) != 0)
		Message("your pcon C_V times step is %g %g %g \n",
				C_DPMS_CONCENTRATION(cp, tp), C_VOLUME(cp, tp),
				CURRENT_TIMESTEP);
	return C_DPMS_CONCENTRATION(cp, tp) * C_VOLUME(cp, tp) * timestepp;
}

int FACE_AREA_C(real area_f[], cell_t cp, Thread *tp, int deposit_face_ID) {
	int flag_flowdirection;
	int n;
	real x[ND_ND];
	flag_flowdirection = -1; /*flow flag*/
	if (flag_flowdirection == -1) {
		c_face_loop(cp, tp, n)
		{
			if (THREAD_ID(C_FACE_THREAD(cp, tp, n)) == deposit_face_ID) {
				F_AREA(area_f, C_FACE(cp, tp, n), C_FACE_THREAD(cp, tp, n));
			}
		}

	}

	return 0;
}

/*count the face dust load in a cell mainly in x y direction*/
int FACE_DUST_LOAD_C(real face_dust_load[], cell_t cp, Thread *tp) {
	int ii;
	real face_of_area[ND_ND];
	real total_mass_c;
	real face_load[ND_ND];
	real f_area_t;
	total_mass_c = TOTAL_MASS_C(cp, tp, CURRENT_TIMESTEP);
	if (total_mass_c != 0)
		Message("your total mass is %g\n", total_mass_c);
	FACE_AREA_C(face_of_area, cp, tp, DEPOSIT_FACE_ID);
	f_area_t = NV_MAG(face_of_area);
	Message("your face_dust_load x y z t is %g %g %g %g\n", face_of_area[0],
			face_of_area[1], face_of_area[2], f_area_t);
	/* if your case is 2D please change your dimension*/
#if RP_3D
	for (ii=0; ii<2; ii++)
	{
		Message("you face area is %g\n",fabs(face_of_area[ii]));
		if(fabs(face_of_area[ii])>1e-20) {
			face_dust_load[ii]= total_mass_c / fabs(face_of_area[ii]);
		}
		else
		{
			Message("please keep your face_area have a non-zero value \n");
		}
		Message("you enter here\n");
		Message("your face load is %g \n", face_dust_load[ii]);
	}
#endif

	return 0;
}

/*
 (theroy detail showed on paper "Compression properties of dust cake of fine fly ashes from a
 fluidized bed coal combustor on a ceramic filter")*/
real TOTALL_POROUSITY(real face_velocity, real face_dust_load) {

	return 1
			- (0.88
					- 0.37 * pow(face_velocity, 0.36)
							* pow(face_dust_load, 0.25))
					* pow(face_velocity, 0.36) * pow(face_dust_load, 0.25);

}

/*your pressure drop function*/
real PRESSURE_DROP_SOURCE(real c_porousity, real vcsity,
		real ad_dynamic_shape_fac, real particle_density,
		real geometric_mean_diameter, real geometric_standard_deviation,
		real face_velocity, real face_dust_load) {
	real acoef;
	real bcoef;
	real ccoef;
	real tpc_drop;
	ccoef = 4 * log(geometric_standard_deviation)
			* log(geometric_standard_deviation);
	Message("your ccoef is %g\n", ccoef);
	acoef = (2970 * vcsity * ad_dynamic_shape_fac * (1 - c_porousity))
			/ pow(c_porousity, 4);
	Message("your acoef is %g\n", acoef);
	bcoef = particle_density * geometric_mean_diameter * geometric_mean_diameter
			* exp(ccoef);
	Message("your bcoef is %g\n", bcoef);
	tpc_drop = acoef * pow(bcoef, -1) * face_velocity * face_dust_load;
	Message("your drop is %g\n", tpc_drop);
	return tpc_drop;
}
