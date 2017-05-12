#include "udf.h"
#define Ax 55e-13			 /*permeability in x direction */
#define Ay 55e-13			 /*permeability in y direction */
#define PAR_DEN 2700		 /*particle real density*/
#define CELL_POROUSITY 0.524 /*porousity of dust cake*/
#define DEPOSIT_FACE_ID 32   /*dust deposit face thread id*/
#define DUST_LAYER_ID 34
#define USER_DOMAIN_ID 1
#define AIR_VISOCOSITY 1.83e-5
#define FACE_VELOCITY 0.02
#define ADJUST_DYNAMIC_SHAPE_FACTOR 1.25
#define MOLE_FREE_PATH 6.9e-8
#define FLAG_USE_FACE_VELOCITY 1
#define GEOMETRIC_STANDARD_DEVIATION 1.6
#define GEOMETRIC_MEAN_DIAMETER 2.5e-6
#define MAX_OUTPUT_UDF_FILENAME 25
#define MAX_OUTPUT_UDF_FILEPATH 250
#define FLAG_OUT_PUT 0
#define UDF_OUT_FILE_PATH "/home/liukan12/fluentudf/udf/useforc/udfoutputfile/"
#define UDF_OUT_NAME "dustoutfile"
#define SAVE_BY_TIMER 0
#define EVERY_N_TIME 10
#define EVERY_N_INTER 20
#define FLAG_TRAN 0
/*for data output*/
real source_x = 0;
real source_y = 0;

FILE *UDF_OUTPUT_FILE(char *fpath, char *fname, int saveflag, long int saventime, int saveniter);

/*count the cunningham correction factor in a certain cell*/
real CUNNINGHAN_CORRECTION(real particle_diam, real molecular_free_path);

/*count the adjust dynamic shape factor factor in a certain cell*/
real ADJUST_DYN_SHA_FAC(real vm_diam, real vm_cunn, real sm_diam, real sm_cunn);

/*count the total dust loading mass for a certain cell*/
real TOTAL_MASS_C(cell_t cp, Thread *tp, real timestepp);

/*count the total dust loading face area for a certain cell*/
int FACE_AREA_C(real area_f[], cell_t cp, Thread *tp, int deposit_face_ID);

/*count the total dust loading mass for a certain cell*/
int FACE_DUST_LOAD_C(real face_dust_load[], real face_load[], cell_t cp,
					 Thread *tp);

/*count the total porousity for a certain cell*/
real TOTALL_POROUSITY(real face_velocity, real face_dust_load);

/*count the total pressure drop on the dintust cake  for a certain cell
 theroy detail showed on paper "Compression properties of dust cake of fine fly ashes from a fluidized bed coal combustor on a ceramic filter" */
real PRESSURE_DROP_SOURCE(real c_porousity, real vcsity,
						  real ad_dynamic_shape_fac, real particle_density,
						  real geometric_mean_diameter, real geometric_standard_deviation,
						  real face_velocity, real face_dust_load);

DEFINE_SOURCE(xmom_source, c, t, dS, eqn)
{

	real porous;
	real air_vc; /*Viscosity*/
	real adsf;   /*adjust_dynamic_shape_factor*/
	real pd;	 /*particle density*/
	real gmd;	/*geometric_mean_diameter*/
	real gsd;	/*geometric_standard_deviation*/
	real fvx;	/*face_velocity*/
	real flx;	/*face_dust_load*/
	real source; /*momentum source on x direction*/
	real face_load_xy[ND_ND];
	real face_deposit_area[ND_ND];
	real face_deposit_area_x;
	real delt_P_x;

	if (FLAG_USE_FACE_VELOCITY && fabs(C_U(c, t)) > FACE_VELOCITY)
		fvx = FACE_VELOCITY;
	else
		fvx = fabs(C_U(c, t));

	FACE_DUST_LOAD_C(face_load_xy, face_deposit_area, c, t);
	flx = face_load_xy[0];
	face_deposit_area_x = face_deposit_area[0];
	porous = TOTALL_POROUSITY(fvx, flx);
	air_vc = AIR_VISOCOSITY;
	adsf = ADJUST_DYNAMIC_SHAPE_FACTOR;
	gmd = GEOMETRIC_MEAN_DIAMETER;
	gsd = GEOMETRIC_STANDARD_DEVIATION;
	pd = PAR_DEN;

	delt_P_x = -PRESSURE_DROP_SOURCE(porous, air_vc, adsf, pd, gmd, gsd, fvx,
									 flx) *
			   C_U(c, t) / fabs(C_U(c, t));
	source = delt_P_x * face_deposit_area_x / C_VOLUME(c, t);
	source_x = source;
	printf(
		"your porous air_vc adsf pd gmd gsd fvx flx source is %g %g %g %g %g %g %g %g %g \n",
		porous, air_vc, adsf, pd, gmd, gsd, fvx, flx, source);
	if (FLAG_OUT_PUT)
	{
		FILE *udfp = NULL;
		char *ufpath = (char *)malloc(MAX_OUTPUT_UDF_FILEPATH * sizeof(char));
		char *ufname = (char *)malloc(MAX_OUTPUT_UDF_FILENAME * sizeof(char));
		strcpy(ufpath, UDF_OUT_FILE_PATH);
		strcpy(ufname, UDF_OUT_NAME);

		udfp = UDF_OUTPUT_FILE(ufpath, ufname, SAVE_BY_TIMER, EVERY_N_TIME, EVERY_N_INTER);
		fprintf(udfp,
				"your porous air_vc adsf pd gmd gsd fvx flx source is %g %g %g %g %g %g %g %g %g \n",
				porous, air_vc, adsf, pd, gmd, gsd, fvx, flx, source);
		fclose(udfp);
		/*		free(udfp);
		free(ufpath);0
		free(ufname);*/
	}
	source_x = source;
	C_UDMI(c, t, 0) = source;
	dS[eqn] = source / fvx;
	return source;
}

DEFINE_SOURCE(ymom_source, c, t, dS, eqn)
{

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
	real face_deposit_areay[ND_ND];
	1 real face_deposit_area_y;
	real delt_P_y;

	if (FLAG_USE_FACE_VELOCITY && fabs(C_V(c, t)) > FACE_VELOCITY)
		fvy = FACE_VELOCITY;
	else
		fvy = fabs(C_V(c, t));

	FACE_DUST_LOAD_C(face_load_xyy, face_deposit_areay, c, t);
	fly = face_load_xyy[1];
	face_deposit_area_y = face_deposit_areay[1];
	porous = TOTALL_POROUSITY(fvy, fly);
	air_vc = AIR_VISOCOSITY;
	adsf = ADJUST_DYNAMIC_SHAPE_FACTOR;
	gmd = GEOMETRIC_MEAN_DIAMETER;
	gsd = GEOMETRIC_STANDARD_DEVIATION;
	pd = PAR_DEN;
	delt_P_y = -PRESSURE_DROP_SOURCE(porous, air_vc, adsf, pd, gmd, gsd, fvy,
									 fly) *
			   C_V(c, t) / fabs(C_V(c, t));
	source = delt_P_y * face_deposit_area_y / C_VOLUME(c, t);
	/*	Message(fullname
	 "your porous air_vc adsf pd gmd gsd fvy fly source is %g %g %g %g %g %g %g %g %g \n",
	 porous, air_vc, adsf, pd, gmd, gsd, fvy, fly, source);*/
	if (FLAG20_OUT_PUT)
	{
		FILE *udfp1 = NULL;
		char *ufpath1 = (char *)malloc(MAX_OUTPUT_UDF_FILEPATH * sizeof(char));
		char *ufname1 = (char *)malloc(MAX_OUTPUT_UDF_FILENAME * sizeof(char));
		strcpy(ufpath1, UDF_OUT_FILE_PATH);
		strcpy(ufname1, UDF_OUT_NAME);

		udfp1 = UDF_OUTPUT_FILE(ufpath1, ufname1, SAVE_BY_TIMER, EVERY_N_TIME, EVERY_N_INTER);
		fprintf(udfp1,
				"your porous air_vc adsf pd gmd gsd fvy fly source is %g %g %g %g %g %g %g %g %g \n",
				porous, air_vc, adsf, pd, gmd, gsd, fvy, fly, source);
		fclose(udfp1);
		/*		free(udfp1);
		free(ufpath1);
		free(ufname1);*/
	}

	source_y = source;
	C_UDMI(c, t, 1) = source;
	dS[eqn] = source / fvy;
	return source;
}

DEFINE_OUTPUT_PARAMETER(pressure_x, n, parlist)
{
	real prex = source_x;
	return prex;
}

DEFINE_OUTPUT_PARAMETER(pressure_y, n, parlist)
{
	real prey = source_y;
	return prey;
}

real CUNNINGHAN_CORRECTION(real particle_diam, real molecular_free_path)
{
	if (particle_diam > 2 * molecular_free_path)
	{
		return (1 + 2.468 * molecular_free_path / particle_diam);
	}
	else
	{
		return (1 + 3.294 * molecular_free_path / particle_diam);
	}
}

real ADJUST_DYN_SHA_FAC(real vm_diam, real vm_cunn, real sm_diam, real sm_cunn)
{

	return (CUNNINGHAN_CORRECTION(vm_diam, MOLE_FREE_PATH) / CUNNINGHAN_CORRECTION(sm_diam, MOLE_FREE_PATH) * (vm_diam / sm_diam) * (vm_diam / sm_diam));
}

real TOTAL_MASS_C(cell_t cp, Thread *tp, real timestepp)
{
	/*	if (FLAG_TRAN)*/
	return (C_DPMS_CONCENTRATION(cp, tp) * C_VOLUME(cp, tp));
	/*	else
	C_UDMI(cp,tp,2)+=C_DPMS_CONCENTRATION(cp,tp);
	Message("your dust load is %g \n",C_DPMS_CONCENTRATION(cp,tp));
	return (C_UDMI(cp,tp,2) * C_VOLUME(cp, tp));*/
}

int FACE_AREA_C(real area_f[], cell_t cp, Thread *tp, int deposit_face_ID)
{
	int flag_flowdirection;
	int n;
	flag_flowdirection = -1; /*flow flag*/
	/*#if RP_NODE*/
	if (flag_flowdirection == -1)
	{
		c_face_loop(cp, tp, n)
		{
			if (THREAD_ID(C_FACE_THREAD(cp, tp, n)) == deposit_face_ID)
			{
				F_AREA(area_f, C_FACE(cp, tp, n), C_FACE_THREAD(cp, tp, n));
				/*Message("your face area is %g %g %g\n",area_f[0],area_f[1],area_f[2]);*/
			}
		}
	}
	/*#endif*/
	return 0;
}

/*count the face dust load in a cell mainly in x y direction*/
int FACE_DUST_LOAD_C(real face_dust_load[], real face_load[], cell_t cp,
					 Thread *tp)
{
	int ii;
	real face_of_area[ND_ND];
	real total_mass_c;
	real f_area_t;
	total_mass_c = TOTAL_MASS_C(cp, tp, CURRENT_TIMESTEP);
	FACE_AREA_C(face_of_area, cp, tp, DEPOSIT_FACE_ID);
	f_area_t = NV_MAG(face_of_area);
	face_load[0] = fabs(face_of_area[0]) / sqrt(
											   face_of_area[0] * face_of_area[0] + face_of_area[1] * face_of_area[1]) *
				   f_area_t;
	face_load[1] = fabs(face_of_area[1]) / sqrt(
											   face_of_area[0] * face_of_area[0] + face_of_area[1] * face_of_area[1]) *
				   f_area_t;
/*	Message("your face_dust_load x y z t is %g %g %g %g\n", face_of_area[0],
	 face_of_area[1], face_of_area[2], f_area_t);*/
/* if your case is 2D please change your dimension*/
/*#if RP_NODE*/
#if RP_3D
	for (ii = 0; ii < 2; ii++)
	{
		/*Message("you face area is %g\n",fabs(face_of_area[ii]));*/
		if (fabs(face_load[ii]) > 1e-20)
		{
			face_dust_load[ii] = total_mass_c / face_load[ii];
			/*Message("you face load i "editor.fontSize": 14,s %g\n",face_dust_load[ii]);*/
		}
		else
		{
			face_dust_load[ii] = 0;
		}
	}
#endif
	/*#endif*/
	return 0;
}

/*
 (theroy detail showed on paper "Compression properties of dust cake of fine fly ashes from a
 fluidized bed coal combustor on a ceramic filter")*/
real TOTALL_POROUSITY(real face_velocity, real face_dust_load)
{

	return (1 - (0.88 - 0.37 * pow(fabs(face_velocity), 0.36) * pow(fabs(face_dust_load), 0.25)) * pow(fabs(face_velocity), 0.36) * pow(fabs(face_dust_load), 0.25));
}

/*your pressure drop function*/
real PRESSURE_DROP_SOURCE(real c_porousity, real vcsity,
						  real ad_dynamic_shape_fac, real particle_density,
						  real geometric_mean_diameter, real geometric_standard_deviation,
						  real face_velocity, real face_dust_load)
{
	real acoef;
	real bcoef;
	real ccoef;
	real tpc_drop;
	ccoef = 4 * log(geometric_standard_deviation) * log(geometric_standard_deviation);
	/*Message("your ccoef is %g\n", ccoef);*/
	acoef = (2970 * vcsity * ad_dynamic_shape_fac * (1 - c_porousity)) / pow(c_porousity, 4);
	/*Message("your acoef is %g\n", acoef);*/
	bcoef = particle_density * geometric_mean_diameter * geometric_mean_diameter * exp(ccoef);
	/*Message("your bcoef is %g\n", bcoef);*/
	tpc_drop = acoef * pow(bcoef, -1) * fabs(face_velocity) * fabs(face_dust_load);
	/*Message("your drop is %g\n", tpc_drop);*/
	return tpc_drop;
}

FILE *UDF_OUTPUT_FILE(char *fpath, char *fname, int saveflag, long int saventime, int saveniter)
{

	FILE *fp;
	char *fullname = (char *)malloc((MAX_OUTPUT_UDF_FILEPATH + MAX_OUTPUT_UDF_FILENAME) * sizeof(char));
	char *fappendfix = (char *)malloc(sizeof(long int));
	if (saveflag == 1)
	{
		if (N_TIME % saventime == 0)
		{
			sprintf(fappendfix, 0 "-%ld", N_TIME);
			strcat(fname, fappendfix);
			Message("your time step now is %d \n", N_ITER);
		}
	}
	else
	{
		if (N_ITER % saveniter == 0)
		{
			sprintf(fappendfix, "-%d", N_ITER);
			strcat(fname, fappendfix);
			Message("your iteration now is %d \n", N_ITER);
		}
	}
	strcpy(fullname, fpath);
	strcat(fullname, fname);
	fp = fopen(fullname, "a+");
	return fp;
}