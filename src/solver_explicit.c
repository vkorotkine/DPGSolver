// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_explicit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // ToBeModified
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "explicit_VOLUME_info.h"
#include "explicit_FACE_info.h"
#include "finalize_RHS.h"
#include "output_to_paraview.h"
#include "explicit_GradW.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Perform explicit time-stepping using:
 *			1) a 3rd order (S)trong (S)tability (P)reserving (R)unge-(K)utta scheme.
 *			2) a low storage 4th order RK scheme.
 *
 *	Comments:
 *		rk4c is only needed if there is a time-dependent term in the residual (e.g. a time-dependent source term).
 *
 *	Notation:
 *
 *	References:
 *		Carpenter(1994)-Fourth-Order_2N-Storage_Runge-Kutta_Schemes
 *		Gottlieb(2001)-Strong_Stability-Preserving_High-Order_Time_Discretization_Methods (eq. (4.2))
 */

static void correct_What_Bezier(double *WhatB, const double t, const unsigned int NvnS, const double *WAvg,
                                const char type)
{
	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	if (type == 'd') { // (d)ensity
		double *UhatB = malloc(NvnS*Nvar * sizeof *UhatB); // free

		convert_variables(WhatB,UhatB,d,d,NvnS,1,'c','p');

		for (size_t n = 0; n < NvnS; n++) {
			UhatB[n] *= t;
			UhatB[n] += (1.0-t)*WAvg[0];
		}

		convert_variables(UhatB,WhatB,d,d,NvnS,1,'p','c');
		free(UhatB);
	} else if (type == 'p') { // (p)ressure
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NvnS; n++) {
			WhatB[var*NvnS+n] *= t;
			WhatB[var*NvnS+n] += (1.0-t)*WAvg[var];
		}}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

void enforce_positivity_highorder(struct S_VOLUME *VOLUME)
{
	if (VOLUME->type == TET || VOLUME->type == WEDGE || VOLUME->type == PYR)
		return; // Bezier basis not yet implemented. (ToBeModified)
	bool PrintOn = 0;

	if (PrintOn)
		printf("\n\n\n *********** indexg: %d ***********\n\n\n\n",VOLUME->indexg);

	/*
	 *	Purpose:
	 *		Enforce positivity of density and pressure by limitting the high-order solution.
	 *
	 *	Comments:
	 *		See Wang(2012) section 3.1 for details of the procedure.
	 *		The positivity is ensured to hold throughout the entire element by performing the check in the Bezier basis.
	 *		The compensation for the orthonormal basis scaling is defined such that TS[0] multiplied by an array of
	 *		the same constant value should return that value.
	 *		The function currently assumes that the conservative variables are being used.
	 *
	 *	References:
	 *		Wang-Shu(2012)-Robust_High_Order_Discontinuous_Galerkin_Schemes_for_Two-Dimensional_Gaseous_Detonations
	 */

	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	struct S_ELEMENT *ELEMENT = get_ELEMENT_type(VOLUME->type);

	unsigned int P, NvnS;
	double       Volume, *TS, *TS_vB, *TInvS_vB;

	P = VOLUME->P;

	Volume   = ELEMENT->Volume;
	NvnS     = ELEMENT->NvnS[P];
	TS       = ELEMENT->TS[P][P][0];
	TS_vB    = ELEMENT->TS_vB[P][P][0];
	TInvS_vB = ELEMENT->TInvS_vB[P][P][0];

	// *** Density *** ///

	// Find average density
	double rhoAvg, *What, *rho_hat;

	What    = VOLUME->What;
	rho_hat = &What[0];

	// Note compensation for orthonormal basis scaling
	mm_d(CBCM,CBT,CBNT,1,1,NvnS,1.0/sqrt(Volume),0.0,&TS[0],rho_hat,&rhoAvg);

//printf("\n\n\n\n\n\nrhoAvg: % .3e\nWhat:\n",rhoAvg);
//array_print_d(NvnS,Nvar,What,'C');

	if (rhoAvg < EPS_PHYS)
		printf("Error: Average density approaching 0.\n"), EXIT_MSG;

	// Convert to the Bezier basis
	double *WhatB, *rho_hatB;

	WhatB = malloc(NvnS*Nvar * sizeof *WhatB); // free
	mm_CTN_d(NvnS,Nvar,NvnS,TS_vB,What,WhatB);

	rho_hatB = &WhatB[0];

	// Find minimum value
	double rhoMin = 1.0/EPS;
	for (size_t n = 0; n < NvnS; n++) {
		if (rho_hatB[n] < rhoMin)
			rhoMin = rho_hatB[n];
	}

//printf("rhoMin: % .3e\nrho_hatB\n",rhoMin);
//array_print_d(NvnS,1,rho_hatB,'C');

	// Correct rho if necessary (Note: momentum terms are corrected as well)
	if (rhoMin < EPS_PHYS) {
if (PrintOn) {
printf("Volume/rhoAvg: % .3e % .3e\n\n",Volume,rhoAvg);
array_print_d(NvnS,NvnS,TS,'R');
array_print_d(NvnS,1,rho_hatB,'C');
array_print_d(NvnS,Nvar,What,'C');
}
		double t = min(1.0,(rhoAvg-EPS_PHYS)/(rhoAvg-rhoMin));
		correct_What_Bezier(WhatB,t,NvnS,&rhoAvg,'d');

		mm_CTN_d(NvnS,Nvar-1,NvnS,TInvS_vB,WhatB,What);
if (PrintOn) {
array_print_d(NvnS,Nvar,What,'C');
}
	}

	// *** Pressure *** ///
	double *p_hatB;

	p_hatB = malloc(NvnS * sizeof *p_hatB); // free
	compute_pressure(WhatB,p_hatB,d,NvnS,1,'c');

/*
	double *W_vS, *p_vS;

	W_vS = malloc(NvnS*Nvar * sizeof *W_vS); // free
	mm_CTN_d(NvnS,Nvar,NvnS,1.0,ELEMENT->ChiS_vS[P][P][0],What,W_vS);

	p_vS = malloc(NvnS * sizeof *p_vS); // free
	compute_pressure(W_vS,p_vS,d,NvnS,1,'c');

	// Convert to the Bezier basis
	double *p_hatB;

	p_hatB = malloc(NvnS * sizeof *p_hatB); // free
	mm_CTN_d(NvnS,1,NvnS,1.0,ELEMENT->ChiBezInvS_vS[P][P][0],p_vS,p_hatB);
*/

	// Correct W if necessary
	double pMin = 1.0/EPS;
	for (size_t n = 0; n < NvnS; n++) {
		if (p_hatB[n] < pMin)
			pMin = p_hatB[n];
	}
//	free(p_hatB); // Uncomment this: ToBeDeleted

	if (pMin < 0.0) {
if (PrintOn) {
array_print_d(NvnS,NvnS,TS_vB,'R');
printf("\n\n\n\n\nVOLUME: %d\nWhat/WhatB:\n",VOLUME->indexg);
array_print_d(NvnS,Nvar,What,'C');
array_print_d(NvnS,Nvar,WhatB,'C');
printf("pMin: % .4e\np_hatB:\n",pMin);
array_print_d(NvnS,1,p_hatB,'C');
}
		// Compute average pressure
		double *WAvg, pAvg;

		WAvg = malloc(Nvar * sizeof *WAvg); // free

		// Note compensation for orthonormal basis scaling
		mm_d(CBCM,CBT,CBNT,1,Nvar,NvnS,1.0/sqrt(Volume),0.0,&TS[0],What,WAvg);
//		mm_CTN_d(1,Nvar,NvnS,&TS[0],What,WAvg);

		compute_pressure(WAvg,&pAvg,d,1,1,'c');
if (PrintOn) {
printf("pAvg: % .4e\nWAvg:\n",pAvg);
array_print_d(1,Nvar,WAvg,'C');
}
		free(WAvg);

		if (pAvg < 0.0)
			printf("Error: Negative average pressure.\n"), EXIT_MSG;

		//     t = min(1.0,(pAvg-0.0)/(pAvg-pMin));
		double t = pAvg/(pAvg-pMin);
		correct_What_Bezier(WhatB,t,NvnS,WAvg,'p');

		// Correct What
		mm_CTN_d(NvnS,Nvar,NvnS,TInvS_vB,WhatB,What);
if (PrintOn) {
printf("% .3e\nWhatB (updated)\n",t);
array_print_d(NvnS,Nvar,WhatB,'C');
printf("What (corrected):\n");
array_print_d(NvnS,Nvar,What,'C');
}
//if (t < 0.4)
//EXIT_MSG;
//if (pAvg > 2)
//EXIT_MSG;
	}
	free(WhatB);
	free(p_hatB);
}

struct S_timestepping {
	double dt, exit_tol, exit_ratio;
};

static void select_timestepping_parameters(struct S_timestepping *data)
{
	char *PDE      = DB.PDE,
	     *TestCase = DB.TestCase;

	if (strstr(PDE,"NavierStokes")) {
		if (strstr(TestCase,"TaylorCouette")) {
			printf("Using default value for timestepping parameters.\n");
			data->dt         = 1e-4;
			data->exit_tol   = EPS;
			data->exit_ratio = 1.0/EPS;
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(PDE,"Euler")) {
		if (strstr(TestCase,"PeriodicVortex")) {
			printf("Using default value for timestepping parameters.\n");
			data->dt         = 1e-3; // Selected for stability of P3 ML5 (where ML0 has 4 QUADs)
			data->exit_tol   = EPS;
			data->exit_ratio = 1.0/EPS;
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void solver_explicit(void)
{
	// Initialize DB Parameters
	bool         Viscous            = DB.Viscous;
	unsigned int OutputInterval     = DB.OutputInterval,
	             Neq                = DB.Neq,
	             ExplicitSolverType = DB.ExplicitSolverType,
	             Adapt              = DB.Adapt;

	unsigned int PrintTesting = 1;

	double       FinalTime = DB.FinalTime;

	// Standard datatypes
	static double rk4a[5] = { 0.0,              -0.417890474499852, -1.192151694642677, -1.697784692471528,
	                         -1.514183444257156 },
	              rk4b[5] = { 0.149659021999229, 0.379210312999627,  0.822955029386982,  0.699450455949122,
	                          0.153057247968152 };
//	              rk4c[5] = { 0.0,               0.149659021999229,  0.370400957364205,  0.622255763134443,
//                            0.958282130674690 };

	char         *dummyPtr_c[2];
	unsigned int i, iMax, tstep, rk,
	             NvnS;
	double       time, dt, maxRHS0, maxRHS, *RES, *RHS, *What, exit_tol, exit_ratio;

	struct S_timestepping *data_time;

	data_time = malloc(sizeof *data_time); // free

	// silence
	dt = maxRHS0 = 0.0;

	for (i = 0; i < 2; i++)
		dummyPtr_c[i] = malloc(STRLEN_MIN * sizeof *dummyPtr_c[i]); // free

	exit_tol   = 10*EPS;
	exit_ratio = 2e4;

// Need to improve how dt is selected! Likely based on characteristic speeds (see nodalDG code for one possibility).  (ToBeDeleted)
/*
	if (!Adapt) {
		dt = pow(0.5,DB.ML+DB.PGlobal+1);
		dt *= 1e-1;
	} else {
		if (Adapt == ADAPT_P)
			dt = pow(0.5,DB.ML+DB.PMax+1);
		else if (Adapt == ADAPT_H)
			dt = pow(0.5,(DB.ML+DB.LevelsMax)+DB.PGlobal+1);
		else if (Adapt == ADAPT_HP) {
			if (TestDB.Active) {
				char *MeshFile = DB.MeshFile;
				unsigned int ML = TestDB.ML,
							 P  = TestDB.PGlobal;

				exit_tol = 4e-7;
				if (strstr(MeshFile,"Subsonic")) {
					// RK3_SSP, Conforming TRI
					double scaling = 0.0;

					if (P <= 6)
						scaling = 2.5;
					else if (P == 7 || P == 8)
						return;
					else if (P <= 8)
						scaling = 1.25;
					else
						scaling = 0.0;

					if (ML > 0)
						return;

					if      (ML == 0) dt = scaling*5.0e-1;
					else if (ML == 1) dt = scaling*1.0e+0;
					else if (ML == 2) dt = scaling*2.0e+0;
					else if (ML == 3) dt = scaling*4.0e+0;
					else
						printf("Add support (%d).\n",ML), EXIT_MSG;
				} else if (strstr(MeshFile,"Supersonic")) {
					exit_tol   = 3e-15;
					exit_ratio = 2e15;
					if      (ML == 0) dt = 1.0e-2;
					else if (ML == 1) dt = 1.0e-1;
					else
						printf("Add support (%d %d).\n",ML,P), EXIT_MSG;
				}

				OutputInterval = 1e4;
			} else { // Standard
				dt = 5e+0*pow(0.5,DB.ML+DB.PGlobal); //P2
			}
		}
	}
	dt = 1e-3;
	exit_tol   = EPS;
	exit_ratio = 1e15;
*/
	select_timestepping_parameters(data_time);
	dt         = data_time->dt;
	exit_tol   = data_time->exit_tol;
	exit_ratio = data_time->exit_ratio;

	// Compute Mass matrix for uncollocated schemes
	update_VOLUME_Ops();
	update_VOLUME_finalize();

	output_to_paraview("ZTest_Sol_Init");

	tstep = 0; time = 0.0;
	while (time < FinalTime) {
		if (Adapt && tstep)
			mesh_update();

		if (time+dt > FinalTime)
			dt = FinalTime-time;

		if (tstep == 0) {
			for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
				enforce_positivity_highorder(VOLUME);
		}

		switch (ExplicitSolverType) {
		default : // RK3_SSP
			// RES is used to store the initial solution at the beginning of the time step.
			for (rk = 0; rk < 3; rk++) {
				// Compute weak gradients (for viscous terms)
				if (Viscous)
					explicit_GradW();

				// Build the RHS (== -Residual)
				printf("V");  explicit_VOLUME_info();
				printf("F");  explicit_FACE_info();
				printf("F "); maxRHS = finalize_RHS();

				// Update What
				for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
					NvnS = VOLUME->NvnS;

					RES  = VOLUME->RES;
					RHS  = VOLUME->RHS;
					What = VOLUME->What;

					if (rk == 0) {
						for (iMax = Neq*NvnS; iMax--; ) {
							*RES++   = *What;
							*What++ += dt*(*RHS++);
						}
					} else if (rk == 1) {
						for (iMax = Neq*NvnS; iMax--; ) {
							*What = 0.25*(3.0*(*RES++) + *What + dt*(*RHS++));
							What++;
						}
					} else if (rk == 2) {
						for (iMax = Neq*NvnS; iMax--; ) {
							*What = (1.0/3.0)*(*RES++ + 2.0*(*What) + 2.0*dt*(*RHS++));
							What++;
						}
					}
					enforce_positivity_highorder(VOLUME);
				}
			}
			break;
		case RK4_LS:
			for (rk = 0; rk < 5; rk++) {
				// Build the RHS (== -Residual)
				printf("V");  explicit_VOLUME_info();
				printf("F");  explicit_FACE_info();
				printf("F "); maxRHS = finalize_RHS();

				// Update What
				for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
					NvnS = VOLUME->NvnS;

					RES  = VOLUME->RES;
					RHS  = VOLUME->RHS;
					What = VOLUME->What;

					for (iMax = Neq*NvnS; iMax--; ) {
						*RES    *= rk4a[rk];
						*RES    += dt*(*RHS++);
						*What++ += rk4b[rk]*(*RES++);
					}
					enforce_positivity_highorder(VOLUME);
				}
			}
			break;
		case EULER:
			// Build the RHS (== -Residual)
			printf("V");  explicit_VOLUME_info();
			printf("F");  explicit_FACE_info();
			printf("F "); maxRHS = finalize_RHS();

			// Update What
			for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				NvnS = VOLUME->NvnS;

				RHS  = VOLUME->RHS;
				What = VOLUME->What;

				for (iMax = Neq*NvnS; iMax--; )
					*What++ += dt*(*RHS++);
				enforce_positivity_highorder(VOLUME);
			}
		}
		time += dt;

		// Output to paraview
		if (PrintTesting && (tstep % OutputInterval == 0 || tstep < 5)) {
			sprintf(dummyPtr_c[1],"%d",tstep);
			strcpy(dummyPtr_c[0],"SolStart");
			strcat(dummyPtr_c[0],dummyPtr_c[1]);
			output_to_paraview(dummyPtr_c[0]);
		}

		// Display solver progress
		if (!tstep)
			maxRHS0 = maxRHS;

		printf("Complete: % 7.2f%%, tstep: %8d, maxRHS (no MInv): % .3e\n",100*time/FinalTime,tstep,maxRHS);

		// Additional exit conditions
//		if ((maxRHS0/maxRHS > 1e3 || maxRHS < 8e-14) && tstep > 2) {
		if (tstep > 2) {
			if (maxRHS0/maxRHS > exit_ratio) {
				printf("Exiting: maxRHS dropped by % .2e orders.\n",log10(exit_ratio));
				break;
			} else if (maxRHS < exit_tol) {
				printf("Exiting: maxRHS is below % .3e.\n",exit_tol);
				break;
			}
		}

		// hp adaptation
//		if (Adapt)
		if (0&&Adapt)
			adapt_hp();

		tstep++;
	}

	// Output to paraview
	if (TestDB.ML <= 1 || (TestDB.PGlobal == 1) || (TestDB.PGlobal+TestDB.ML) <= 8) {
		char *fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut), // free
		     *string   = malloc(STRLEN_MAX * sizeof *string);   // free

		strcpy(fNameOut,"SolFinal_");
		sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
		                               strcat(fNameOut,DB.MeshType);
		if (DB.Adapt == ADAPT_0) {
			sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
			sprintf(string,"P%d_",DB.PGlobal); strcat(fNameOut,string);
		} else {
			sprintf(string,"_ML%d",TestDB.ML); strcat(fNameOut,string);
			sprintf(string,"P%d_",TestDB.PGlobal); strcat(fNameOut,string);
		}
		output_to_paraview(fNameOut);

		free(fNameOut);
		free(string);
	}

	for (i = 0; i < 2; i++)
		free(dummyPtr_c[i]);
	free(data_time);
}
