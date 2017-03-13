// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_explicit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // ToBeModified
#include <string.h>
#include <unistd.h>

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

static void enforce_positivity_highorder(struct S_VOLUME *VOLUME)
{
//	return;
	/*
	 *	Purpose:
	 *		Enforce positivity of density and pressure by limitting the high-order solution.
	 *
	 *	Comments:
	 *		See Wang(2012) section 3.1 for details of the procedure.
	 *		The positivity is ensured to hold throughout the entire element by performing the check in the Bezier basis.
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
	mm_d(CBCM,CBT,CBNT,1,1,NvnS,sqrt(Volume),0.0,&TS[0],rho_hat,&rhoAvg);

//printf("\n\n\n\n\n\nrhoAvg: % .3e\nWhat:\n",rhoAvg);
//array_print_d(NvnS,Nvar,What,'C');

	if (rhoAvg < EPS_PHYS)
		printf("Error: Average density approaching 0.\n"), EXIT_MSG;

	// Convert to the Bezier basis
	double *rho_hatB;

	rho_hatB = malloc(NvnS * sizeof *rho_hatB); // free
	mm_CTN_d(NvnS,1,NvnS,TS_vB,rho_hat,rho_hatB);

	// Find minimum value
	double rhoMin = 1.0/EPS;
	for (size_t n = 0; n < NvnS; n++) {
		if (rho_hatB[n] < rhoMin)
			rhoMin = rho_hatB[n];
	}

//printf("rhoMin: % .3e\nrho_hatB\n",rhoMin);
//array_print_d(NvnS,1,rho_hatB,'C');

	// Correct rho if necessary
	if (rhoMin < EPS_PHYS) {
printf("% .3e\n",rhoAvg);
array_print_d(NvnS,1,rho_hatB,'C');
array_print_d(NvnS,Nvar,What,'C');
		double t = min(1.0,(rhoAvg-EPS_PHYS)/(rhoAvg-rhoMin));
		for (size_t n = 0; n < NvnS; n++) {
			rho_hatB[n] *= t;
			rho_hatB[n] += (1.0-t)*rhoAvg;
		}

		// Correct rho in What
		mm_CTN_d(NvnS,1,NvnS,TInvS_vB,rho_hatB,rho_hat);
array_print_d(NvnS,Nvar,What,'C');
	}
	free(rho_hatB);

	// *** Pressure *** ///
	double *WhatB, *p_hatB;

	WhatB = malloc(NvnS*Nvar * sizeof *WhatB); // free
	mm_CTN_d(NvnS,Nvar,NvnS,TS_vB,What,WhatB);

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
array_print_d(NvnS,NvnS,TS_vB,'R');
printf("\n\n\n\n\nVOLUME: %d\nWhat/WhatB:\n",VOLUME->indexg);
array_print_d(NvnS,Nvar,What,'C');
array_print_d(NvnS,Nvar,WhatB,'C');
printf("pMin: % .4e\np_hatB:\n",pMin);
array_print_d(NvnS,1,p_hatB,'C');
		// Compute average pressure
		double *WAvg, pAvg;

		WAvg = malloc(Nvar * sizeof *WAvg); // free

		// Note compensation for orthonormal basis scaling
		mm_d(CBCM,CBT,CBNT,1,Nvar,NvnS,sqrt(Volume),0.0,&TS[0],What,WAvg);

		compute_pressure(WAvg,&pAvg,d,1,1,'c');
printf("pAvg: % .4e\nWAvg:\n",pAvg);
array_print_d(1,Nvar,WAvg,'C');
		free(WAvg);

		if (pAvg < 0.0)
			printf("Error: Negative average pressure.\n"), EXIT_MSG;

		//     t = min(1.0,(pAvg-0.0)/(pAvg-pMin));
		double t = pAvg/(pAvg-pMin);

		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NvnS; n++) {
			WhatB[var*NvnS+n] *= t;
			WhatB[var*NvnS+n] += (1.0-t)*WAvg[var];
		}}

printf("% .3e\nWhatB (updated)\n",t);
array_print_d(NvnS,Nvar,WhatB,'C');

p_hatB = malloc(NvnS * sizeof *p_hatB); // free
compute_pressure(WhatB,p_hatB,d,NvnS,1,'c');
array_print_d(NvnS,1,p_hatB,'C');

		// Correct What
		mm_CTN_d(NvnS,Nvar,NvnS,TInvS_vB,WhatB,What);
printf("What (corrected):\n");
array_print_d(NvnS,Nvar,What,'C');
if (t < 0.5)
EXIT_MSG;
if (pAvg > 2)
EXIT_MSG;
	}
	free(WhatB);
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

	struct S_VOLUME *VOLUME;

	// silence
	dt = maxRHS0 = 0.0;

	for (i = 0; i < 2; i++)
		dummyPtr_c[i] = malloc(STRLEN_MIN * sizeof *dummyPtr_c[i]); // free

	exit_tol   = 10*EPS;
	exit_ratio = 1e4;

// Need to improve how dt is selected! Likely based on characteristic speeds (see nodalDG code for one possibility).  (ToBeDeleted)
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
				char         MType[20];
				unsigned int ML = TestDB.ML,
							 P  = TestDB.PGlobal;

//				strcpy(MType,"supersonic");
				strcpy(MType,"subsonic");

				exit_tol = 4e-8;
				if (strstr(MType,"subsonic")) {
					// RK3_SSP, Conforming TRI
					if      (ML == 0) dt = 5.0e-1;
					else if (ML == 1) dt = 2.0e+0;
					else if (ML == 2) dt = 4.0e+0;
					else if (ML == 3) dt = 8.0e+0;
					else
						printf("Add support (%d).\n",ML), EXIT_MSG;
				} else if (strstr(MType,"supersonic")) {
					if      (ML == 0) dt = 1.0e-3;
					else
						printf("Add support (%d).\n",ML), EXIT_MSG;
				}

				if        (P == 1) {
				} else if (P == 2) {
				} else if (P == 3) {
				} else if (P == 4) {
				} else if (P == 5) {
				} else if (P == 6) {
				} else if (P == 7) {
				} else if (P == 8) {
				} else {
					printf("Add support.\n"), EXIT_MSG;
				}
				OutputInterval = 1e4;
			} else { // Standard
				dt = 5e+0*pow(0.5,DB.ML+DB.PGlobal); //P2
			}
		}
	}

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
				for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
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
				for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
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
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
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
		if ((maxRHS0/maxRHS > exit_ratio || maxRHS < exit_tol) && tstep > 2) {
			printf("Exiting: maxRHS dropped by 10 orders or is below 8e-14.\n");
			break;
		}

		// hp adaptation
//		if (Adapt)
		if (0&&Adapt)
			adapt_hp();

		tstep++;
	}

	for (i = 0; i < 2; i++)
		free(dummyPtr_c[i]);
}
