// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_explicit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // ToBeModified
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "explicit_VOLUME_info.h"
#include "explicit_FACE_info.h"
#include "finalize_RHS.h"
#include "output_to_paraview.h"
#include "explicit_GradW.h"

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
	static double rk4a[5] = { 0.0,              -0.417890474499852, -1.192151694642677, -1.697784692471528, -1.514183444257156 },
	              rk4b[5] = { 0.149659021999229, 0.379210312999627,  0.822955029386982,  0.699450455949122,  0.153057247968152 };
//	              rk4c[5] = { 0.0,               0.149659021999229,  0.370400957364205,  0.622255763134443,  0.958282130674690 };

	char         *dummyPtr_c[2];
	unsigned int i, iMax, tstep, rk,
	             NvnS;
	double       time, dt, maxRHS0, maxRHS, *RES, *RHS, *What;

	struct S_VOLUME *VOLUME;

	// silence
	dt = maxRHS0 = 0.0;

	for (i = 0; i < 2; i++)
		dummyPtr_c[i] = malloc(STRLEN_MIN * sizeof *dummyPtr_c[i]); // free

// Need to improve how dt is selected! Likely based on characteristic speeds (see nodalDG code for one possibility).  (ToBeDeleted)
	if (!Adapt) {
		dt = pow(0.5,DB.ML+DB.PGlobal+1);
		dt *= 1e-1;
	} else {
		if (Adapt == ADAPT_P)
			dt = pow(0.5,DB.ML+DB.PMax+1);
		else if (Adapt == ADAPT_H)
			dt = pow(0.5,(DB.ML+DB.LevelsMax)+DB.PGlobal+1);
		else if (Adapt == ADAPT_HP)
//			dt = 1e4*pow(0.5,max(DB.ML,DB.LevelsMax)+DB.PMax+1);
			dt = 1e+1*pow(0.5,DB.ML+DB.PGlobal);
			OutputInterval = 2e3;
	}
	printf("%d %d\n",DB.ML,DB.PGlobal);

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
				}
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
		if ((maxRHS0/maxRHS > 1e10 || maxRHS < 8e-14) && tstep > 2) {
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
