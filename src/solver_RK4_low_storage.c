// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

// ToBeDeleted
#include <math.h>

/*
 *	Purpose:
 *		Perform time-stepping using a low storage 4th order (R)unge-(K)utta scheme.
 *
 *	Comments:
 *		rk4c is only needed if there is a time-dependent term in the residual (e.g. a time-dependent source term).
 *
 *	Notation:
 *
 *	References:
 *		Carpenter(1994)-Fourth-Order_2N-Storage_Runge-Kutta_Schemes
 */

void solver_RK4_low_storage(void)
{
	// Initialize DB Parameters
	unsigned int OutputInterval = DB.OutputInterval,
	             Neq            = DB.Neq;

	double       FinalTime = DB.FinalTime;

	// Standard datatypes
	static double rk4a[5] = { 0.0,              -0.417890474499852, -1.192151694642677, -1.697784692471528, -1.514183444257156 },
	              rk4b[5] = { 0.149659021999229, 0.379210312999627,  0.822955029386982,  0.699450455949122,  0.153057247968152 },
	              rk4c[5] = { 0.0,               0.149659021999229,  0.370400957364205,  0.622255763134443,  0.958282130674690 };

	char         *dummyPtr_c[2];
	unsigned int i, iMax, tstep, rk,
	             NvnS;
	double       time, dt, maxRHS0, maxRHS, *RES, *RHS, *What;

	struct S_VOLUME *VOLUME;

	for (i = 0; i < 2; i++)
		dummyPtr_c[i] = malloc(STRLEN_MIN * sizeof *dummyPtr_c[i]); // free

// Need to improve how dt is selected! Likely based on characteristic speeds (see nodalDG code for one possibility).
	dt = pow(0.5,DB.ML+DB.PGlobal)*0.5;
//	dt = pow(0.5,10.0);

	tstep = 0; time = 0.0;
	while (time < FinalTime) {
		update_VOLUME_Ops();


		if (time+dt > FinalTime)
			dt = FinalTime-time;

		for (rk = 0; rk < 5; rk++) {
			// Build the RHS (== -Residual)
			explicit_VOLUME_info();
//exit(1);
			explicit_FACET_info();
			maxRHS = finalize_RHS();

			// Update What
			for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
				NvnS = VOLUME->NvnS;

				RES  = VOLUME->RES;
				RHS  = VOLUME->RHS;
				What = VOLUME->What;

				for (iMax = Neq*NvnS; iMax--; ) {
					*RES    *= rk4a[rk];
					*RES    += dt*(*RHS++);
					*What++ += rk4b[rk]*(*RES++);
				}
				free(VOLUME->RHS);
//printf("%d %d\n",VOLUME->indexg,NvnS);
//array_print_d(NvnS,Neq,VOLUME->What,'C');
			}
//			exit(1);
		}
//exit(1);
		time += dt;

		// Output to paraview
		if (tstep % OutputInterval == 0 || tstep < 3) {
			sprintf(dummyPtr_c[1],"%d",tstep);
			strcpy(dummyPtr_c[0],"SolStart");
			strcat(dummyPtr_c[0],dummyPtr_c[1]);
			output_to_paraview(dummyPtr_c[0]);
		}

		// Display solver progress
		if (!tstep)
			maxRHS0 = maxRHS;

		printf("Complete: % 7.2f%%, tstep: %8d, maxRHS: % .3e\n",100*time/FinalTime,tstep,maxRHS);

		// Additional exit conditions
		if (maxRHS0/maxRHS > 1e6 && tstep) {
			printf("Exiting: maxRHS dropped by 6 orders.\n");
			break;
		}

		tstep++;
	}

	for (i = 0; i < 2; i++)
		free(dummyPtr_c[i]);
}
