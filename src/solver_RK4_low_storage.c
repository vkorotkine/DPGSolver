// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
//#include <string.h>

//#include "database.h"
//#include "parameters.h"
//#include "functions.h"

//#include "petscsys.h"

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
	double FinalTime = DB.FinalTime;

	// Standard datatypes
	static double rk4a[5] = { 0.0,              -0.417890474499852, -1.192151694642677, -1.697784692471528, -1.514183444257156 },
	              rk4b[5] = { 0.149659021999229, 0.379210312999627,  0.822955029386982,  0.699450455949122,  0.153057247968152 },
	              rk4c[5] = { 0.0,               0.149659021999229,  0.370400957364205,  0.622255763134443,  0.958282130674690 };
	
	unsigned int i, iMax, tstep, rk;
	double time, dt, *RES, *RHS;

	struct S_VOLUME *VOLUME;

// Need to improve how dt is selected! Likely based on characteristic speeds (see nodalDG code for one possibility).
	dt = 1e-4;


	tstep = 0; time = 0.0;
	while (time < FinalTime) {
		if (time+dt > FinalTIme)
			dt = FinalTime-time;

		for (rk = 0; rk < 5; rk++) {
			// Build the RHS (== -Residual)
//			explicit_VOLUME_info();
//			explicit_FACET_info();
//			finalize_RHS();

			// Update What
			for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
				NvnS = VOLUME->NvnS;

				RES  = VOLUME->RES;
				RHS  = VOLUME->RHS;
				What = VOLUME->What;

				for (i = 0, iMax = Neq*NvnS; i < iMax; i++) {
					RES[i] = rk4a[rk]*RES[i] + dt*RHS[i];
					What[i] += rk4b[rk]*RES[i];
				}
			}
		}
		time += dt;

		// Output to paraview
		if (tstep % OutputInterval
		tstep++;
	}
}
