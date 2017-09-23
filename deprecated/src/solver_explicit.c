// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_explicit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "finalize_RHS.h"
#include "output_to_paraview.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"
#include "array_print.h"

#include "solver.h"

static void correct_What_Bezier(double *const WhatB, double const t, unsigned int const NvnS, double const *const WAvg,
                                char const type)
{
	unsigned int const d    = DB.d,
	                   Nvar = d+2;

	if (type == 'd') { // (d)ensity
		double *const UhatB = malloc(NvnS*Nvar * sizeof *UhatB); // free
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
		EXIT_UNSUPPORTED;
	}
}

void enforce_positivity_highorder(struct S_VOLUME *VOLUME)
{
	if (VOLUME->type == TET || VOLUME->type == WEDGE || VOLUME->type == PYR)
		return; // Bezier basis not yet implemented. (ToBeModified)

	/*
	 *	Purpose:
	 *		Enforce positivity of density and pressure by limitting the high-order solution.
	 *
	 *	Comments:
	 *		See Wang(2012) section 3.1 for details of the procedure.
	 *		The positivity is ensured to hold throughout the entire element by performing the check in the Bezier basis.
	 *		The compensation for the orthonormal basis scaling is defined such that TS[0] multiplied by an array equal
	 *		to the same constant value should return that value.
	 *		The function currently assumes that the conservative variables are being used.
	 *
	 *	References:
	 *		Wang-Shu(2012)-Robust_High_Order_Discontinuous_Galerkin_Schemes_for_Two-Dimensional_Gaseous_Detonations
	 */

	struct S_ELEMENT const *const ELEMENT = get_ELEMENT_type(VOLUME->type);

	unsigned int const d    = DB.d,
	                   Nvar = d+2,
	                   P    = VOLUME->P,
	                   NvnS = ELEMENT->NvnS[P];
	double const volume         = ELEMENT->volume,
	            *const TS       = ELEMENT->TS[P][P][0],
	            *const TS_vB    = ELEMENT->TS_vB[P][P][0],
	            *const TInvS_vB = ELEMENT->TInvS_vB[P][P][0];

	// *** Density *** ///

	// Find average density
	double        rhoAvg,
	       *const What    = VOLUME->What,
	       *const rho_hat = &What[0];

	// Note compensation for orthonormal basis scaling
	mm_d(CBCM,CBT,CBNT,1,1,NvnS,1.0/sqrt(volume),0.0,&TS[0],rho_hat,&rhoAvg);

	if (rhoAvg < EPS_PHYS)
		printf("Error: Average density approaching 0.\n"), EXIT_MSG;

	// Convert to the Bezier basis
	double *const WhatB    = malloc(NvnS*Nvar * sizeof *WhatB); // free
	mm_CTN_d(NvnS,Nvar,NvnS,TS_vB,What,WhatB);

	double const *const rho_hatB = &WhatB[0];

	// Find minimum value
	double rhoMin = 1.0/EPS;
	for (size_t n = 0; n < NvnS; n++) {
		if (rho_hatB[n] < rhoMin)
			rhoMin = rho_hatB[n];
	}

	// Correct rho (and related terms) if necessary
	if (rhoMin < EPS_PHYS) {
		double const t = min(1.0,(rhoAvg-EPS_PHYS)/(rhoAvg-rhoMin));
		correct_What_Bezier(WhatB,t,NvnS,&rhoAvg,'d');

		mm_CTN_d(NvnS,Nvar,NvnS,TInvS_vB,WhatB,What);
	}

	// *** Pressure *** ///
	double *const p_hatB = calloc(NvnS , sizeof *p_hatB); // free
	compute_pressure(WhatB,p_hatB,d,NvnS,1,'c');

	// Correct W if necessary
	double pMin = 1.0/EPS;
	for (size_t n = 0; n < NvnS; n++) {
		if (p_hatB[n] < pMin)
			pMin = p_hatB[n];
	}
	free(p_hatB);

	if (pMin < 0.0) {
		// Compute average pressure
		double *const WAvg = malloc(Nvar * sizeof *WAvg); // free

		// Note compensation for orthonormal basis scaling
		mm_d(CBCM,CBT,CBNT,1,Nvar,NvnS,1.0/sqrt(volume),0.0,&TS[0],What,WAvg);

		double pAvg;
		compute_pressure(WAvg,&pAvg,d,1,1,'c');

		if (pAvg < 0.0) {
			printf("\n\n\np: %e\nWAvg:\n",pAvg);
			array_print_d(1,Nvar,WAvg,'R');
			printf("Error: Negative average pressure.\n"), EXIT_MSG;
		}

		double const t = pAvg/(pAvg-pMin);
		correct_What_Bezier(WhatB,t,NvnS,WAvg,'p');
		free(WAvg);

		// Correct What
		mm_CTN_d(NvnS,Nvar,NvnS,TInvS_vB,WhatB,What);
	}
	free(WhatB);
}

struct S_timestepping {
	double dt, exit_tol, exit_ratio;
};

static void select_timestepping_parameters(struct S_timestepping *data)
{
	// Need to improve how dt is selected! Likely based on characteristic speeds (see nodalDG code for one possibility).
	// (ToBeDeleted)
	// The parameters set below were used for the Euler ellipsoidal section case. ToBeDeleted.
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
*/
	char const *const PDE      = DB.PDE,
	           *const TestCase = DB.TestCase;

	if (strstr(PDE,"NavierStokes")) {
		if (strstr(TestCase,"TaylorCouette")) {
			if (TestDB.Active) {
				unsigned int const ML = TestDB.ML,
				                   P  = TestDB.PGlobal;
				if (P == 1) { // Nodal basis (16 TRIs on ML = 0)
					if      (ML <= 0) { data->dt = 1e-1; }
					else if (ML <= 1) { data->dt = 1e-2; }
					else              { EXIT_UNSUPPORTED; }
				} else if (P <= 2) {
					if      (ML <= 0) { data->dt = 5e-2; }
					else if (ML <= 1) { data->dt = 1e-2; }
					else              { EXIT_UNSUPPORTED; }
				} else if (P == 3) {
					if      (ML <= 0) { data->dt = 2e-2; }
					else if (ML <= 1) { data->dt = 5e-3; }
					else              { EXIT_UNSUPPORTED; }
				} else {
					EXIT_UNSUPPORTED;
				}
				if (DB.mu <= 1e-3) {
					data->dt *= pow(0.5,1.0);
					data->exit_tol = 4e-5;
				} else if (DB.mu == 1e-0) {
					data->dt *= pow(0.5,9.0);
					data->exit_tol = 1e-3;
				} else {
					EXIT_UNSUPPORTED;
				}
			} else {
				EXIT_UNSUPPORTED;
			}
			data->exit_ratio = 1.0/EPS;
		} else if (strstr(TestCase,"PlaneCouette")) {
			if (TestDB.Active) {
				unsigned int const ML = TestDB.ML,
				                   P  = TestDB.PGlobal;
				if (P == 1) { // Nodal basis (16 TRIs on ML = 0)
					if      (ML <= 1) { data->dt = 1e-0; }
					else              { EXIT_UNSUPPORTED; }
				} else {
					EXIT_UNSUPPORTED;
				}
				data->dt *= pow(0.5,4.0);
			} else {
				EXIT_UNSUPPORTED;
			}
			data->exit_tol   = 1e-6;
			data->exit_ratio = 1.0/EPS;
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(PDE,"Euler")) {
		if (strstr(TestCase,"PeriodicVortex")) {
//			printf("Using default value for timestepping parameters.\n");
			data->dt         = 1e-3; // Selected for stability of P3 ML5 (where ML0 has 4 QUADs)
			data->exit_tol   = EPS;
			data->exit_ratio = 1.0/EPS;
		} else if (strstr(TestCase,"SupersonicVortex")) {
			printf("Using default value for timestepping parameters.\n");
			data->dt         = 1e-0;
			data->exit_tol   = 1e-6;
			data->exit_ratio = 1.0/EPS;
		} else if (strstr(TestCase,"ParabolicPipe")) {
			printf("Using default value for timestepping parameters.\n");
			data->dt         = 1e-3;
			data->exit_tol   = 1e-3;
			data->exit_ratio = 1.0/EPS;
		} else if (strstr(TestCase,"EllipticPipe")) {
			printf("Using default value for timestepping parameters.\n");
			data->dt         = 1e-3;
			data->exit_tol   = 1e-3;
		} else if (strstr(TestCase,"GaussianBump")) {
			data->dt         = 1e-3;
			data->exit_tol   = 1e-3;
			data->exit_ratio = 1.0/EPS;
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void update_RHS(double *maxRHS, bool const PrintEnabled)
{
	/*
	 *	Purpose:
	 *		Update the RHS terms for each VOLUME.
	 *
	 *	Comments:
	 *		RHS as used here is defined as all terms of the discretization which are not associated with the unsteady
	 *		term (i.e. d/dt What == RHS. With this notation, RHS thus corresponds to the negative of the steady
	 *		residual.
	 */

	struct S_solver_info solver_info = constructor_solver_info(PrintEnabled,false,false,'E',DB.Method);

	compute_RLHS(&solver_info);

	if (PrintEnabled) { printf("F "); } *maxRHS = finalize_RHS();

}

static void perform_timestepping(double const dt, double *maxRHS, bool const PrintEnabled)
{
	/*
	 *	Purpose:
	 *		Perform the timestepping for the current time step using the selected scheme.
	 *
	 *	Comments:
	 *		An positivity enforcement is included here based on the method of Wang(2012). Based on the discussion of
	 *		section 3.2, the time step must be dynamically adapted such that the average density and pressure values are
	 *		not negative at any RK stage. They recommend achieving this by reducing the timestep by two whenever such an
	 *		occurence is encountered.
	 *		For improved efficiency, it would be advantageous to reorder the VOLUMEs in order of increasing stable
	 *		timestep limit such that the potentially problematic VOLUMEs are treated first and unnecessary computations
	 *		are not performed. Ideally, the limiting VOLUMEs would be driven all the way to the end of the timestep
	 *		using the minimal possible stencil. This would,however, require the storage of the intermediate states for
	 *		many VOLUMEs.
	 *		Currently, non of this functionality is implemented but it would likely be very advantageous for difficult
	 *		unsteady cases. (ToBeModified)
	 *
	 *		The currently supported options are:
	 *			1) (R)unge-(K)utta (3)rd order (S)trong (S)tability (P)reserving
	 *			2) (R)unge-(K)utta (4)th order (L)ow (S)torage
	 *				rk4c is only needed if there is a time-dependent term in the residual (e.g. a time-dependent source
	 *				term).
	 *			3) Forward (EULER) - 1st order
	 *
	 *	References:
	 *		Gottlieb(2001)-Strong Stability-Preserving High-Order Time Discretization Methods (eq. (4.2))
	 *		Carpenter(1994)-Fourth-Order 2N-Storage Runge-Kutta Schemes
	 *		Wang-Shu(2012)-Robust High Order Discontinuous Galerkin Schemes for Two-Dimensional Gaseous Detonations
	 */

	unsigned int const d    = DB.d,
	                   Nvar = d+2;

	switch (DB.ExplicitSolverType) {
	case RK3_SSP: {
		// RES is used to store the initial solution at the beginning of the time step.
		for (size_t rk = 0; rk < 3; rk++) {
			update_RHS(maxRHS,PrintEnabled);

			// Update What
			for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				unsigned int const NvnS = VOLUME->NvnS;

				double *RES       = VOLUME->RES,
					   *What      = VOLUME->What;
				double const *RHS = VOLUME->RHS;

				if (rk == 0) {
					for (size_t iMax = Nvar*NvnS; iMax--; ) {
						*RES++   = *What;
						*What++ += dt*(*RHS++);
					}
				} else if (rk == 1) {
					for (size_t iMax = Nvar*NvnS; iMax--; ) {
						*What = 0.25*(3.0*(*RES++) + *What + dt*(*RHS++));
						What++;
					}
				} else if (rk == 2) {
					for (size_t iMax = Nvar*NvnS; iMax--; ) {
						*What = (1.0/3.0)*(*RES++ + 2.0*(*What) + 2.0*dt*(*RHS++));
						What++;
					}
				}
				enforce_positivity_highorder(VOLUME);
			}
		}
		break;
	} case RK4_LS: {
		static double const rk4a[5] = { 0.0,              -0.417890474499852, -1.192151694642677, -1.697784692471528,
		                               -1.514183444257156 },
		                    rk4b[5] = { 0.149659021999229, 0.379210312999627,  0.822955029386982,  0.699450455949122,
		                                0.153057247968152 },
		                    rk4c[5] = { 0.0,               0.149659021999229,  0.370400957364205,  0.622255763134443,
		                                0.958282130674690 };

		// silence
		if (0)
			printf("%f\n",rk4c[0]);

		for (size_t rk = 0; rk < 5; rk++) {
			update_RHS(maxRHS,PrintEnabled);

			// Update What
			for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				unsigned int const NvnS = VOLUME->NvnS;

				double *RES       = VOLUME->RES,
					   *What      = VOLUME->What;
				double const *RHS = VOLUME->RHS;

				for (size_t iMax = Nvar*NvnS; iMax--; ) {
					*RES    *= rk4a[rk];
					*RES    += dt*(*RHS++);
					*What++ += rk4b[rk]*(*RES++);
				}
				enforce_positivity_highorder(VOLUME);
			}
		}
		break;
	} case EULER: {
		update_RHS(maxRHS,PrintEnabled);

		// Update What
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int const NvnS = VOLUME->NvnS;

			double       *What = VOLUME->What;
			double const *RHS  = VOLUME->RHS;

			for (size_t iMax = Nvar*NvnS; iMax--; )
				*What++ += dt*(*RHS++);
			enforce_positivity_highorder(VOLUME);
		}
		break;
	} default: {
		EXIT_UNSUPPORTED;
		break;
	}}
}

void solver_explicit(bool const PrintEnabled)
{
	struct S_timestepping *const data_time = malloc(sizeof *data_time); // free
	select_timestepping_parameters(data_time);
	double const exit_tol   = data_time->exit_tol,
	             exit_ratio = data_time->exit_ratio;
	double       dt         = data_time->dt,
	             time       = 0.0;

	// Compute Mass matrix for uncollocated schemes
	update_VOLUME_Ops();
	update_VOLUME_finalize();

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
		enforce_positivity_highorder(VOLUME);

	output_to_paraview("ZTest_Sol_Init");

	char *dummyPtr_c[2];
	for (size_t i = 0; i < 2; i++)
		dummyPtr_c[i] = malloc(STRLEN_MIN * sizeof *dummyPtr_c[i]); // free

	unsigned int const OutputInterval = DB.OutputInterval,
	                   Adapt          = DB.Adapt,
	                   PrintTesting   = 1;
	double const FinalTime = DB.FinalTime;

	unsigned int tstep   = 0;
	double       maxRHS0 = 0.0;
	while (time < FinalTime) {
		if (Adapt && tstep)
			mesh_update();

		if (time+dt > FinalTime)
			dt = FinalTime-time;

		double maxRHS;
		perform_timestepping(dt,&maxRHS,PrintEnabled);
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

		if (PrintEnabled)
			printf("Complete: % 7.2f%%, tstep: %8d, maxRHS (no MInv): % .3e\n",100*time/FinalTime,tstep,maxRHS);

		// Additional exit conditions
		if (tstep > 2) {
			if (maxRHS0/maxRHS > exit_ratio) {
				if (PrintEnabled)
					printf("Exiting: maxRHS dropped by % .2e orders.\n",log10(exit_ratio));
				break;
			} else if (maxRHS < exit_tol) {
				if (PrintEnabled)
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
	for (size_t i = 0; i < 2; i++)
		free(dummyPtr_c[i]);
	free(data_time);

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
}
