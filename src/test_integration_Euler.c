// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "adaptation.h"
#include "output_to_paraview.h"
#include "solver_explicit.h"
#include "solver_implicit.h"
#include "compute_errors.h"
#include "array_free.h"
#include "test_integration_Poisson.h"
#include "element_functions.h"
#include "array_norm.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Test optimal convergence of the Euler solver implementation.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

static void h_adapt(void)
{
	/*
	 * Purpose:
	 *		Perform local mesh refinement around specified XYZref coordinate locations.
	 */

	// Initialize DB Parameters
	char         *Geometry = DB.Geometry;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int NrefMax = 3, MLMax = 5;

	unsigned int Nref, NML[NrefMax], CurvedOnly[NrefMax];
	double       *XYZref;

	XYZref = malloc(NrefMax*DMAX * sizeof *XYZref); // free

	if (TestDB.ML > 0)
		printf("Error: Only enter for ML == 0.\n"), EXIT_MSG;

	if (strstr(Geometry,"JoukowskiSymmetric")) {
		double a  = DB.JSa,
		       xL = DB.JSxL;

		Nref = 2;

		NML[0] = 0;
		NML[1] = 3;

		XYZref[0+0*DMAX] = xL;  XYZref[1+0*DMAX] = 0.0; XYZref[2+0*DMAX] = 0.0;
		XYZref[0+1*DMAX] = 2*a; XYZref[1+1*DMAX] = 0.0; XYZref[2+1*DMAX] = 0.0;
	} else if (strstr(Geometry,"Ellipsoidal_Section")) {
		Nref = 3;

		unsigned int i = 0;
		NML[i] = 0; CurvedOnly[i] = 0; i++;
		NML[i] = 2; CurvedOnly[i] = 1; i++;
		NML[i] = 0; CurvedOnly[i] = 0; i++;

		i = 0;
		XYZref[0+i*DMAX] = DB.aIn;  XYZref[1+i*DMAX] = 0.0;     XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] = 0.0;     XYZref[1+i*DMAX] = DB.bOut; XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] = 0.0;     XYZref[1+i*DMAX] = DB.bIn;  XYZref[2+i*DMAX] = 0.0; i++;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	// Store indices of VOLUMEs to be updated and the adaptation type (HREFINE here) in hp_update
	unsigned int ML;
	for (ML = 0; ML < MLMax; ML++) {
		unsigned int NVglobal = DB.NVglobal;

		struct S_VOLUME *VOLUME;

		unsigned int *hp_update;
		hp_update = calloc(NVglobal , sizeof *hp_update); // free

		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int n, ve, dim, Nve, indexg;
			double       *XYZ_vV;

			struct S_ELEMENT *ELEMENT;

			ELEMENT = get_ELEMENT_type(VOLUME->type);

			Nve = ELEMENT->Nve;

			indexg = VOLUME->indexg;

			// Store XYZ_vV as a row vector
			XYZ_vV = malloc(Nve*d * sizeof *XYZ_vV); // free
			for (ve = 0; ve < Nve; ve++) {
			for (dim = 0; dim < d; dim++) {
				XYZ_vV[ve*d+dim] = VOLUME->XYZ_vV[ve+dim*Nve];
			}}

			for (n = 0; n < Nref; n++) {
				if (VOLUME->level < NML[n]) {
					if (CurvedOnly[n] && !VOLUME->curved)
						continue;

					// Check if one of the XYZ_vV matches any of the specified XYZref
					for (ve = 0; ve < Nve; ve++) {
						if (array_norm_diff_d(d,&XYZref[n*DMAX],&XYZ_vV[ve*d],"Inf") < EPS) {
							hp_update[indexg] = HREFINE;
							break;
						}
					}
				}
			}
		}

		// Add additional indices to the list of hp_update to ensure that the resulting mesh will not be more than
		// 1-irregular.
		ensure_1irregular(hp_update);

		// Mark VOLUMEs for refinement and perform the mesh update
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (hp_update[VOLUME->indexg] == 0) {
				// Do nothing
			} else if (hp_update[VOLUME->indexg] == HREFINE) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = HREFINE;
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}
		}
		free(hp_update);
		mesh_update();
	}
	free(XYZref);
}

void test_integration_Euler(int nargc, char **argv)
{
	unsigned int pass = 0;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *		Meshes for a curved Poisson problem.
	 *
	 *	Expected Output:
	 *		Optimal convergence orders in L2 for the solution (P+1).
	 *
	 */

	unsigned int P, ML, PMin, PMax, MLMin, MLMax, Adapt, AdaptiveRefine;
	double       *mesh_quality;

	strcpy(argvNew[1],"test/Test_Euler_2D_TRI");


	TestDB.PG_add = 1;
	TestDB.IntOrder_add  = 0; // > 1 for non-zero error for L2 projection on TP elements
	TestDB.IntOrder_mult = 2;

	// Convergence orders
	PMin  = 1; PMax  = 5;
	MLMin = 0; MLMax = 2;
TestDB.PGlobal = PMin;

	mesh_quality = malloc((MLMax-MLMin+1) * sizeof *mesh_quality); // free

	AdaptiveRefine = 0;
//	Adapt = ADAPT_0;
	Adapt = ADAPT_HP;
	if (Adapt != ADAPT_0) {
		TestDB.ML = DB.ML;
		code_startup(nargc,argvNew,0,2);
	}

	for (P = PMin; P <= PMax; P++) {
	for (ML = MLMin; ML <= MLMax; ML++) {
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		if (Adapt != ADAPT_0) {
			if (ML == MLMin) {
				mesh_to_level(ML);
				if (AdaptiveRefine)
					h_adapt();
			} else {
				mesh_h_adapt(1,'r');
			}
			mesh_to_order(P);
		} else {
			code_startup(nargc,argvNew,0,1);
		}

		// Output mesh edges to paraview
		if (TestDB.PGlobal == 3 && TestDB.ML <= 2) {
			char *fNameOut = get_fNameOut("MeshEdges_");
			output_to_paraview(fNameOut);
			free(fNameOut);
		}
//		if (ML <= MLMin+1 || P == 1)
		if (ML <= MLMin+1 || P > 1)
			solver_explicit();
//		if (!(ML == MLMax && P == 3) && P != 1)
			solver_implicit();

		compute_errors_global();

		printf("ML, P, dof: %d %d %d\n",ML,P,DB.dof);

		if (P == PMin)
			evaluate_mesh_regularity(&mesh_quality[ML-MLMin]);

		if (P == PMax && ML == MLMax) {
			check_convergence_orders(MLMin,MLMax,PMin,PMax,&pass);
			check_mesh_regularity(mesh_quality,MLMax-MLMin+1,&pass);
		}

		if (Adapt == ADAPT_0)
			code_cleanup();
	}}
	if (Adapt != ADAPT_0)
		code_cleanup();
	free(mesh_quality);

	printf("Convergence Orders - Poisson (3D - TRI  ):       ");
	test_print(pass);

	array_free2_c(2,argvNew);
}
