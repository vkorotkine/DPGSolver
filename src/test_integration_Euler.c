// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "S_DB.h"
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
	PMin  = 1; PMax  = 4;
	MLMin = 0; MLMax = 3;
TestDB.PGlobal = 1;

	mesh_quality = malloc((MLMax-MLMin+1) * sizeof *mesh_quality); // free

	AdaptiveRefine = 1;
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

//		DB.PGlobal = P;
//		DB.ML = ML;

		if (Adapt != ADAPT_0) {
			if (ML == MLMin) {
				mesh_to_level(TestDB.ML);
				if (AdaptiveRefine) {
					struct S_VOLUME *VOLUME;
					unsigned int Ind00, Ind01, Ind10, Ind11;
					Ind00 = 100; Ind01 = 100; Ind10 = 100; Ind11 = 100;
					if (strstr(DB.MeshType,"ToBeCurvedTRI")) {
						Ind00 = 0; Ind01 = 1; Ind10 = 1;
					} else if (strstr(DB.MeshType,"CurvedTRI")) {
//						Ind00 = 5; Ind01 = 6; Ind10 = 9;// Ind11 = 6; // Refine trailing edge
//						Ind00 = 2; Ind01 = 1; Ind10 = 5;// Ind11 = 6; // Refine leading edge
						Ind00 = 2; Ind01 = 0; Ind10 = 5;// Ind11 = 6; // Refine leading edge (alternate)
					} else if (strstr(DB.MeshType,"ToBeCurvedQUAD")) {
						Ind00 = 0; Ind10 = 3; Ind11 = 3;
					} else if (strstr(DB.MeshType,"CurvedQUAD")) {
						Ind00 = 1; Ind10 = 1;
					}

					for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
						if (VOLUME->indexg == Ind00 || VOLUME->indexg == Ind01) {
							VOLUME->Vadapt = 1;
							VOLUME->adapt_type = HREFINE;
						}
					}
					mesh_update();
					for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
						if (VOLUME->indexg == Ind10 || VOLUME->indexg == Ind11) {
							VOLUME->Vadapt = 1;
							VOLUME->adapt_type = HREFINE;
						}
					}
					mesh_update();
				}
			} else {
				mesh_h_adapt(1,'r');
			}
			mesh_to_order(TestDB.PGlobal);
		} else {
			code_startup(nargc,argvNew,0,1);
		}

		// Output mesh edges to paraview
		if (TestDB.PGlobal == 4 && TestDB.ML <= 2) {
			char *fNameOut = get_fNameOut("MeshEdges_");
			output_to_paraview(fNameOut);
			free(fNameOut);
		}
//		if (ML == MLMin || P == PMin)
//		if (ML == MLMin)
		if (0)
			solver_explicit();
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
