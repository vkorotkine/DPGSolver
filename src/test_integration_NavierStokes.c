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
#include "test_integration_Euler.h"

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

struct S_convorder {
	bool         PrintEnabled, SolveExplicit, SolveImplicit, AdaptiveRefine, TestTRI;
	unsigned int PMin, PMax, MLMin, MLMax, Adapt, PG_add, IntOrder_add, IntOrder_mult;
	char         **argvNew, *PrintName;
};


static void set_test_convorder_data(struct S_convorder *data, const char *TestName)
{
	// default values
	data->PrintEnabled   = 1;
	data->SolveExplicit  = 1;
	data->SolveImplicit  = 1;
	data->AdaptiveRefine = 1;
	data->Adapt = ADAPT_HP;

	data->PMin  = 1;
	data->PMax  = 3;
	data->MLMin = 0;
	data->MLMax = 3;

	data->PG_add        = 1;
	data->IntOrder_add  = 0;
	data->IntOrder_mult = 2;


	if (strstr(TestName,"n-Cylinder_Hollow")) {
		data->SolveImplicit = 0;
		if (strstr(TestName,"ToBeCurved")) {
			if (strstr(TestName,"TRI")) {
				data->PrintEnabled = 1;
				strcpy(data->argvNew[1],"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedTRI");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void test_convorder(int nargc, char **argvNew, const char *TestName, struct S_convorder *data)
{
	/*
	 *	Comments:
	 *		This function could be used for testing the Euler equation solver as well by setting DB.Viscous = 0.
	 */

	unsigned int pass = 0;

	bool         PrintEnabled, SolveExplicit, SolveImplicit, AdaptiveRefine;
	unsigned int Adapt, PMin, PMax, MLMin, MLMax;
	double       *mesh_quality;

	set_test_convorder_data(data,TestName);

	PrintEnabled   = data->PrintEnabled;
	SolveExplicit  = data->SolveExplicit;
	SolveImplicit  = data->SolveImplicit;
	AdaptiveRefine = data->AdaptiveRefine;
	Adapt          = data->Adapt;

	PMin  = data->PMin;  PMax  = data->PMax;
	MLMin = data->MLMin; MLMax = data->MLMax;

	TestDB.PGlobal       = 3;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	mesh_quality = malloc((MLMax-MLMin+1) * sizeof *mesh_quality); // free

	if (Adapt != ADAPT_0) {
		TestDB.ML = DB.ML;
		code_startup(nargc,argvNew,0,2);
	}

	for (size_t P = PMin; P <= PMax; P++) {
	for (size_t ML = MLMin; ML <= MLMax; ML++) {
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		if (Adapt != ADAPT_0) {
			if (ML == MLMin) {
				mesh_to_level(ML);
				if (AdaptiveRefine)
					h_adapt_test();
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

		if (SolveExplicit)
			solver_explicit();

		if (SolveImplicit)
			solver_implicit(PrintEnabled);

		compute_errors_global();

		if (PrintEnabled)
			printf("ML, P, dof: %zu %zu %d\n",ML,P,DB.dof);

		if (P == PMin)
			evaluate_mesh_regularity(&mesh_quality[ML-MLMin]);

		if (P == PMax && ML == MLMax) {
			check_convergence_orders(MLMin,MLMax,PMin,PMax,&pass,PrintEnabled);
			check_mesh_regularity(mesh_quality,MLMax-MLMin+1,&pass,PrintEnabled);
		}

		if (Adapt == ADAPT_0)
			code_cleanup();
	}}
	set_PrintName("conv_orders",data->PrintName,&data->TestTRI);

	if (Adapt != ADAPT_0)
		code_cleanup();
	free(mesh_quality);

	test_print2(pass,data->PrintName);
}

void test_integration_NavierStokes(int nargc, char **argv)
{
	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	PrintName  = malloc(STRLEN_MAX * sizeof *PrintName); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *		Meshes for NavierStokes cases.
	 *
	 *	Expected Output:
	 *		Optimal convergence orders in L2 for the solution (P+1).
	 *
	 */

	struct S_convorder *data_c;

	data_c = calloc(1 , sizeof *data_c); // free
	data_c->argvNew   = argvNew;
	data_c->PrintName = PrintName;

	test_convorder(nargc,argvNew,"n-Cylinder_Hollow_ToBeCurvedTRI",data_c);
	test_convorder(nargc,argvNew,"n-Cylinder_Hollow_ToBeCurvedQUAD",data_c);
	test_convorder(nargc,argvNew,"n-Cylinder_Hollow_ToBeCurvedMIXED2D",data_c);

	// Add tests for:
	// 1) Analytical Solution for simple shearing motion between parallel flat plates (Illingworth(1950), problem 3)
	// 2) Stokes flow over a sphere (analytical solution for compressible?)


	array_free2_c(2,argvNew);
	free(PrintName);
	free(data_c);
}
