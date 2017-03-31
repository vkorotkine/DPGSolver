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
 *		Test various aspects of the Euler solver implementation:
 *			1) Equivalence between real and complex versions of functions
 *			2) Optimal convergence orders
 *
 *	Comments:
 *		Complex versions of functions are used for complex step verification of the linearization.
 *
 *	Notation:
 *
 *	References:
 */

struct S_equivalence_rc {
	char **argvNew, *PrintName;
};

struct S_convorder {
	bool         PrintEnabled, SolveExplicit, AdaptiveRefine, TestTRI;
	unsigned int PMin, PMax, MLMin, MLMax, Adapt, PG_add, IntOrder_add, IntOrder_mult;
	char         **argvNew, *PrintName;
};

void h_adapt_test(void)
{
	/*
	 * Purpose:
	 *		Perform local mesh refinement around specified XYZref coordinate locations.
	 */

	// Initialize DB Parameters
	char         *Geometry = DB.Geometry;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int NrefMax = 5, MLMax = 5;

	unsigned int Nref, NML[NrefMax], CurvedOnly[NrefMax];
	double       *XYZref;

	XYZref = malloc(NrefMax*DMAX * sizeof *XYZref); // free

	if (TestDB.ML > 0)
		printf("Error: Only enter for ML == 0.\n"), EXIT_MSG;

	if (strstr(Geometry,"n-Cylinder")) {
		Nref = 0;
	} else if (strstr(Geometry,"JoukowskiSymmetric")) {
		double a  = DB.JSa,
		       xL = DB.JSxL;

		Nref = 2;

		NML[0] = 0;
		NML[1] = 3;

		XYZref[0+0*DMAX] = xL;  XYZref[1+0*DMAX] = 0.0; XYZref[2+0*DMAX] = 0.0;
		XYZref[0+1*DMAX] = 2*a; XYZref[1+1*DMAX] = 0.0; XYZref[2+1*DMAX] = 0.0;
	} else if (strstr(Geometry,"n-Ellipsoid")) {
		Nref = 5;

		unsigned int i = 0;
		NML[i] = 2; CurvedOnly[i] = 0; i++;
		NML[i] = 1; CurvedOnly[i] = 0; i++;
		NML[i] = 2; CurvedOnly[i] = 0; i++;
		NML[i] = 1; CurvedOnly[i] = 0; i++;
		NML[i] = 1; CurvedOnly[i] = 0; i++;

		i = 0;
		XYZref[0+i*DMAX] =  DB.aIn;  XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] =  DB.aOut; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] = -DB.aIn;  XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] = -DB.aOut; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] =  0.0;     XYZref[1+i*DMAX] = DB.bIn; XYZref[2+i*DMAX] = 0.0; i++;
	} else if (strstr(Geometry,"EllipsoidalBump")) {
		Nref = 2;

		unsigned int i = 0;
		NML[i] = 2; CurvedOnly[i] = 1; i++;
		NML[i] = 2; CurvedOnly[i] = 1; i++;

		i = 0;
		XYZref[0+i*DMAX] =  DB.aIn; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
		XYZref[0+i*DMAX] = -DB.aIn; XYZref[1+i*DMAX] = 0.0;    XYZref[2+i*DMAX] = 0.0; i++;
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
			free(XYZ_vV);
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

static void set_test_equivalence_rc_data(struct S_equivalence_rc *data, const char *TestName)
{
	// default values
	data->P     = 3;
	data->ML    = 2;
	data->Adapt = ADAPT_HP;

	data->PG_add        = 1;
	data->IntOrder_add  = 0;
	data->IntOrder_mult = 2;
}

static void test_equivalence_rc(int nargc, char **argvNew, const char *TestName, struct S_equivalence_rc *data)
{
	unsigned int Nvar = DB.Nvar;

	unsigned int Adapt;

	set_test_equivalence_rc_data(data,TestName);

	Adapt = data->Adapt;

	TestDB.PGlobal = data->P;
	TestDB.ML      = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	code_startup(nargc,argvNew,0,2);
	mesh_to_level(ML);
	mesh_to_order(P);

	// Copy What to What_c
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int NvnS = VOLUME->NvnS;

		if (VOLUME->What_c)
			free(VOLUME->What_c);
		VOLUME->What_c = malloc(NvnS*Nvar * sizeof *(VOLUME->What_c));

		for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
			VOLUME->What_c[i] = VOLUME->What[i];
		}
	}

	// Compute RHS terms using the real and complex functions
	explicit_VOLUME_info();
	explicit_VOLUME_info_c();

	// Check for equivalence
	unsigned int pass = 0;
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int NvnS = VOLUME->NvnS;

		array_print_d(NvnS,Nvar,VOLUME->RHS,'C');
		array_print_cmplx(NvnS,Nvar,VOLUME->RHS_c,'C');
		EXIT_UNSUPPORTED;
	}

//	set_PrintName_ConvOrders(data->PrintName,&data->TestTRI);
	code_cleanup();

	test_print2(pass,data->PrintName);

	// Ensure that the RHS terms are not equal to 0.0 as well
}

static void set_test_convorder_data(struct S_convorder *data, const char *TestName)
{
	// default values
	data->PrintEnabled   = 1;
	data->SolveExplicit  = 1;
	data->AdaptiveRefine = 1;
	data->Adapt = ADAPT_HP;

	data->PMin  = 1;
	data->PMax  = 3;
	data->MLMin = 0;
	data->MLMax = 3;

	data->PG_add        = 1;
	data->IntOrder_add  = 0;
	data->IntOrder_mult = 2;


	if (strstr(TestName,"n-Cylinder_HollowSection")) {
		data->SolveExplicit = 0;
		if (strstr(TestName,"ToBeCurved")) {
			if (strstr(TestName,"MIXED2D")) {
				data->PrintEnabled = 0;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED2D");
			} else if (strstr(TestName,"TET")) {
				// Starting with a coarser initial mesh than ML=2 lead to blow-up in solver_implicit. This is
				// potentially a result of poor element quality using the h-refinement of the initial mesh. Can try
				// using a refined mesh sequence from gmsh or with SolverExplicit = 1 for initial convergence.
				data->MLMax = 1;
				data->PMax  = 2;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedTET");
			} else if (strstr(TestName,"HEX")) {
				data->MLMax = 2;
				data->PMax  = 2;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedHEX");
			} else if (strstr(TestName,"WEDGE")) {
				data->MLMax = 2;
				data->PMax  = 2;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedWEDGE");
			} else if (strstr(TestName,"MIXED_TP")) {
				data->MLMax = 2;
				data->PMax  = 2;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED3D_TP");
			} else if (strstr(TestName,"MIXED_HW")) {
				data->MLMax = 2;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED3D_HW");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"CurvedMIXED2D")) {
//			data->PrintEnabled = 0;
			printf("Modified Parameters.\n");
			data->PMax = 2;
			data->MLMax = 2;
			strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void test_convorder(int nargc, char **argvNew, const char *TestName, struct S_convorder *data)
{
	unsigned int pass = 0;

	bool         PrintEnabled, SolveExplicit, AdaptiveRefine;
	unsigned int Adapt, PMin, PMax, MLMin, MLMax;
	double       *mesh_quality;

	set_test_convorder_data(data,TestName);

	PrintEnabled   = data->PrintEnabled;
	SolveExplicit  = data->SolveExplicit;
	AdaptiveRefine = data->AdaptiveRefine;
	Adapt          = data->Adapt;

	PMin  = data->PMin;  PMax  = data->PMax;
	MLMin = data->MLMin; MLMax = data->MLMax;

	TestDB.PGlobal       = 1;
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
	set_PrintName_ConvOrders(data->PrintName,&data->TestTRI);

	if (Adapt != ADAPT_0)
		code_cleanup();
	free(mesh_quality);

	test_print2(pass,data->PrintName);
}

void test_integration_Euler(int nargc, char **argv)
{
	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	PrintName  = malloc(STRLEN_MAX * sizeof *PrintName); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *		Meshes for Euler cases.
	 *
	 *	Expected Output:
	 *
	 *		Correspondence between RHS terms computed using real and complex functions.
	 *		Optimal convergence orders in L2 for the solution (P+1).
	 *
	 */

	// **************************************************************************************************** //
	// Real/Complex Equivalence
	// **************************************************************************************************** //
	struct S_equivalence_rc *data_rc;

	data_rc = calloc(1 , sizeof *data_rc); // free
	data_rc->argvNew   = argvNew;
	data_rc->PrintName = PrintName;

	test_equivalence_rc(nargc,argvNew,"n-Cylinder_HollowSection_CurvedMIXED2D",data_rc);

	free(data_rc);

	// **************************************************************************************************** //
	// Convergence Order
	// **************************************************************************************************** //
	struct S_convorder *data_c;

	data_c = calloc(1 , sizeof *data_c); // free
	data_c->argvNew   = argvNew;
	data_c->PrintName = PrintName;

	test_convorder(nargc,argvNew,"n-Cylinder_HollowSection_CurvedMIXED2D",data_c);
//	test_convorder(nargc,argvNew,"n-Cylinder_HollowSection_ToBeCurvedMIXED2D",data_c);

bool test_3D = 0;
if (test_3D) {
	test_convorder(nargc,argvNew,"n-Cylinder_HollowSection_ToBeCurvedTET",data_c);
	test_convorder(nargc,argvNew,"n-Cylinder_HollowSection_ToBeCurvedHEX",data_c); // Optimal
//	test_convorder(nargc,argvNew,"n-Cylinder_HollowSection_ToBeCurvedWEDGE",data_c); // Need to implement operators
} else {
	printf("\nWarning: 3D SupersonicVortex testing is currently disabled.\n\n"); TestDB.Nwarnings++;
}

	printf("\n\n***Add integration tests for PeriodicVortex case (Stationary and moving).***\n\n"); TestDB.Nwarnings++;



	array_free2_c(2,argvNew);
	free(PrintName);
	free(data_c);
}
