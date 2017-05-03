// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_integration_conv_order.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"

#include "adaptation.h"
#include "output_to_paraview.h"
#include "initialize_test_case.h"
#include "solver_Poisson.h"
#include "solver_explicit.h"
#include "solver_implicit.h"
#include "compute_errors.h"
#include "element_functions.h"
#include "finalize_LHS.h"
#include "array_norm.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for (conv)ergence (order) integration testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static char *get_fNameOut (char *output_type);
static void h_adapt_test  (void);

static void set_test_convorder_data(struct S_convorder *const data, char const *const TestName)
{
	/*
	 *	Comments:
	 *		As nodes having integration strength of order 2*P form a basis for the polynomial space of order P for TP
	 *		elements and SI elements (P <= 2), non-trivial L2 projection error requires:
	 *			TP         : InOrder_add > 1
	 *			SI (P <= 2): InOrder_add > 0
	 */

	// default values
	data->PrintEnabled   = 0;
	data->Compute_L2proj = 0;
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

	if (strstr(TestName,"Poisson")) {
		data->AdaptiveRefine = 0;
		data->MLMax = 4;
		data->PG_add        = 0;
		data->IntOrder_add  = 2; // See comments
		if (strstr(TestName,"n-Ellipsoid_HollowSection")) {
			if (strstr(TestName,"TRI")) {
				strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_CurvedTRI");
			} else if (strstr(TestName,"QUAD")) {
				strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ellipsoid_HollowSection_CurvedQUAD");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			// ToBeModified (Include all relevant tests)
//	strcpy(argvNew[1],"test/Test_Poisson_dm1-Spherical_Section_2D_mixed");
//	strcpy(argvNew[1],"test/Test_Poisson_dm1-Spherical_Section_2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_Ellipsoidal_Section_2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_Ellipsoidal_Section_2D_QUAD");
//	strcpy(argvNew[1],"test/Test_Poisson_Ringleb2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_Ringleb2D_QUAD");
//	strcpy(argvNew[1],"test/Test_Poisson_HoldenRamp2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_GaussianBump2D_TRI");
//	strcpy(argvNew[1],"test/Test_Poisson_GaussianBump2D_mixed");
//	strcpy(argvNew[1],"test/Test_Poisson_3D_TET");
//	strcpy(argvNew[1],"test/Test_Poisson_3D_HEX");
//	strcpy(argvNew[1],"test/Test_Poisson_mixed3D_TP");
//	strcpy(argvNew[1],"test/Test_Poisson_mixed3D_HW");
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"Euler")) {
//		data->PrintEnabled   = 1;
		if (strstr(TestName,"n-Cylinder_HollowSection")) {
			data->SolveExplicit = 0;
			if (strstr(TestName,"ToBeCurved")) {
				if (strstr(TestName,"MIXED2D")) {
//					data->PrintEnabled = 0;
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
//				data->PrintEnabled = 0;
				printf("Modified Parameters.\n"); PRINT_FILELINE;
				data->PMax = 2;
				data->MLMax = 2;
				strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-Cube")) {
			data->SolveImplicit = 0;
			if (strstr(TestName,"Curved")) {
				EXIT_UNSUPPORTED;
			} else {
				if (strstr(TestName,"QUAD")) {
					data->PMin = 1;
					data->PMax = 3;
					data->MLMax = 5;
					strcpy(data->argvNew[1],"test/Euler/Test_Euler_PeriodicVortex_QUAD");
				} else {
					EXIT_UNSUPPORTED;
				}
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"NavierStokes")) {
		data->PG_add = 0;
		data->PrintEnabled = 1;
data->MLMax = 0;
data->PMax  = 1;

		strcpy(data->argvNew[1],"test/NavierStokes/Test_NavierStokes_");
		if (strstr(TestName,"n-Cylinder_Hollow")) {
			strcat(data->argvNew[1],"TaylorCouette_");
			if (strstr(TestName,"ToBeCurved")) {
				strcat(data->argvNew[1],"ToBeCurved");
				if (strstr(TestName,"TRI")) {
					strcat(data->argvNew[1],"TRI");
				} else if (strstr(TestName,"QUAD")) {
					strcat(data->argvNew[1],"QUAD");
				} else {
					EXIT_UNSUPPORTED;
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(TestName,"n-Cube")) {
			strcat(data->argvNew[1],"PlaneCouette_");
			if (strstr(TestName,"Straight")) {
				if (strstr(TestName,"QUAD")) { strcat(data->argvNew[1],"QUAD"); }
				else { EXIT_UNSUPPORTED; }
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		printf("%s\n",TestName); EXIT_UNSUPPORTED;
	}
}

void test_conv_order(struct S_convorder *const data, char const *const TestName)
{
	/*
	 *	Expected Output:
	 *		Optimal convergence orders in L2 for the solution (P+1) and its gradients (P).
	 */

	set_test_convorder_data(data,TestName);

	int  const               nargc   = data->nargc;
	char const *const *const argvNew = (char const *const *const) data->argvNew;

	bool const PrintEnabled   = data->PrintEnabled,
	           Compute_L2proj = data->Compute_L2proj,
	           SolveExplicit  = data->SolveExplicit,
	           SolveImplicit  = data->SolveImplicit,
	           AdaptiveRefine = data->AdaptiveRefine;

	unsigned int const Adapt = data->Adapt,
	                   PMin  = data->PMin,
	                   PMax  = data->PMax,
	                   MLMin = data->MLMin,
	                   MLMax = data->MLMax;

	TestDB.PGlobal       = 1; // ToBeModified (Not important for Adapt != ADAPT_0)
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	double *mesh_quality = malloc((MLMax-MLMin+1) * sizeof *mesh_quality); // free

	if (!SolveExplicit && !SolveImplicit)
		EXIT_UNSUPPORTED; // Enabled at least one

	unsigned int pass = 0;
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

			if (!SolveImplicit)
				initialize_test_case(0);
		} else {
			code_startup(nargc,argvNew,0,1);
		}

		if (strstr(TestName,"Poisson")) {
			if (Compute_L2proj) { // Compute errors of L2 projection of the exact solution
				initialize_test_case(0);

				// Output to paraview
				if (TestDB.ML <= 1 || (TestDB.PGlobal == 1) || (TestDB.PGlobal == 5 && TestDB.ML <= 4)) {
					char *const fNameOut = get_fNameOut("SolFinal_");
					output_to_paraview(fNameOut);
					free(fNameOut);
				}
			} else {
				solver_Poisson(PrintEnabled);
			}
		} else if (strstr(TestName,"Euler") || strstr(TestName,"NavierStokes")) {
//			if (SolveExplicit)
//			if (!SolveImplicit || (ML == MLMin && P == PMin))
			if (SolveExplicit && (!SolveImplicit || (ML <= MLMin+1)))
				solver_explicit(PrintEnabled);

			if (SolveImplicit)
				solver_implicit(PrintEnabled);
		}

		// Output mesh edges to paraview
		if (TestDB.PGlobal == 3 && TestDB.ML <= 2) {
			char *const fNameOut = get_fNameOut("MeshEdges_");
			output_to_paraview(fNameOut);
			free(fNameOut);
		}

		compute_errors_global();

		if (PrintEnabled) {
			compute_dof();
			printf("ML, P, dof: %zu %zu %d\n",ML,P,DB.dof);
		}

		if (P == PMin)
			evaluate_mesh_regularity(&mesh_quality[ML-MLMin]);

		if (P == PMax && ML == MLMax) {
			check_convergence_orders(MLMin,MLMax,PMin,PMax,&pass,PrintEnabled);
			check_mesh_regularity(mesh_quality,MLMax-MLMin+1,&pass,PrintEnabled);
		}

		if (Adapt == ADAPT_0) {
			set_PrintName("conv_orders",data->PrintName,&data->TestTRI);
			code_cleanup();
		}
	}}
	if (Adapt != ADAPT_0) {
		set_PrintName("conv_orders",data->PrintName,&data->TestTRI);
		code_cleanup();
	}
	free(mesh_quality);

	test_print2(pass,data->PrintName);
}

char *get_fNameOut(char *output_type)
{
	char string[STRLEN_MIN], *fNameOut;

	fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut); // keep

	strcpy(fNameOut,output_type);
	sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
								   strcat(fNameOut,DB.MeshType);
	// Choose one of the two options below to clean this up (probably DB and not TestDB) (ToBeDeleted)
	if (DB.Adapt == ADAPT_0) {
		EXIT_UNSUPPORTED;
		sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
		sprintf(string,"P%d_",DB.PGlobal); strcat(fNameOut,string);
	} else {
		sprintf(string,"_ML%d",TestDB.ML); strcat(fNameOut,string);
		sprintf(string,"P%d_",TestDB.PGlobal); strcat(fNameOut,string);
	}

	return fNameOut;
}

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

	if (strstr(Geometry,"n-Cylinder") || strstr(Geometry,"n-Cube")) {
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
