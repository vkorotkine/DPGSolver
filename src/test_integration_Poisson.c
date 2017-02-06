// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_integration_Poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "test_integration_linearization.h"
#include "compute_errors.h"
#include "array_norm.h"
#include "solver_Poisson.h"
#include "finalize_LHS.h"
#include "adaptation.h"
#include "initialize_test_case.h"
#include "output_to_paraview.h"

/*
 *	Purpose:
 *		Test various aspects of the Poisson solver implementation:
 *			1) Linearization
 *			2) Optimal convergence orders
 *
 *	Comments:
 *
 *		It was very difficult to find a case where it was clear that blending of a curved boundary was leading to a loss
 *		of optimal convergence, despite the potentially unbounded mapping derivatives with h-refinement required for the
 *		optimal error estimate in Ciarlet(1972) (Theorem 5). As noted in Scott(1973) (p. 54), the mapping function and
 *		all of its derivatives are bounded IN TERMS OF THE BOUNDARY CURVATUVE AND ITS DERIVATIVES, which motivated the
 *		implementation of cases with geometry possessing high element curvatuve on coarse meshes. For these cases,
 *		optimal convergence is lost until the mesh has been refined "enough" (such that the element curvature is small)
 *		at which point it is recovered. Perhaps even more significant, the L2 error of the projection of the exact
 *		solution sometimes increased with mesh refinement on coarse meshes (run with Compute_L2proj = 1), a trend which
 *		was not observed for the computed solution. There was no modification found which could serve to fix this issue
 *		(such as the use of alternate blending functions or generalizations of the high-order blending used to fix the
 *		NIELSON blending as proposed by Lenoir(1986)). The mathematical and numerical support for this observation can
 *		be found in Zwanenburg(2017).
 *
 *		As a results of the conclusions of the study performed in Zwanenburg(2017) (and based on my current
 *		understanding), it seems that any of the optimal blending functions should give analogous results (SCOTT,
 *		SZABO_BABUSKA, LENOIR) and that the NORMAL surface parametrization is best (certainly in the 2D case). For the
 *		extension to 3D blending, the approach to be taken would be to extend the SB blending to PYR elements (using a
 *		combination of GH and SB for the LINE and TRI parts of the WEDGE element). The surface parametrization for the
 *		non-edge geometry nodes could be use the NORMAL parametrization. An investigation into the optimal EDGE
 *		parametrization is still required.
 *
 *		It may also be noted that, despite converging at the same rate, the L2error of the computed solution is
 *		significantly higher than that of the L2 projected exact solution for certain polynomial orders. Comparison with
 *		results obtained based on the DPG solver may be interesting here. (ToBeModified)
 *
 *
 *	*** IMPORTANT ***   Convergence Order Testing   *** IMPORTANT ***
 *
 *		It was found that optimal convergence was not possible to obtain using a series of uniformly refined TET meshes
 *		based on the "refine by splitting" algorithm in gmsh. However, optimal orders were recovered when a series of
 *		unstructured meshes consisting of TETs of decreasing volume was used.
 *
 *		Gmsh(2.14)'s "refine by splitting" algorithm splits each TET into 8 TETs, in a manner identical to the TET8
 *		h-refinement algorithm implemented in the code. Assuming an initially regular TET, with all edges having length
 *		= 2.0, is refined in this manner, the result will be 4 regular TETs with edge length = 1.0 and 4 other TETs with
 *		5 edges having length = 1.0 and 1 edge having length = sqrt(2.0). Taking a measure of the regularity of the mesh
 *		to be the ratio of the spheres enclosing and enclosed by the TETs, the regularity bound is violated through this
 *		refinement process if the splitting direction of the internal octohedron is taken at random, but not if taken
 *		consistently along the same axis (See Lenoir(1986) for the regularity requirement). Hence an alternative method
 *		of generating the refined mesh sequence must be employed (as compared to gmsh) to achieve the optimal
 *		convergence:
 *
 *			1) A sequence of unstructured (non-nested) meshes (gmsh) with diminishing volume gave optimal orders.
 *			2) TET8 refinement along a consistent axis (this code) gave optimal orders.
 *			3) TET6 and TET12 refinement algorithms both gave sub-optimal orders to date.
 *				- Continued investigation may be made after a mesh quality improvement algorithm is implemented (e.g.
 *				  Peiro(2013)-Defining_Quality_Measures_for_Validation_and_Generation_of_High-Order_Tetrahedral_Meshes).
 *				  ToBeModified.
 *
 *	*** IMPORTANT ***   Convergence Order Testing   *** IMPORTANT ***
 *
 *		When testing with FLUX_IP, the iterative solver seems much more likely to fail based on several tests,
 *		regardless of the value selected for tau. The tolerance on the symmetry test might also need to be increased to
 *		show a passing result for this flux.
 *
 *	Notation:
 *
 *	References:
 *		Ciarlet(1972)-Interpolation_Theory_Over_Curved_Elements,_with_Applications_to_Finite_Element_Methods
 *		Scott(1973)-Finite_Element_Techniques_for_Curved_Boundaries
 *		Lenoir(1986)-Optimal_Isoparametric_Finite_Elements_and_Error_Estimates_for_Domains_Involving_Curved_Boundaries
 *		Zwanenburg(2017)-A_Necessary_High-Order_Meshing_Constraint_when_using_Polynomial_Blended_Geometry_Elements_for_
 *		                 Curved_Boundary_Representation
 */

struct S_linearization {
	Mat A, A_cs, A_csc;
	Vec b, b_cs, b_csc, x, x_cs, x_csc;
};

static void test_linearization(int nargc, char **argvNew, const unsigned int Nref, const unsigned int update_argv,
                               const char *TestName, struct S_linearization *data)
{
	unsigned int pass;

	PetscBool Symmetric;

	Mat A, A_cs, A_csc;
	Vec b, b_cs, b_csc, x, x_cs, x_csc;

	A = data->A; A_cs = data->A_cs; A_csc = data->A_csc;
	b = data->b; b_cs = data->b_cs; b_csc = data->b_csc;
	x = data->x; x_cs = data->x_cs; x_csc = data->x_csc;

	code_startup(nargc,argvNew,Nref,update_argv);

	implicit_info_Poisson();

	finalize_LHS(&A,&b,&x,0);
	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

//	MatView(A,PETSC_VIEWER_STDOUT_SELF);
//	MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);

	MatIsSymmetric(A,1e3*EPS,&Symmetric);
//	MatIsSymmetric(A_cs,1e5*EPS,&Symmetric);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf")     < 1e2*EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf") < 1e2*EPS &&
	    Symmetric)
		pass = 1, TestDB.Npass++;
	else
		printf("%e %e %d\n",PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf"),
		                    PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf"),Symmetric);

	printf("%s",TestName);
	test_print(pass);

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();
}

static char *get_fNameOut(char *output_type)
{
	char string[STRLEN_MIN], *fNameOut;

	fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut);

	strcpy(fNameOut,output_type);
	sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
								   strcat(fNameOut,DB.MeshType);
	if (DB.Adapt == ADAPT_0) {
		sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
		sprintf(string,"P%d_",DB.PGlobal); strcat(fNameOut,string);
	} else {
		sprintf(string,"_ML%d",TestDB.ML); strcat(fNameOut,string);
		sprintf(string,"P%d_",TestDB.PGlobal); strcat(fNameOut,string);
	}

	return fNameOut;
}

void test_integration_Poisson(int nargc, char **argv)
{
	unsigned int pass = 0;
	char         **argvNew, *TestName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	TestName   = malloc(STRLEN_MAX * sizeof *TestName); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *
	 *		Meshes for a curved Poisson problem. (ToBeModified: Potentially just coarsest mesh then h-refinement)
	 *
	 *	Expected Output:
	 *
	 *		Correspondence of LHS matrices computed using complex step and exact linearization.
	 *		Optimal convergence orders in L2 for the solution (P+1) and its gradients (P).
	 *
	 */

	unsigned int P, ML, PMin, PMax, MLMin, MLMax, Adapt, Compute_L2proj, AdaptiveRefine;
	double       *mesh_quality;
	struct S_linearization *data;

	data = calloc(1 , sizeof *data); // free


	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	// **************************************************************************************************** //
	// 2D (Mixed TRI/QUAD mesh)
	TestDB.PGlobal = 2;
	TestDB.ML      = 0;

	//              0         10        20        30        40        50
	strcpy(TestName,"Linearization Poisson (2D - Mixed):              ");
	strcpy(argvNew[1],"test/Test_Poisson_dm1-Spherical_Section_2D_mixed");

//if (0)
	test_linearization(nargc,argvNew,0,1,TestName,data);


	// **************************************************************************************************** //
	// 3D (TET mesh)
	TestDB.PGlobal = 2;
	TestDB.ML      = 0;

	//              0         10        20        30        40        50
	strcpy(TestName,"Linearization Poisson (3D - TET):                ");
	strcpy(argvNew[1],"test/Test_Poisson_3D_TET");

if (0) // The 3D testing needs to be updated (ToBeDeleted)
	test_linearization(nargc,argvNew,0,1,TestName,data);


	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
//	strcpy(argvNew[1],"test/Test_Poisson_dm1-Spherical_Section_2D_mixed");
	strcpy(argvNew[1],"test/Test_Poisson_dm1-Spherical_Section_2D_TRI");
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

	TestDB.PG_add = 0;
	TestDB.IntOrder_add  = 2; // > 1 for non-zero error for L2 projection on TP elements
	TestDB.IntOrder_mult = 2;

	// Convergence orders
	PMin  = 1; PMax  = 4;
	MLMin = 0; MLMax = 3;
TestDB.PGlobal = 1;

	mesh_quality = malloc((MLMax-MLMin+1) * sizeof *mesh_quality); // free

	Compute_L2proj = 0; // Use IntOrder_add > 0 for non-trivial P1-P2 results for L2proj
	AdaptiveRefine = 0;
//	Adapt = ADAPT_0;
	Adapt = ADAPT_HP;
	if (Adapt != ADAPT_0) {
		TestDB.ML = DB.ML;
		code_startup(nargc,argvNew,0,2);
	}

	struct S_VOLUME *VOLUME;

	for (P = PMin; P <= PMax; P++) {
	for (ML = MLMin; ML <= MLMax; ML++) {
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		if (Adapt != ADAPT_0) {
			if (ML == MLMin) {
				mesh_to_level(TestDB.ML);
				if (AdaptiveRefine) {
					unsigned int Ind00, Ind01, Ind10, Ind11;
					Ind00 = 100; Ind01 = 100; Ind10 = 100; Ind11 = 100;
					if (strstr(DB.MeshType,"ToBeCurvedTRI")) {
						Ind00 = 0; Ind01 = 1; Ind10 = 1;
					} else if (strstr(DB.MeshType,"CurvedTRI")) {
						Ind00 = 3; Ind10 = 4; Ind11 = 4;
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

		if (Compute_L2proj) { // Compute errors of L2 projection of the exact solution
			char *string, *fNameOut;

			fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut); // free
			string   = malloc(STRLEN_MIN * sizeof *string);   // free

			initialize_test_case(0);
			// Output to paraview
			if (TestDB.ML <= 1 || (TestDB.PGlobal == 1) || (TestDB.PGlobal == 5 && TestDB.ML <= 4)) {
				fNameOut = get_fNameOut("SolFinal_");
				output_to_paraview(fNameOut);

				if (TestDB.PGlobal == 5 && TestDB.ML <= 2) {
					free(fNameOut); fNameOut = get_fNameOut("MeshEdges_");
					output_to_paraview(fNameOut);
				}
			}

			free(fNameOut);
			free(string);
		} else {
			solver_Poisson();
		}
		compute_errors_global();

		printf("dof: %d\n",DB.dof);

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



	free(argvNew[0]); free(argvNew[1]); free(argvNew);
	free(TestName);
	free(data);
}
