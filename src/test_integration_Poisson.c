// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_Poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACET.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "test_integration_linearization.h"
#include "compute_errors.h"
#include "array_norm.h"
#include "solver_Poisson.h"
#include "finalize_LHS.h"
#include "adaptation.h"

/*
 *	Purpose:
 *		Test various aspects of the Poisson solver implementation:
 *			1) Linearization
 *			2) Optimal convergence orders
 *
 *	Comments:
 *
 *	*** IMPORTANT ***   Convergence Order Testing   *** IMPORTANT ***
 *
 *		It was found that optimal convergence was not possible to obtain using a series of uniformly refined TET meshes
 *		based on the "refine by splitting" algorithm in gmsh. However, optimal orders were recovered when a series of
 *		unstructured meshes consisting of TETs of decreasing volume was used.
 *
 *		Gmsh(2.14)'s "refine by splitting" algorithm splits each TET into 8 TETs, in a manner identical to the original
 *		h refinement algorithm for TETs implemented in the code. Assuming an initially regular TET, with all edges having
 *		length = 2.0, is refined in this manner, the result will be 4 regular TETs with edge length = 1.0 and 4 other
 *		TETs with 5 edges having length = 1.0 and 1 edge having length = sqrt(2.0). Taking a measure of the regularity
 *		of the mesh to be the ratio of the spheres enclosing and enclosed by the TETs, the regularity bound is violated
 *		through this refinement process. (See Lenoir(1986) for the regularity requirement). Hence an alternative method
 *		of generating the refined mesh sequence must be employed to achieve the optimal convergence:
 *
 *			1) A sequence of unstructured (non-nested) meshes with diminishing volume gave optimal orders.
 *			2) Refinement based on splitting into 12 TETs or 4 TETs and 2 PYRs must be investigated. (ToBeModified)
 *
 *		The second condition above is necessary as the h-refinement algorithm will not provide optimal convergence if
 *		both of these refinement alternatives result in similar behaviour to that discussed above.
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

	MatIsSymmetric(A,1e2*EPS,&Symmetric);
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

void test_integration_Poisson(int nargc, char **argv)
{
	unsigned int pass = 0;
	char         **argvNew, *TestName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	TestName   = malloc(STRLEN_MAX * sizeof *TestName); // free

	strcpy(argvNew[0],argv[0]);

	TestDB.TestCase = malloc(STRLEN_MAX * sizeof *(TestDB.TestCase)); // free

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

	unsigned int P, ML, PMin, PMax, MLMin, MLMax, Adapt;
	double       *mesh_quality;
	struct S_linearization *data;

	data = calloc(1 , sizeof *data); // free


	strcpy(TestDB.TestCase,"Poisson");

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	// **************************************************************************************************** //
	// 2D (Mixed TRI/QUAD mesh)
	TestDB.PGlobal = 3;
	TestDB.ML      = 0;

	//              0         10        20        30        40        50
	strcpy(TestName,"Linearization Poisson (2D - Mixed):              ");
	strcpy(argvNew[1],"test/Test_Poisson_mixed2D");

TestDB.PGlobal = 1;
	test_linearization(nargc,argvNew,0,1,TestName,data);

	// **************************************************************************************************** //
	// 3D (TET mesh)
	TestDB.PGlobal = 2;
	TestDB.ML      = 0;

	//              0         10        20        30        40        50
	strcpy(TestName,"Linearization Poisson (3D - TET):                ");
	strcpy(argvNew[1],"test/Test_Poisson_3D_TET");

//	test_linearization(nargc,argvNew,0,1,TestName,data);


	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
//	strcpy(argvNew[1],"test/Test_Poisson_mixed2D");
	strcpy(argvNew[1],"test/Test_Poisson_3D_TET");
//	strcpy(argvNew[1],"test/Test_Poisson_mixed3D_TP");
//	strcpy(argvNew[1],"test/Test_Poisson_mixed3D_HW");

/*
 *	To Do:
 *		Implement the curved geometry treatment based on blending.
 *			Working in 2D. => Optimal orders (identical results to gmsh's refine by splitting).
 *			Working in 3D using straight elements and isotropic TET refinement (TET -> 8 TETs) => Optimal orders (P1) or
 *			when project_to_sphere not used (All P). Note that the mesh regularity is bounded in this case for
 *			h-refinement implemented here but not for that of the "refine by splitting" meshes generated by gmsh.
 *			PYR refinement should also be investigated while the HEX and WEDGE refinements would not pose any problems
 *			with regards to this issue.
 *		Check literature (Mavriplis (FV), Galbraith) for mention of this issue in TET verification studies.
 *		Investigate cube_to_sphere effects on mesh regularity; problems are only being observed in 3D.
 */
	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	// Convergence orders
	PMin = 1;  PMax = 2;
	MLMin = 0; MLMax = 3;

	mesh_quality = malloc((MLMax-MLMin+1) * sizeof *mesh_quality); // free

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
			mesh_to_level(TestDB.ML);
			mesh_to_order(TestDB.PGlobal);
		} else {
			code_startup(nargc,argvNew,0,1);
		}

		solver_Poisson();
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
	// test with all fluxes as they are implemented (ToBeModified)

	printf("Convergence Orders - Poisson (3D - TET  ):       ");
	test_print(pass);



	free(argvNew[0]); free(argvNew[1]); free(argvNew);
	free(TestName);
	free(data);
	free(TestDB.TestCase);
}
