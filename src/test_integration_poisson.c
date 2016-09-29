// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_poisson.h"

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
#include "solver_poisson.h"
#include "finalize_LHS.h"
// ToBeDeleted: Check for unnecessary includes above

/*
 *	Purpose:
 *		Test various aspects of the Poisson solver implementation:
 *			1) Linearization
 *			2) Optimal convergence orders
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_poisson(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	TestDB.TestCase = malloc(STRLEN_MAX * sizeof *(TestDB.TestCase)); // free

	/*
	 *	Input:
	 *
	 *		ToBeModified
	 *
	 *	Expected Output:
	 *
	 *		ToBeModified
	 *
	 */

	unsigned int P, ML, PMin, PMax, MLMin, MLMax;

	Mat A = NULL, A_cs = NULL, A_csc = NULL;
	Vec b = NULL, b_cs = NULL, b_csc = NULL,
	    x = NULL, x_cs = NULL, x_csc = NULL;

	// **************************************************************************************************** //
	// TRIs (change to Mixed) ToBeModified
	strcpy(argvNew[1],"test/Test_poisson_TRI");
	strcpy(TestDB.TestCase,"Poisson");

	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	// Linearization
// ToBeModified
	TestDB.PGlobal = 1;
	TestDB.ML      = 0;

	code_startup(nargc,argvNew,0,1);

	implicit_info_Poisson();
//	finalize_LHS(&A,&b,&x,1);
//	finalize_LHS(&A,&b,&x,2);
//	finalize_LHS(&A,&b,&x,3);
//	finalize_Mat(&A,1);

//	compute_A_cs(&A_cs,&b_cs,&x_cs,1);
//	compute_A_cs(&A_cs,&b_cs,&x_cs,2);
//	compute_A_cs(&A_cs,&b_cs,&x_cs,3);
//	finalize_Mat(&A_cs,1);

//	MatView(A,PETSC_VIEWER_STDOUT_SELF);
//	MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);
//	EXIT_MSG;

	finalize_LHS(&A,&b,&x,0);
	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

MatView(A,PETSC_VIEWER_STDOUT_SELF);
MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);
	MatView(A_csc,PETSC_VIEWER_STDOUT_SELF);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A,A_cs,"Inf")  < EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A,A_csc,"Inf") < EPS)
//	if (PetscMatAIJ_norm_diff_d(DB.dof,A,A_cs,"Inf")  < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("Linearization Poisson (2D - TRI  ):              ");
	test_print(pass);
EXIT_MSG;

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();

	// Convergence orders
	PMin = 0;  PMax = 4;
	MLMin = 0; MLMax = 4;

	for (P = PMin; P <= PMax; P++) {
	for (ML = MLMin; ML <= MLMax; ML++) {
/*
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		code_startup(nargc,argvNew,0,1);

		solver_Poisson();
		compute_errors_global();

		code_cleanup();
*/
	}}
	// test with various boundary conditions and fluxes (ToBeDeleted)


	free(argvNew[0]); free(argvNew[1]); free(argvNew);
	free(TestDB.TestCase);
}
