// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_integration_L2_projection_errors.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Test.h"
#include "S_DB.h"

#include "test_code_integration.h"
#include "compute_errors.h"

/*
 *	Purpose:
 *		Test L2 error convergence orders using L2 projections of the exact solution to polynomial spaces of order P.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_L2_projection_errors(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	DB.TestL2projection = 1;

	/*
	 *	Input:
	 *
	 *		Mesh files on which to compute errors.
	 *
	 *	Expected Output:
	 *
	 *		Optimal convergence orders using a sequence of uniformly refined meshes. (ToBeModified)
	 *
	 */

	unsigned int P, ML, PMin, PMax, MLMin, MLMax;

	// **************************************************************************************************** //
	// LINEs


	// **************************************************************************************************** //
	// TRIs
	strcpy(argvNew[1],"test/Test_L2_proj_errors_TRI");

	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	PMin = 0;  PMax = 4;
	MLMin = 0; MLMax = 4;

	for (P = PMin; P <= PMax; P++) {
	for (ML = MLMin; ML <= MLMax; ML++) {
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		code_startup(nargc,(char const *const *const) argvNew,0,1);
		compute_errors_global();
		code_cleanup();
	}}

pass = 0;
printf("pass: %d\n",pass);


	free(argvNew[0]); free(argvNew[1]); free(argvNew);
}
