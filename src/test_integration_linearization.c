// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_linearization.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "array_norm.h"
#include "implicit_VOLUME_info.h"
#include "finalize_LHS.h"

/*
 *	Purpose:
 *		Test correctness of implementation of code linearization.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_linearization(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *
	 *		Mesh file on which to run the test.
	 *
	 *	Expected Output:
	 *
	 *		The analytical linearization should match that computed using complex step.
	 *
	 */

//	PetscInt m, n, M, N, d_nz, *d_nnz, o_nz, *o_nnz;
	Mat      A = NULL;
	Vec      b = NULL;
//	MPI_Comm comm = MPI_COMM_WORLD;

	// **************************************************************************************************** //
	// LINEs


	// **************************************************************************************************** //
	// TRIs
	strcpy(argvNew[1],"test/Test_linearization_TRI");

	code_startup(nargc,argvNew,2);

	implicit_VOLUME_info();
	compute_dof();
	finalize_LHS(A,b,1);



//	if(CHKERRQ(MatCreateSeqAIJ(comm,m,n,0,nnz[],*A))) EXIT_MSG;
	if (A == NULL)
		printf("A NULL\n");

pass = nargc;
printf("%d\n",pass);

	code_cleanup(0);


	free(argvNew[0]); free(argvNew[1]); free(argvNew);
}
