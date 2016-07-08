// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <petscksp.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"
#include "database.h"

/*
 *	Purpose:
 *		Test correctness of implementation of update_connectivity.
 *
 *	Comments:
 *		Ensure to test with BCs other than periodic as well. (ToBeDeleted)
 *		Fix memory leaks. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void code_startup(int nargc, char **argv)
{
	int  MPIrank, MPIsize;

	// Start MPI and PETSC
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	// Test memory leaks only from Petsc and MPI using valgrind
	//PetscFinalize(), exit(1);

	DB.MPIsize = MPIsize;
	DB.MPIrank = MPIrank;

	// Initialization
	initialization(nargc,argv);
	setup_parameters();
	setup_mesh();
	setup_operators();
	setup_structures();
}

void test_imp_update_connectivity(int nargc, char **argv)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		ToBeModified
	 *
	 *	Expected Output (d = 1, Trivial):
	 *
	 *		ToBeModified
	 *
	 */

	char **argvNew;

	// LINEs
// Start with d = 1? Can ensure that the code set up is working for this case at the same time.


	// TRIs
	argvNew = malloc(2 * sizeof *argvNew);
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew);
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew);

	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/TestTRI");

	for (int i = 0; i < 2; i++)
		printf("%s\n",argvNew[i]);

	code_startup(nargc,argvNew);



exit(1);

/*
	pass = 0;
	if (array_norm_diff_ui(Nn,nOrd,nOrd10,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("get_facet_ordering (d = 1, case 0):              ");
	test_print(pass);
*/

}
