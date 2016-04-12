#ifndef TEST

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petscksp.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Solve the (N)avier-(S)tokes equations (or a subset of the NS equations) using the (D)iscontinuous
 *		(P)etrov-(G)alerkin method.
 *
 *	Comments:
 *
 *	Notation:
 *
 *		DB : (D)ata(B)ase
 *
 *	References:
 *		Demkowicz(2010)_A class of discontinuous Petrov–Galerkin methods. Part I - The transport equation
 *		ToBeModified: Add significant references.
 *
 */

struct S_DB DB;

int main(int nargc, char **argv)
{
	int MPIrank, MPIsize;

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

	// Preprocessing
	if (!DB.MPIrank)
		printf("Preprocessing:\n\n");

	if (!DB.MPIrank)
		printf("  Set up Parameters\n");
	setup_parameters();

	if (!DB.MPIrank)
		printf("  Set up Mesh\n");
	setup_mesh();

	if (!DB.MPIrank)
		printf("  Set up Operators\n");
	setup_operators();

	if (!DB.MPIrank)
		printf("  Set up Structures\n");
	setup_structures();

	if (!DB.MPIrank)
		printf("  Set up Geometry\n");
	setup_geometry();


	memory_free();

	// End MPI and PETSC
	PetscFinalize();

	return 0;
}

#else // Run if -DTEST is passed as a compilation flag

#include <stdio.h>
#include <time.h>

#include "test.h"
#include "database.h"
#include "functions.h"

/*
 *	Purpose:
 *		Run test functions:
 *			1) Speed comparisons
 *			2) Correctness of implementation
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_DB DB;
struct S_TEST TestDB;

#define TOL 1e-15

int main(void)
{
	clock_t ts, te;

	TestDB.Ntest = 0;
	TestDB.Npass = 0;
	TestDB.Nwarnings = 0;


	printf("\n\nRunning Tests:\n\n\n");
	ts = clock();

	// Implementation tests
	test_imp_array_find_index();
	test_imp_array_norm();
	test_imp_array_sort();
	test_imp_array_swap();

	test_imp_math_factorial();
	test_imp_math_gamma();

	test_imp_matrix_diag();
	test_imp_matrix_identity();
	test_imp_matrix_inverse();
	test_imp_matrix_mm();

	test_imp_find_periodic_connections();

	test_imp_cubature_TP();
	test_imp_basis_TP();
	test_imp_grad_basis_TP();

	// Speed tests
	test_speed_array_swap();
//	test_speed_mm_d();


	te = clock();


	printf("\n\nRan %d test(s) in %.4f seconds.\n",TestDB.Ntest,(te-ts)/(float)CLOCKS_PER_SEC);

	unsigned int Nfail = TestDB.Ntest - TestDB.Npass;
	if (Nfail > 0) {
		printf("\n\n******** FAILED %d TEST(S) ********\n\n",Nfail);
	} else {
		printf("\nAll tests passed.\n\n");

		if (TestDB.Nwarnings)
			if (TestDB.Nwarnings == 1)
			printf("Warnings (%d) were generated while running tests. "
				   "Scroll through test passing list and verify that all is OK.\n\n",TestDB.Nwarnings);
	}


	return 0;
}


#endif // End TEST
