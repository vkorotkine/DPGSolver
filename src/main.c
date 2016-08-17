// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "main.h"

#include "S_DB.h"
#include "Test.h"

struct S_DB   DB;
struct S_TEST TestDB;


#ifndef TEST

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h> // ToBeModified: Likely not use system headers for mpi/petsc
#include <petscksp.h>
 
#include "Parameters.h"
#include "Macros.h"

#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "output_to_paraview.h"
#include "solver_explicit.h"
//#include "solver_implicit.h"
#include "compute_errors.h"
#include "memory_free.h"

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

int main(int nargc, char **argv)
{
	printf("\n\n\n*** Test to see when unrolled mv multiplications break even with BLAS on Guillimin before running"
	       " any large jobs. ***\n\n\n");
	printf("\n\n\n*** Test to see whether the use of floats initially then transferring to doubles results in speed-up."
	       " ***\n\n\n");
	printf("\n\n\n*** Combine explicit and implicit info for implicit runs (ToBeDeleted) ***\n\n\n");
	// Note: The same recombination can be done for the flux and boundary condition functions as well (ToBeDeleted).

	char *fNameOut, *string;
	int  MPIrank, MPIsize;

	struct S_TIME {
		clock_t ts, te;
		double  tt;
	} total, preproc, solving, postproc;

	total.ts = clock();

	fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut); // free
	string   = malloc(STRLEN_MIN * sizeof *string);   // free

	// Start MPI and PETSC
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	// Test memory leaks only from Petsc and MPI using valgrind
	//PetscFinalize(), EXIT_MSG;

	DB.MPIsize = MPIsize;
	DB.MPIrank = MPIrank;

	// Initialization
	preproc.ts = clock();
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

	preproc.te = clock();

	// Solving
	solving.ts = clock();
	if (!DB.MPIrank)
		printf("\n\nSolving:\n\n");

	if (!DB.MPIrank)
		printf("  Initializing\n");
	initialize_test_case(DB.LevelsMax+1);

	// Output initial solution to paraview
	strcpy(fNameOut,"SolInitial_");
	                             strcat(fNameOut,DB.TestCase);
	sprintf(string,"%dD",DB.d); strcat(fNameOut,string);
	output_to_paraview(fNameOut);

	if (DB.Restart >= 0) {
		if (!DB.MPIrank)
			printf("  Initializing restarted solution if enabled.\n");
		// Need to ensure that the same proc distribution is used as for lower order solutions.
//		restart_read();
	}

	if (!DB.MPIrank)
		printf("  Nonlinear Iterative Solve\n\n");

	if (strstr(DB.SolverType,"Explicit")) {
		solver_explicit();
	} else if (strstr(DB.SolverType,"Implicit")) {
		; //solver_implicit();
	} else {
		printf("Error: Unsupported SolverType in dpg_solver.\n"), EXIT_MSG;
	}
	solving.te = clock();

	// Postprocessing
	postproc.ts = clock();
	if (!DB.MPIrank)
		printf("\n\nPostprocessing:\n\n");

	// Output final solution to paraview
	printf("  Output final solution to paraview\n");
	strcpy(fNameOut,"SolFinal_");
	sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
	                               strcat(fNameOut,DB.MeshType);
	sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
	if (DB.Adapt == ADAPT_0)
		sprintf(string,"P%d_",DB.PGlobal), strcat(fNameOut,string);
	output_to_paraview(fNameOut);

	// Compute errors
	if (!DB.MPIrank)
		printf("  Computing errors\n");
	compute_errors_global();

	postproc.te = clock();

	free(fNameOut);
	free(string);

	memory_free();

	// End MPI and PETSC
	PetscFinalize();

	total.te = clock();

	printf("\n\n\nTotal time       : % .2f s\n\n",(total.te-total.ts)/(double) CLOCKS_PER_SEC);
	printf("  Preprocessing  : % .2f s\n",(preproc.te-preproc.ts)/(double) CLOCKS_PER_SEC);
	printf("  Solving        : % .2f s\n",(solving.te-solving.ts)/(double) CLOCKS_PER_SEC);
	printf("  Postprocessing : % .2f s\n",(postproc.te-postproc.ts)/(double) CLOCKS_PER_SEC);
	printf("\n\n\n");

	return 0;
}

#else // Run if -DTEST is passed as a compilation flag

#include <stdio.h>
#include <time.h>

#include "test_unit_array_find_index.h"
#include "test_unit_array_norm.h"
#include "test_unit_array_sort.h"
#include "test_unit_array_swap.h"
#include "test_unit_math_functions.h"
#include "test_unit_matrix_functions.h"
#include "test_unit_bases.h"
#include "test_unit_grad_bases.h"
#include "test_unit_cubature.h"
#include "test_unit_find_periodic_connections.h"
#include "test_unit_sum_factorization.h"
#include "test_unit_plotting.h"
#include "test_unit_fluxes_inviscid.h"
#include "test_unit_jacobian_fluxes_inviscid.h"
#include "test_unit_jacobian_boundary.h"
#include "test_unit_get_facet_ordering.h"
#include "test_unit_equivalence_real_complex.h"

#include "test_integration_L2_projections.h"
#include "test_integration_update_h.h"
#include "test_integration_linearization.h"

#include "test_speed_array_swap.h"
#include "test_speed_mm_CTN.h"

/*
 *	Purpose:
 *		Run test functions:
 *			1) Speed comparisons
 *			2) Correctness of implementation (Individual functions as well as overall code)
 *
 *	Comments:
 *		Get some kind of code coverage figure as well (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

int main(int nargc, char **argv)
{
	struct S_RunTest {
		unsigned int unit, integration, speed;
	} RunTest;

	clock_t ts, te;

	TestDB.Ntest = 0;
	TestDB.Npass = 0;
	TestDB.Nwarnings = 0;

	RunTest.unit        = 1;
	RunTest.integration = 1;
	RunTest.speed       = 0;


	printf("\n\nRunning Tests:\n\n\n");
	ts = clock();

	// Unit tests
	if (RunTest.unit) {
		test_unit_array_find_index();
		test_unit_array_norm();
		test_unit_array_sort();
		test_unit_array_swap();

		test_unit_math_factorial();
		test_unit_math_gamma();

		test_unit_matrix_diag();
		test_unit_matrix_identity();
		test_unit_matrix_inverse();
		test_unit_matrix_mm();
		test_unit_convert_to_CSR();

		test_unit_find_periodic_connections();

		test_unit_cubature_TP();
		test_unit_cubature_SI();
		test_unit_cubature_PYR();

		test_unit_basis_TP();
		test_unit_basis_SI();
		test_unit_basis_PYR();
		test_unit_grad_basis_TP();
		test_unit_grad_basis_SI();
		test_unit_grad_basis_PYR();

		test_unit_sum_factorization();
		test_unit_plotting();

		test_unit_fluxes_inviscid();
		test_unit_jacobian_fluxes_inviscid();
		test_unit_jacobian_boundary();
		test_unit_get_facet_ordering();

		test_unit_equivalence_real_complex();
		printf("\nAdd tests for real/complex equivalence for all relevant functions.\n\n");
		printf("\nFor the VOLUME/FACET info functions, test that all 'versions' give identical results.\n\n");
		TestDB.Nwarnings++;
	}
test_unit_jacobian_fluxes_inviscid();

	// Integration tests
	if (RunTest.integration) {
		test_integration_update_h(nargc,argv);
		test_integration_L2_projections(nargc,argv);
	}
	test_integration_linearization(nargc,argv);


	te = clock();

	// Speed tests
	if (RunTest.speed) {
		test_speed_array_swap();
		test_speed_mm_CTN();
	}


	printf("\n\nRan %d test(s) in %.4f seconds.\n",TestDB.Ntest,(te-ts)/(float)CLOCKS_PER_SEC);

	unsigned int Nfail = TestDB.Ntest - TestDB.Npass;
	if (Nfail > 0) {
		printf("\n\n******** FAILED %d TEST(S) ********\n\n",Nfail);
	} else {
		printf("\nAll tests passed.\n\n");

		if (TestDB.Nwarnings)
			printf("Warnings (%d) were generated while running tests. "
			       "Scroll through test passing list and verify that all is OK.\n\n",TestDB.Nwarnings);
	}

	return 0;
}


#endif // End TEST
