// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
#include "solver_implicit.h"
#include "solver_Poisson.h"
#include "compute_errors.h"
#include "memory_free.h"

/*

NOTES:
	- The ordering used is the same as the nodes for the elements. That is,
		the node at 0 is the one at -1,-1 on the reference space, and the third 
		node in the msh file will be at the 1,1 point on the reference space.
		Using this knowledge, the control points were ordered accordingly when
		creating the connecitivity for each volume.


File Name Convention: (For test files)
This holds the naming convention to be used for the ctrl and mesh files for the Bezier and
NURBS meshes. The difficulty primarily is that the order P is needed for the Bezier mesh
and so this has to be in the mesh file name. For this reason, the ctrl file will not have the 
P value in its name but the mesh files will have the P values in them.

- Mesh File:
	n-GaussianBump2D_ToBeCurved_BezierP3_QUAD0x.msh

	Build this name in the initialization using the control file
	and the basis type.


To Do:

- Create a Python B Spline Mesh Generator
	- Refine the existing algorithm, and make this a much better generator
		that should be able to handle the case of multiple B Spline patches
	- Add the Bezier Extraction step using Numpy properly.
		- Research what the forms of the operators will be to make it easiest to
			do. 
		- Use this now to generate the meshes. That is, create a general mesh 
			without knot repetition lets say, and perform the refinement to create the
			Bezier elements.

	- This should be a seperate module that will take the elements and output the 
		Gmsh file in the proper format.
		- Will only need lines and Quads
		- Also, BCs are simple. 
		- Check if this Gmsh file works in the code.
			- Use the CurvedQUAD keyword for the mesh. Then, when we need to use the Bezier
				basis, we will read 
			- Add a keyword $NURBS_Weights$ and place the weights of all node points in 
				this section. (For now, use B-Splines so all weights are 1)
			- Add a keyword $NURBS_ControlPoints$ which will contain the xyz_coeff
				values for each element (have used Bezier extraction so each element 
				can be treated separately with their own control points).

- Find issue with optimal order convergence
	- Try the Bezier mesh using the polynomial basis functions. See the convergence order
		for this case.
		- This will check whether the cause of the issues is the mesh itself.

	(Done) - Try Bezier Mesh with very small curvature
		- Orders still do not converge. Since orders converged for the Periodic Vortex
			case but not for this there must be an issue with the mesh when it is
			generated.
	
	(Done) - Try the Bezier basis functions using the periodic vortex case
		(Done) - Check the periodic vortex order convergence (mesh should not be an issue
			since the mesh is rectangular).
			- Used same mesh in the code with HP refinement (even though using Bezier
				basis). Orders worked out.
				- Therefore, the basis has been implemented perfectly
	
	(Done) - Try CURVED mesh (Gaussian Bump) with Bezier basis
		- Orders passed for the Bezier basis. Used the HP adaptation in place with the Bezier 
			basis to do this test.


- Create the NURBS Mesh Generator
	- First work with one patch cases
		(Done) - Create the a simply mesh with no Bezier Extraction
		(Done) - Use weights on the first mesh in order to be able to test the NURBS
			capabilities.
		- Use Bezier extraction to extract the Bezier volumes from the simple mesh
		- Add h and p refinement (on an element level after the extraction process)
			in order to be able to add more volumes
		- When plotting, plot using the proper NURBS definition (do not just connect
			points).
		- Add Boundary Condition capabilities easily like in the other mesh
			case.
			- Set the edges (set i,j range) and set the BC on this range
		- Output now Control Points and Weights for each control point also.

	(Done) - Create a semi-circular bump case
		(Done) - The unrefined mesh will be able to set the outer lines for the 
			mesh using a set of control points (i.e. n knots along x direction
			used to define the curve and maybe just 2 knot in y direction for 
			the y case (only top and bottom))
			- Then, perform uniform refinement to be able to create the mesh 
				levels.

	- Handle multiple patch case
		- Define element connectivity for the unrefined patch (knot sequence should
			be same between connected faces and specify which face are connected).
		- Refinement will need to be uniform for all patches
			- Once refined, set which elements are connected.

- Create the NURBS basis functions in own code
	- Work first with the B-spline mesh with all weights as one (validation
		will be done to make sure all outputs are identical as before for 
		the explicit form of the equations).
		- Add the option to read in the weights now from the mesh file 
			for the IGA-DG code. Each control point will now have a weight 
			associated with it.	
		- Build a special operator called W_vI which is the weight function
			evaluated at the integration nodes. This should not change unless
			adaptivity is used. 
		-  


 - Create a tecplot output for visualizing the solution using zones.
 	Do this because tecplot is easier to use
 


*/












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
	initialize_test_case_parameters();
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

	if (strstr(DB.TestCase,"Poisson")) {
		solver_Poisson(1);
	} else {
		if (strstr(DB.SolverType,"Explicit")) {
			solver_explicit(1);
		} else if (strstr(DB.SolverType,"Implicit")) {
			solver_implicit(1);
		} else {
			printf("Error: Unsupported SolverType in dpg_solver.\n"), EXIT_MSG;
		}
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
	DB.TestL2projection = 0;
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
#include <stdbool.h>

#include "petscsys.h"

#include "Macros.h"

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
#include "test_regression_fluxes_inviscid.h"
#include "test_unit_jacobian_fluxes.h"
#include "test_unit_jacobian_boundary.h"
#include "test_unit_get_face_ordering.h"
#include "test_unit_equivalence_real_complex.h"

#include "test_integration_update_h.h"
#include "test_integration_L2_projections.h"
#include "test_integration_Advection.h"
#include "test_integration_Poisson.h"
#include "test_integration_Euler.h"
#include "test_integration_NavierStokes.h"

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
 *		Linearizations are tested using the complex step method (Squire(1998), Martins(2003)).
 *
 *	Notation:
 *
 *	References:
 *		Squire(1998)-Using Complex Variables to Estimate Derivatives of Real Functions
 *		Martins(2003)-The Complex-Step Derivative Approximation
 */

int main(int nargc, char **argv)
{
	struct S_RunTest {
		unsigned int unit, integration, speed;
	} RunTest;

	clock_t ts, te;

	TestDB.Active = 1;
	TestDB.Ntest = 0;
	TestDB.Npass = 0;
	TestDB.Nwarnings = 0;

	RunTest.unit        = 0;
	RunTest.integration = 0;
	RunTest.speed       = 0;

	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);

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

		test_unit_matrix_functions();

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

		test_regression_fluxes_inviscid();
		test_unit_jacobian_fluxes();
		test_unit_jacobian_boundary();
		test_unit_get_face_ordering();

		test_unit_equivalence_real_complex();
	}

	// Integration tests
	if (RunTest.integration) {
		test_integration_update_h(nargc,argv);
		test_integration_L2_projections(nargc,argv);
	}
if (1) {
//	test_integration_Advection(nargc,argv);
//	test_integration_Poisson(nargc,argv);
	test_integration_Euler(nargc,argv);
//	test_integration_NavierStokes(nargc,argv);
}

	PetscFinalize();

	te = clock();

	// Speed tests
	if (RunTest.speed) {
		test_speed_array_swap();
		test_speed_mm_CTN();
	}


	printf("\n\nRan %d test(s) in %.4f seconds.\n",TestDB.Ntest,(te-ts)/(float)CLOCKS_PER_SEC);

	int Nfail = TestDB.Ntest - TestDB.Npass;
	if (Nfail != 0) {
		printf("\n******** FAILED %d TEST(S) ********\n\n",Nfail);
	} else {
		printf("\nAll tests passed.\n\n");

		if (TestDB.Nwarnings)
			printf("Warnings (%d) were generated while running tests. "
			       "Scroll through test passing list and verify that all is OK.\n\n",TestDB.Nwarnings);
	}

	return 0;
}


#endif // End TEST
