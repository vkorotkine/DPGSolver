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
	if (!DB.MPIrank) printf("Preprocessing:\n\n");

	if (!DB.MPIrank) printf("  Set up Parameters\n");
	setup_parameters();

	if (!DB.MPIrank) printf("  Set up Mesh\n");
	setup_mesh();

	if (!DB.MPIrank) printf("  Set up Operators\n");
	setup_operators();

	if (!DB.MPIrank) printf("  Set up Structures\n");
	setup_structures();

	if (!DB.MPIrank) printf("  Set up Geometry\n");
	setup_geometry();


	memory_free();

	// End MPI and PETSC
	PetscFinalize();

	return 0;
}
