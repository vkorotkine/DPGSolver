// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_code_integration.h"

#include <string.h>

#include "petscsys.h"

#include "Parameters.h"
#include "Test.h"
#include "S_DB.h"

#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "adaptation.h"
#include "memory_free.h"

/*
 *	Purpose:
 *		Provide functions for integration testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void update_MeshFile(void)
{
	// Standard datatypes
	char *d, *ML;

	d  = malloc(STRLEN_MIN * sizeof *d);  // free
	ML = malloc(STRLEN_MIN * sizeof *ML); // free

	sprintf(d,"%d",DB.d);
	sprintf(ML,"%d",DB.ML);

	strcpy(DB.MeshFile,"");
	strcat(DB.MeshFile,DB.MeshPath);
	strcat(DB.MeshFile,DB.TestCase);
	strcat(DB.MeshFile,"/");
	strcat(DB.MeshFile,DB.TestCase);
	strcat(DB.MeshFile,strcat(d,"D_"));
	strcat(DB.MeshFile,DB.MeshType);
	strcat(DB.MeshFile,strcat(ML,"x.msh"));

	free(d);
	free(ML);
}

void code_startup(int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv)
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
	if (update_argv) {
		strcpy(DB.TestCase,TestDB.TestCase);
		DB.PGlobal = TestDB.PGlobal;
		DB.ML      = TestDB.ML;
		update_MeshFile();
	}

	setup_parameters();
	if (update_argv)
		setup_parameters_L2proj();

	setup_mesh();
	setup_operators();
	setup_structures();
	setup_geometry();

	initialize_test_case(Nref);
}

void code_cleanup(void)
{
	mesh_to_level(0);
	memory_free();
}
