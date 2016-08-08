// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_mesh.h"

/*
 *	Purpose:
 *		Set up mesh related parameters.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void setup_mesh()
{
	initialize_ELEMENTs();

	// Read mesh file
	if (!DB.MPIrank && !TEST)
		printf("    Read MeshFile\n");
	gmsh_reader();

	finalize_ELEMENTs();

	// Build Connectivity
	if (!DB.MPIrank && !TEST)
		printf("    Set up connectivity\n");
	setup_connectivity();

	// Modify connectivity if periodic
	if (DB.NPVe != 0) {
		if (!DB.MPIrank && !TEST)
			printf("    Modify connectivity for periodic\n");
		setup_periodic();
	}
}
