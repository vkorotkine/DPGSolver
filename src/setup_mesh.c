// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "setup_mesh.h"

#include <stdlib.h>
#include <stdio.h>

#include "S_DB.h"

#include "element_functions.h"
#include "gmsh_reader.h"
#include "setup_connectivity.h"
#include "setup_periodic.h"

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
// Move this to be immediately adjacent to finalize_ELEMENTs when finished (ToBeDeleted)
	initialize_ELEMENTs();

	// Read mesh file
	if (!DB.MPIrank && !DB.Testing)
		printf("    Read MeshFile\n");
	gmsh_reader();

	finalize_ELEMENTs();

	// Build Connectivity
	if (!DB.MPIrank && !DB.Testing)
		printf("    Set up connectivity\n");
	setup_connectivity();

	// Modify connectivity if periodic
	if (DB.NPVe != 0) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    Modify connectivity for periodic\n");
		setup_periodic();
	}
}
