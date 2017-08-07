// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "test_integration_linearization.h"

#include <string.h>

#include "Simulation.h"

void test_integration_linearization (const char*const ctrl_name)
{
	struct Simulation*const simulation = constructor_Simulation(); // destructed

	set_simulation_core(simulation,ctrl_name);
	// set up parameters when needed (geometry/solver/postprocessing)
	set_up_mesh(simulation);
	// re-write gmsh reader such that global variables are not used. Also do not include VeCGmsh in Element struct.


	destructor_Simulation(simulation);
}
