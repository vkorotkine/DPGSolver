// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "test_integration_linearization.h"

#include <string.h>

#include "Simulation.h"
#include "Intrusive.h"
#include "Element.h"
#include "set_up_mesh.h"

#include "Macros.h"

void test_integration_linearization (const char*const ctrl_name)
{
	struct Simulation*const simulation = constructor_Simulation(); // destructed

	set_simulation_core(simulation,ctrl_name);
	// set up parameters when needed (geometry/solver/postprocessing)

	struct Intrusive_List* elements = constructor_Element_List(simulation->d);

	struct Mesh* mesh = set_up_mesh(simulation->mesh_name_full,simulation->d,elements);
DO_NOTHING_P(mesh);

	destructor_Simulation(simulation);
}
