// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "test_integration_linearization.h"

#include <string.h>

#include "Simulation.h"
#include "Intrusive.h"
#include "Element.h"
#include "Mesh.h"
#include "Volume.h"
#include "Face.h"

#include "Macros.h"

void test_integration_linearization (const char*const ctrl_name)
{
	struct Simulation*const simulation = constructor_Simulation(); // destructed

	set_simulation_core(simulation,ctrl_name);
// set up parameters when needed (geometry/solver/postprocessing)

	set_Simulation_elements(simulation,constructor_Element_List(simulation->d));

	struct Mesh* mesh = set_up_mesh(simulation->mesh_name_full,simulation->d,simulation->elements);

	simulation->volumes = constructor_Volume_List(simulation,mesh);
	simulation->faces   = constructor_Face_List(simulation,mesh);

	destructor_Mesh(mesh);
	destructor_Simulation(simulation);
}
