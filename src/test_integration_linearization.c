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
	struct Simulation*const sim = constructor_Simulation(); // destructed

	set_simulation_core(sim,ctrl_name);
// set up parameters when needed (geometry/solver/postprocessing)

	set_Simulation_elements(sim,constructor_Element_List(sim->d));

	struct Mesh_Input mesh_input = { .d              = sim->d,
	                                 .domain_type    = sim->domain_type,
	                                 .mesh_name_full = sim->mesh_name_full,
	                                 .geom_name      = sim->geom_name,
	                                 .geom_spec      = sim->geom_spec,
	                                 .input_path     = sim->input_path, };

	struct Mesh* mesh = set_up_mesh(&mesh_input,sim->elements);

	sim->volumes = constructor_Volume_List(sim,mesh);
	sim->faces   = constructor_Face_List(sim,mesh);

	destructor_Mesh(mesh);




	destructor_Simulation(sim);
}
