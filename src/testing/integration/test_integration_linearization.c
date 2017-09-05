// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_integration_linearization.h"

#include <string.h>

#include "Macros.h"

#include "simulation.h"
#include "intrusive.h"
#include "element.h"
#include "mesh.h"
#include "volume.h"
#include "face.h"

void test_integration_linearization (const char*const ctrl_name)
{
	struct Simulation*const sim = constructor_Simulation(ctrl_name); // destructed
// set up parameters when needed (geometry/solver/postprocessing)

	struct Mesh_Input mesh_input = set_Mesh_Input(sim);
	struct Mesh* mesh = constructor_Mesh(&mesh_input,sim->elements);

	sim->volumes = constructor_Volume_List(sim,mesh);
	sim->faces   = constructor_Face_List(sim,mesh);

	destructor_Mesh(mesh);




	destructor_Simulation(sim);
}
