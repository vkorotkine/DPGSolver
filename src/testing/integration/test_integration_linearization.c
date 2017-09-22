/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/**	\file
 */

#include "test_integration_linearization.h"

#include <string.h>

#include "macros.h"

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

	sim->volumes = constructor_Volumes(sim,mesh);
	sim->faces   = constructor_Faces(sim,mesh);

	destructor_Mesh(mesh);




	destructor_Simulation(sim);
}
