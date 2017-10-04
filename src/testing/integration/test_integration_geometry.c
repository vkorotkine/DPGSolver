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
/** \file
 */

#include "test_integration_geometry.h"

#include <stdlib.h>
#include <string.h>

#include "test_base.h"

#include "macros.h"

#include "simulation.h"
#include "mesh.h"
#include "volume.h"
#include "solver_volume.h"
#include "face.h"
#include "file_processing.h"
#include "geometry.h"


// Static function declarations ************************************************************************************* //

/** \brief Compare \ref Solver_Volume::geom_coef and \todo [ref here] Face::normal_f_i finite members with their expected values.
 *  \return `true` if tests passed. */
static bool compare_members_geom
	(struct Test_Info*const test_info, ///< Defined in \ref test_integration_mesh.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void test_integration_geometry (struct Test_Info*const test_info, const char*const ctrl_name)
{
	struct Simulation*const sim = constructor_Simulation(ctrl_name); // destructed

	struct Mesh_Input mesh_input = set_Mesh_Input(sim);
	struct Mesh* mesh = constructor_Mesh(&mesh_input,sim->elements);

	sim->volumes = constructor_Volumes(sim,mesh);
	sim->faces   = constructor_Faces(sim,mesh);

	destructor_Mesh(mesh);

	sim->volumes = constructor_Solver_Volumes(sim);
	sim->faces   = constructor_Solver_Faces(sim);

	set_up_solver_geometry(sim);


	const bool pass = compare_members_geom(test_info,sim);

	char* test_name_end = extract_name(ctrl_name,false); // free

	char test_name[STRLEN_MAX];
	strcpy(test_name,"Geom initialization - ");
	strcat(test_name,test_name_end);

	free(test_name_end);

	test_increment_and_print(test_info,pass,test_name);

	destructor_Simulation(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static bool compare_members_geom (struct Test_Info*const test_info, const struct Simulation*const sim)
{
	bool pass = true;
UNUSED(test_info);
UNUSED(sim);

pass = false;

//set up volume xyz coordinates and normal vectors
//	compute metrics using:
//check free-stream preservation.
//start with the wavy periodic vortex case.
//	add blending and verify implementation by correcting the wavy periodic vortex; check:
//		edges:   no change (computed from parametrized map exactly as in the standard mapping)
//		faces:   no change (computed from parametrized map exactly as in the standard mapping)
//		volumes: minor change in position of internal nodes.
//also try with polynomial geometry based on bezier patch/patches.


	return pass;
}
