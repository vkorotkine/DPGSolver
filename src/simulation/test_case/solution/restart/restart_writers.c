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
///	\file

#include "restart_writers.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "definitions_mesh.h"

#include "multiarray.h"

#include "volume_solver.h"

#include "file_processing.h"
#include "intrusive.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

///\{ \name Supported restart file types.
#define RESTART_GMSH 1
///\}

/// \brief Write a restart file with the mesh component in gmsh format.
static void output_restart_gmsh
	(const struct Simulation*const sim,           ///< \ref Simulation.
	 const struct Restart_Info*const restart_info ///< \ref Restart_Info.
	);

// Interface functions ********************************************************************************************** //

void output_restart (const struct Simulation*const sim, const struct Restart_Info*const restart_info)
{
	switch (sim->domain_type) {
	case DOM_STRAIGHT:
	case DOM_PARAMETRIC:
		break; // Do nothing
	default:
		// Only affine restart files are currently supported.
/// \todo Think about whether blended domains could fit in this framework as only the boundary volumes are curved.
		EXIT_ERROR("Unsupported: %d\n",sim->domain_type);
		break;
	}

	if (strstr(sim->mesh_name_full,".msh"))
		output_restart_gmsh(sim,restart_info);
	else
		EXIT_ERROR("Unsupported: %s\n",sim->mesh_name_full);

EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for information related to the restart mesh nodes.
struct Nodes_Info {
	const struct const_Multiarray_d* nodes; ///< The xyz node coordinates.
};

/** \brief Constructor for the \ref Nodes_Info container.
 *  \return See brief. */
static const struct Nodes_Info* constructor_Nodes_Info
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Destructor for the \ref Nodes_Info container.
static void destructor_Nodes_Info
	(const struct Nodes_Info*const nodes_info ///< Standard.
	);

/** \brief Return the name of the restart file.
 *  \return See brief. */
static const char* compute_restart_name
	(const int restart_type,            ///< The type of restart file. See options above.
	 const struct Simulation*const sim, ///< \ref Simulation.
	 const struct Restart_Info*const ri ///< \ref Restart_Info.
	);

static void output_restart_gmsh (const struct Simulation*const sim, const struct Restart_Info*const restart_info)
{
	const struct Nodes_Info*const ni = constructor_Nodes_Info(sim); // destructed

	const char*const f_name = compute_restart_name(RESTART_GMSH,sim,restart_info);

	FILE* file = fopen_create_dir(f_name); // closed
	fprintf(file,"$MeshFormat");
	fprintf(file,"2.2 0 8");
	fprintf(file,"$EndMeshFormat");

	destructor_Nodes_Info(ni);
EXIT_UNSUPPORTED;
	fclose(file);
}

// Level 1 ********************************************************************************************************** //

static const struct Nodes_Info* constructor_Nodes_Info (const struct Simulation*const sim)
{
	struct Nodes_Info* nodes_info = calloc(1,sizeof(*nodes_info)); // free

	struct Multiarray_d*const xyz_ve_red = constructor_empty_Multiarray_d('R',2,(ptrdiff_t[]){0,DIM}); // destructed
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Volume*const vol = (struct Volume*) curr;

		push_back_Multiarray_d(xyz_ve_red,vol->xyz_ve);
print_const_Multiarray_d(vol->xyz_ve);
	}
print_Multiarray_d(xyz_ve_red);
	/// \todo sort by approximate nearest neighbor.

	destructor_Multiarray_d(xyz_ve_red);

//print_const_Multiarray_d(ni->nodes);
EXIT_UNSUPPORTED; UNUSED(sim);
	return nodes_info;
}

static void destructor_Nodes_Info (const struct Nodes_Info*const nodes_info)
{
	destructor_const_Multiarray_d(nodes_info->nodes);
	free((void*)nodes_info);
}

static void destructor_Nodes_Info
	(const struct Nodes_Info*const nodes_info ///< Standard.
	);
static const char* compute_restart_name
	(const int restart_type, const struct Simulation*const sim, const struct Restart_Info*const ri)
{
	static const char* name_path = "../restart/";

	char ml_p_spec[STRLEN_MIN];
	sprintf(ml_p_spec,"%s%d%s%d","__ml",ri->ml,"__p",ri->p);

	char name_ext[STRLEN_MIN];
	switch (restart_type) {
		case RESTART_GMSH: strcpy(name_ext,".msh");                      break;
		default:           EXIT_ERROR("Unsupported: %d\n",restart_type); break;
	}

	static char output_name[5*STRLEN_MAX];
	sprintf(output_name,"%s%s%c%s%c%s%s%s",
	        name_path,sim->pde_name,'/',sim->pde_spec,'/',"restart",ml_p_spec,name_ext);
	return output_name;
}
