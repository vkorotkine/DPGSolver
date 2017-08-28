// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "simulation.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#include "Macros.h"
#include "element.h"
#include "volume.h"
#include "face.h"

#include "file_processing.h"
#include "const_cast.h"

#include "constants_mesh.h"

// Static function declarations ************************************************************************************* //

///< \brief Holds data relating to the mesh as read from the control file.
struct Mesh_Ctrl_Data {
	char mesh_path[STRLEN_MAX],      ///< Relative path to the meshes directory.
	     mesh_elem_type[STRLEN_MIN], ///< Type of elements present in the mesh; used to set the mesh file name.
	     mesh_curving[STRLEN_MAX],   /**< Type of the mesh.
	                                  *   	Options: Straight, Curved, Parametric. */
	     mesh_generator[STRLEN_MAX], /**< Type of mesh generator used.
	                                  *   	Options: gmsh */
	     mesh_extension[STRLEN_MIN], ///< File extension (set based on \ref mesh_generator.
	     mesh_domain[STRLEN_MIN];    ///< The domain type. \todo Should replace mesh_curving.
};

/// \brief Set full control file name (including path and file extension).
static void set_ctrl_name_full
	(struct Simulation*const sim, /// Standard.
	 const char*const ctrl_name   /// Defined in \ref set_simulation_core.
	);

/// \brief Set the mesh parameters.
static void set_mesh_parameters
	(struct Simulation*const sim ///< Standard.
	);

/** \brief Set associations between `char*` and `int` variables.
 *
 *	This is done such that if/switch conditions are simplified when these variables are used.
 *
 *	\todo This should be moved to the Solver context.
 */
static void set_string_associations
	(struct Simulation*const sim ///< Standard.
	);

/// \brief Set the path to the input files.
static void set_input_path
	(struct Simulation*const sim ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Simulation* constructor_Simulation ()
{
	struct Simulation* sim = malloc(sizeof *sim); // returned;

	return sim;
}

void destructor_Simulation (struct Simulation* sim)
{
	destructor_Elements((struct Intrusive_List*) sim->elements);
	destructor_Volumes(sim->volumes);
	destructor_Faces(sim->faces);
	FREE_NULL(sim);
}

void set_simulation_core (struct Simulation*const sim, const char*const ctrl_name)
{
/// \todo Remove commented variables below.

	set_ctrl_name_full(sim,ctrl_name);
	FILE *ctrl_file = fopen_checked(sim->ctrl_name_full);

	// Read information
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"PDEName"))  read_skip_const_c(line,sim->pde_name);
		if (strstr(line,"PDESpec"))  read_skip_const_c(line,sim->pde_spec);
		if (strstr(line,"Geometry")) read_skip_const_c(line,sim->geom_name);
		if (strstr(line,"GeomSpec")) read_skip_const_c(line,sim->geom_spec);

		if (strstr(line,"MeshLevel")) read_skip_const_i(line,&sim->ml);

//		if (strstr(line,"NodeType"))  read_skip_const_c(line,sim->node_type);
//		if (strstr(line,"BasisType")) read_skip_const_c(line,sim->basis_type);

//		if (strstr(line,"Vectorized")) read_skip_const_b(line,&sim->vectorized);
//		if (strstr(line,"Collocated")) read_skip_const_b(line,&sim->collocated);

		if (strstr(line,"Dimension")) read_skip_const_i(line,&sim->d);
//		if (strstr(line,"Method"))    read_skip_const_i(line,&sim->method);
//		if (strstr(line,"Adapt"))     read_skip_const_i(line,&sim->adapt_type);
//		if (strstr(line,"PGlobal"))   read_skip_const_i(line,&sim->p);
//		if (strstr(line,"PMax"))      read_skip_const_i(line,&sim->p_max);

//		if (strstr(line,"LevelsMax")) read_skip_const_i(line,&sim->ml_max);
	}
	fclose(ctrl_file);

	set_mesh_parameters(sim);
	set_input_path(sim);

if (0)
	set_string_associations(sim);
}

void set_Simulation_elements (struct Simulation*const sim, struct const_Intrusive_List* elements)
{
	*(struct const_Intrusive_List**)&sim->elements = elements;
}





///{ \todo Remove these functions if unused.
// Setters/Getters

void set_Simulation_flags
	(struct Simulation*const sim, const bool collocated)
{
	*(bool*)&sim->collocated = collocated;
}

void set_Simulation_parameters
	(struct Simulation*const sim, const int d, const int n_var, const int n_eq)
{
	*(int*)&sim->d     = d;
	*(int*)&sim->n_var = n_var;
	*(int*)&sim->n_eq  = n_eq;
}
///}


// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Set the mesh extension based on the mesh generator used.
static void set_mesh_extension
	(const char*const mesh_generator, ///< The mesh generator.
	 char*const mesh_extension        ///< The mesh file extension.
	);

/// \brief Assemble the mesh name with full path and extension.
static void mesh_name_assemble
	(struct Simulation*const sim,                     ///< Standard.
	 const struct Mesh_Ctrl_Data*const mesh_ctrl_data ///< The \ref Mesh_Ctrl_Data.
	);

/** \brief Set the domain type based on the `MeshDomain` input.
 *	\todo Remove the redundant `MeshCurving` variable.
 */
static void set_domain_type
	(struct Simulation*const sim,                     ///< Standard.
	 const struct Mesh_Ctrl_Data*const mesh_ctrl_data ///< The \ref Mesh_Ctrl_Data.
	);

static void set_ctrl_name_full (struct Simulation*const sim, const char*const ctrl_name)
{
	strcpy((char*)sim->ctrl_name_full,"control_files/");
	strcat((char*)sim->ctrl_name_full,ctrl_name);
	strcat((char*)sim->ctrl_name_full,".ctrl");
}

static void set_string_associations (struct Simulation*const sim)
{
// use const_cast.
	// pde_index
	if (strstr(sim->pde_name,"Advection"))
		*(int*)&sim->pde_index = PDE_ADVECTION;
	else if (strstr(sim->pde_name,"Poisson"))
		*(int*)&sim->pde_index = PDE_POISSON;
	else if (strstr(sim->pde_name,"Euler"))
		*(int*)&sim->pde_index = PDE_EULER;
	else if (strstr(sim->pde_name,"NavierStokes"))
		*(int*)&sim->pde_index = PDE_NAVIERSTOKES;
	else
		EXIT_UNSUPPORTED;
}

static void set_mesh_parameters (struct Simulation*const sim)
{
	struct Mesh_Ctrl_Data mesh_ctrl_data;

	// Read ctrl info
	FILE *ctrl_file = fopen_checked(sim->ctrl_name_full); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"mesh_generator")) read_skip_c(line,mesh_ctrl_data.mesh_generator);
		if (strstr(line,"MeshPath"))       read_skip_c(line,mesh_ctrl_data.mesh_path);
		if (strstr(line,"MeshType"))       read_skip_c(line,mesh_ctrl_data.mesh_elem_type);
		if (strstr(line,"MeshCurving"))    read_skip_c(line,mesh_ctrl_data.mesh_curving);
		if (strstr(line,"MeshDomain"))     read_skip_c(line,mesh_ctrl_data.mesh_domain);
	}
	fclose(ctrl_file);

	set_mesh_extension(mesh_ctrl_data.mesh_generator,mesh_ctrl_data.mesh_extension);
	mesh_name_assemble(sim,&mesh_ctrl_data);

	set_domain_type(sim,&mesh_ctrl_data);
}

static void set_input_path (struct Simulation*const sim)
{
	char* input_path = (char*) sim->input_path;

	strcpy(input_path,"");
	strcat_path_c(input_path,"input_files",true);
	strcat_path_c(input_path,sim->pde_name,true);
	strcat_path_c(input_path,sim->pde_spec,true);
}

// Level 1 ********************************************************************************************************** //

static void set_mesh_extension (const char*const mesh_generator, char*const mesh_extension)
{
	if (strstr(mesh_generator,"gmsh"))
		strcpy(mesh_extension,".msh");
	else
		EXIT_UNSUPPORTED;
}

static void mesh_name_assemble (struct Simulation*const sim, const struct Mesh_Ctrl_Data*const mesh_ctrl_data)
{
	char* mesh_name_full = (char*)sim->mesh_name_full;

	strcpy(mesh_name_full,"");
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_path,false);
	strcat_path_c(mesh_name_full,sim->geom_name,true);
	strcat_path_c(mesh_name_full,sim->pde_name,true);
	strcat_path_c(mesh_name_full,sim->pde_spec,true);
	strcat_path_c(mesh_name_full,sim->geom_spec,true);
	strcat_path_c(mesh_name_full,sim->geom_name,false);
	strcat_path_i(mesh_name_full,sim->d);
	strcat_path_c(mesh_name_full,"D_",false);
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_curving,false);
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_elem_type,false);
	strcat_path_i(mesh_name_full,sim->ml);
	strcat_path_c(mesh_name_full,"x",false);
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_extension,false);
}

static void set_domain_type (struct Simulation*const sim, const struct Mesh_Ctrl_Data*const mesh_ctrl_data)
{
	if (strstr(mesh_ctrl_data->mesh_domain,"Straight")) {
		const_cast_i(&sim->domain_type,DOM_STRAIGHT);
		if (!strstr(mesh_ctrl_data->mesh_curving,"Straight"))
			EXIT_UNSUPPORTED;
	} else if (strstr(mesh_ctrl_data->mesh_domain,"Curved")) {
		const_cast_i(&sim->domain_type,DOM_CURVED);
		if (!strstr(mesh_ctrl_data->mesh_curving,"Curved"))
			EXIT_UNSUPPORTED;
	} else if (strstr(mesh_ctrl_data->mesh_domain,"Mapped")) {
		const_cast_i(&sim->domain_type,DOM_MAPPED);
		if (!strstr(mesh_ctrl_data->mesh_curving,"ToBeCurved"))
			EXIT_UNSUPPORTED;
	} else {
		EXIT_UNSUPPORTED;
	}
}
