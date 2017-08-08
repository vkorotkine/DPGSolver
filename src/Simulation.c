// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "Simulation.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#include "Macros.h"
#include "file_processing.h"

///< \brief Holds data relating to the mesh as read from the control file.
struct Mesh_Ctrl_Data {
	char mesh_path[STRLEN_MAX],      ///< Relative path to the meshes directory.
	     mesh_elem_type[STRLEN_MIN], ///< Type of elements present in the mesh; used to set the mesh file name.
	     mesh_curving[STRLEN_MAX],   /**< Type of the mesh.
	                                  *   	Options: Straight, Curved, Parametric. */
	     mesh_generator[STRLEN_MAX], /**< Type of mesh generator used.
	                                  *   	Options: gmsh */
	     mesh_extension[STRLEN_MIN]; ///< File extension (set based on \ref mesh_generator.
};

static void set_ctrl_name_full      (struct Simulation*const sim, const char*const ctrl_name);
static void set_string_associations (struct Simulation*const sim);
static void set_mesh_extension      (const char*const mesh_generator, char*const mesh_extension);
static void mesh_name_assemble
	(const struct Simulation*const sim, const struct Mesh_Ctrl_Data*const mesh_ctrl_data);

void set_mesh_name_full (const struct Simulation*const sim);


struct Simulation* constructor_Simulation ()
{
	struct Simulation* sim = malloc(sizeof *sim); // returned;

	return sim;
}

void destructor_Simulation (struct Simulation* sim)
{
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

		if (strstr(line,"MeshLevel")) read_skip_const_ui(line,&sim->ml);

//		if (strstr(line,"NodeType"))  read_skip_const_c(line,sim->node_type);
//		if (strstr(line,"BasisType")) read_skip_const_c(line,sim->basis_type);

//		if (strstr(line,"Vectorized")) read_skip_const_b(line,&sim->vectorized);
//		if (strstr(line,"Collocated")) read_skip_const_b(line,&sim->collocated);

		if (strstr(line,"Dimension")) read_skip_const_ui(line,&sim->d);
//		if (strstr(line,"Method"))    read_skip_const_ui(line,&sim->method);
//		if (strstr(line,"Adapt"))     read_skip_const_ui(line,&sim->adapt_type);
//		if (strstr(line,"PGlobal"))   read_skip_const_ui(line,&sim->p);
//		if (strstr(line,"PMax"))      read_skip_const_ui(line,&sim->p_max);

//		if (strstr(line,"LevelsMax")) read_skip_const_ui(line,&sim->ml_max);
	}
	fclose(ctrl_file);

	set_mesh_name_full(sim);

if (0)
	set_string_associations(sim);
}

void set_mesh_name_full (const struct Simulation*const sim)
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
	}
	fclose(ctrl_file);

	// Set the mesh name
	set_mesh_extension(mesh_ctrl_data.mesh_generator,mesh_ctrl_data.mesh_extension);
	mesh_name_assemble(sim,&mesh_ctrl_data);
}






// Setters/Getters

void set_Simulation_flags
	(struct Simulation*const sim, const bool collocated)
{
	*(bool*)&sim->collocated = collocated;
}

void set_Simulation_parameters
	(struct Simulation*const sim, const unsigned int d, const unsigned int n_var, const unsigned int n_eq)
{
	*(unsigned int*)&sim->d     = d;
	*(unsigned int*)&sim->n_var = n_var;
	*(unsigned int*)&sim->n_eq  = n_eq;
}

void set_Simulation_element (struct Simulation*const sim, const struct Element*const e_head)
{
	*(const struct Element**)&sim->element_head = e_head;
}

void set_Simulation_volume (struct Simulation*const sim, struct Volume* v_head)
{
	sim->volume_head = v_head;
}

void set_Simulation_face (struct Simulation*const sim, struct Face* f_head)
{
	sim->face_head = f_head;
}



// Static functions ************************************************************************************************* //

/// \brief Set full control file name (including path and file extension).
static void set_ctrl_name_full
	(struct Simulation*const sim, /// Standard.
	 const char*const ctrl_name   /// Defined in \ref set_simulation_core.
	)
{
	strcpy((char*)sim->ctrl_name_full,"control_files/");
	strcat((char*)sim->ctrl_name_full,ctrl_name);
	strcat((char*)sim->ctrl_name_full,".ctrl");
}

/** \brief Set associations between `char*` and `unsigned int` variables.
 *
 *	This is done such that if/switch conditions are simplified when these variables are used.
 *
 *	\todo This should be moved to the Solver context.
 */
static void set_string_associations (struct Simulation*const sim)
{
	// pde_index
	if (strstr(sim->pde_name,"Advection"))
		*(unsigned int*)&sim->pde_index = PDE_ADVECTION;
	else if (strstr(sim->pde_name,"Poisson"))
		*(unsigned int*)&sim->pde_index = PDE_POISSON;
	else if (strstr(sim->pde_name,"Euler"))
		*(unsigned int*)&sim->pde_index = PDE_EULER;
	else if (strstr(sim->pde_name,"NavierStokes"))
		*(unsigned int*)&sim->pde_index = PDE_NAVIERSTOKES;
	else
		EXIT_UNSUPPORTED;
}

/// \brief Set the mesh extension based on the mesh generator used.
static void set_mesh_extension (const char*const mesh_generator, char*const mesh_extension)
{
	if (strstr(mesh_generator,"gmsh"))
		strcpy(mesh_extension,".msh");
	else
		EXIT_UNSUPPORTED;
}

/// \brief Assemble the mesh name with full path and extension.
static void mesh_name_assemble (const struct Simulation*const sim, const struct Mesh_Ctrl_Data*const mesh_ctrl_data)
{
	char* mesh_name_full = (char*)sim->mesh_name_full;

	strcpy(mesh_name_full,"");
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_path,false);
	strcat_path_c(mesh_name_full,sim->geom_name,true);
	strcat_path_c(mesh_name_full,sim->pde_name,true);
	strcat_path_c(mesh_name_full,sim->pde_spec,true);
	strcat_path_c(mesh_name_full,sim->geom_spec,true);
	strcat_path_c(mesh_name_full,sim->geom_name,false);
	strcat_path_ui(mesh_name_full,sim->d);
	strcat_path_c(mesh_name_full,"D_",false);
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_curving,false);
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_elem_type,false);
	strcat_path_ui(mesh_name_full,sim->ml);
	strcat_path_c(mesh_name_full,"x",false);
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_extension,false);
}

