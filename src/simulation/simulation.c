// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "simulation.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

#include "macros.h"
#include "definitions_mesh.h"
#include "definitions_intrusive.h"

#include "element.h"
#include "mesh.h"
#include "volume.h"
#include "solver_volume.h"
#include "face.h"
#include "file_processing.h"
#include "const_cast.h"
#include "intrusive.h"

// Static function declarations ************************************************************************************* //

/** \brief Set core parameters for the simulation as specified in the control file.
 *
 *  Requires `sim` to be dynamically allocated. This allows for the definition of `const` members after the declaration
 *  which would otherwise be undefined behaviour.
 */
static void set_simulation_core
	(struct Simulation*const sim, ///< Standard.
	 const char*const ctrl_name   ///< Control file name (excluding the file extension).
	);

// Interface functions ********************************************************************************************** //

struct Simulation* constructor_Simulation (const char*const ctrl_name)
{
	struct Simulation* sim = calloc(1,sizeof *sim); // returned;

	set_simulation_core(sim,ctrl_name);
	set_Simulation_elements(sim,constructor_Elements(sim->d));

	return sim;
}

void destructor_Simulation (struct Simulation* sim)
{
// Add function pointers here when this gets bigger.

	switch (sim->elements->name) {
		case IL_ELEMENT: destructor_const_Elements(sim->elements); break;
		default:         EXIT_UNSUPPORTED;                         break;
	}


	switch (sim->volumes->name) {
		case IL_VOLUME:        destructor_Volumes(sim->volumes);        break;
		case IL_SOLVER_VOLUME: destructor_Solver_Volumes(sim->volumes); break;
		default:               EXIT_UNSUPPORTED;                        break;
	}

	switch (sim->faces->name) {
		case IL_FACE:        destructor_Faces(sim->faces);        break;
//		case IL_SOLVER_FACE: destructor_Solver_Faces(sim->faces); break;
		default:             EXIT_UNSUPPORTED;                    break;
	}

	free(sim);
}

struct Mesh_Input set_Mesh_Input (const struct Simulation*const sim)
{
	struct Mesh_Input mesh_input =
		{ .d                = sim->d,
		  .domain_type      = sim->domain_type,
		  .mesh_unrealistic = sim->mesh_unrealistic,
		  .mesh_name_full   = sim->mesh_name_full,
		  .geom_name        = sim->geom_name,
		  .geom_spec        = sim->geom_spec,
		  .input_path       = sim->input_path, };

	return mesh_input;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

///< \brief Holds data relating to the mesh as read from the control file.
struct Mesh_Ctrl_Data {
	char mesh_generator[STRLEN_MAX], ///< The name of the file used to generate the mesh.
	     mesh_format[STRLEN_MAX],    ///< Format of the input mesh. Options: gmsh.
	     mesh_domain[STRLEN_MIN],    ///< Type of the mesh. Options: straight, curved, parametric.
	     mesh_elem_type[STRLEN_MIN], ///< Type of elements present in the mesh; used to set the mesh file name.
	     mesh_path[STRLEN_MAX],      ///< Relative path to the meshes directory.

	     mesh_extension[STRLEN_MIN]; ///< File extension (set based on \ref mesh_format.
};

/** \brief Set full control file name (including path and file extension).
 *
 *	If "TEST" is not included as part of the name (default option):
 *	- it is assumed that `ctrl_name` includes the full name and path from CMAKE_PROJECT_DIR/control_files;
 *	- `ctrl_name_full` -> ../control_files/`ctrl_name`.
 *	otherwise:
 *	- it is assumed that `ctrl_name` includes the full name and path from CMAKE_PROJECT_DIR/testing/control_files;
 *	- `ctrl_name_full` -> ../testing/control_files/`ctrl_name`.
 */
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

static void set_simulation_core (struct Simulation*const sim, const char*const ctrl_name)
{
	set_ctrl_name_full(sim,ctrl_name);
	FILE *ctrl_file = fopen_checked(sim->ctrl_name_full);

	// Read information
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"pde_name"))  read_skip_const_c_1(line,sim->pde_name);
		if (strstr(line,"pde_spec"))  read_skip_const_c_1(line,sim->pde_spec);
		if (strstr(line,"geom_name")) read_skip_const_c_1(line,sim->geom_name);
		if (strstr(line,"geom_spec")) read_skip_const_c_1(line,sim->geom_spec);

		if (strstr(line,"dimension"))        read_skip_const_i_1(line,1,&sim->d,1);
		if (strstr(line,"mesh_level"))       read_skip_const_i_1(line,1,sim->ml,2);
		if (strstr(line,"mesh_unrealistic")) read_skip_const_b(line,&sim->mesh_unrealistic);

		if (strstr(line,"basis_geom")) read_skip_const_c_1(line,sim->basis_geom);
		if (strstr(line,"basis_sol"))  read_skip_const_c_1(line,sim->basis_sol);
		if (strstr(line,"geom_rep"))   read_skip_const_c_1(line,sim->geom_rep);

		if (strstr(line,"p_s_v"))    read_skip_const_i_1(line,1,sim->p_s_v,2);
		if (strstr(line,"p_s_f"))    read_skip_const_i_1(line,1,sim->p_s_f,2);
		if (strstr(line,"p_sg_v"))   read_skip_const_i_1(line,1,sim->p_sg_v,2);
		if (strstr(line,"p_sg_f"))   read_skip_const_i_1(line,1,sim->p_sg_f,2);
		if (strstr(line,"p_cub_x"))  read_skip_const_i_1(line,1,&sim->p_c_x,1);
		if (strstr(line,"p_cub_p"))  read_skip_const_i_1(line,1,&sim->p_c_p,1);
		if (strstr(line,"p_test_p")) read_skip_const_i_1(line,1,&sim->p_t_p,1);
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

// Level 1 ********************************************************************************************************** //

/// \brief Set the mesh extension based on the mesh generator used.
static void set_mesh_extension
	(const char*const mesh_generator, ///< The mesh generator.
	 char*const mesh_extension        ///< The mesh file extension.
	);

/** \brief Assemble the mesh name with full path and extension.
 *	If the standard mesh path (../meshes/) is used, assemble the full mesh name based on data included in the control
 *	file. Otherwise, it is assumed that the full mesh name (including the path) was provided in the control file.
 */
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
	strcpy((char*)sim->ctrl_name_full,"../");
	if (strstr(ctrl_name,"TEST"))
		strcat((char*)sim->ctrl_name_full,"testing/");
	strcat((char*)sim->ctrl_name_full,"control_files/");
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
		if (strstr(line,"mesh_format"))    read_skip_c(line,mesh_ctrl_data.mesh_format);
		if (strstr(line,"mesh_domain"))    read_skip_c(line,mesh_ctrl_data.mesh_domain);
		if (strstr(line,"mesh_type"))      read_skip_c(line,mesh_ctrl_data.mesh_elem_type);
		if (strstr(line,"mesh_path"))      read_skip_c(line,mesh_ctrl_data.mesh_path);
	}
	fclose(ctrl_file);

	set_mesh_extension(mesh_ctrl_data.mesh_format,mesh_ctrl_data.mesh_extension);
	mesh_name_assemble(sim,&mesh_ctrl_data);

	set_domain_type(sim,&mesh_ctrl_data);
}

static void set_input_path (struct Simulation*const sim)
{
	char* input_path = (char*) sim->input_path;

	strcpy(input_path,"");
	strcat_path_c(input_path,"input_files","/");
	strcat_path_c(input_path,sim->pde_name,"/");
	strcat_path_c(input_path,sim->pde_spec,"/");
}

// Level 2 ********************************************************************************************************** //

static void set_mesh_extension (const char*const mesh_format, char*const mesh_extension)
{
	if (strstr(mesh_format,"gmsh"))
		strcpy(mesh_extension,".msh");
	else
		EXIT_UNSUPPORTED;
}

static void mesh_name_assemble (struct Simulation*const sim, const struct Mesh_Ctrl_Data*const mesh_ctrl_data)
{
	char* mesh_name_full = (char*)sim->mesh_name_full;

	strcpy(mesh_name_full,mesh_ctrl_data->mesh_path);
	if (strstr(mesh_ctrl_data->mesh_path,"../meshes/")) {
		strcat_path_c(mesh_name_full,sim->geom_name,"/");
		strcat_path_c(mesh_name_full,sim->pde_name,"/");
		strcat_path_c(mesh_name_full,sim->pde_spec,"/");
		strcat_path_c(mesh_name_full,sim->geom_spec,"/");
		strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_domain,"__");

		const char*const mesh_gen_name = extract_name(mesh_ctrl_data->mesh_generator,true); // free
		strcat_path_c(mesh_name_full,mesh_gen_name,"__");
		free((void*)mesh_gen_name);

		strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_elem_type,"_ml");
		strcat_path_i(mesh_name_full,sim->ml[0]);
		strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_extension,NULL);
	}
}

static void set_domain_type (struct Simulation*const sim, const struct Mesh_Ctrl_Data*const mesh_ctrl_data)
{
	if (strstr(mesh_ctrl_data->mesh_domain,"straight"))
		const_cast_i(&sim->domain_type,DOM_STRAIGHT);
	else if (strstr(mesh_ctrl_data->mesh_domain,"curved"))
		const_cast_i(&sim->domain_type,DOM_CURVED);
	else if (strstr(mesh_ctrl_data->mesh_domain,"parametric"))
		const_cast_i(&sim->domain_type,DOM_PARAMETRIC);
	else
		EXIT_UNSUPPORTED;
}
