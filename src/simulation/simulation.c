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
	char mesh_generator[STRLEN_MAX], ///< The name of the file used to generate the mesh.
	     mesh_format[STRLEN_MAX],    ///< Format of the input mesh. Options: gmsh.
	     mesh_domain[STRLEN_MIN],    ///< Type of the mesh. Options: straight, curved, parametric.
	     mesh_elem_type[STRLEN_MIN], ///< Type of elements present in the mesh; used to set the mesh file name.
	     mesh_path[STRLEN_MAX],      ///< Relative path to the meshes directory.

	     mesh_extension[STRLEN_MIN]; ///< File extension (set based on \ref mesh_format.
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
	free(sim);
}

void set_simulation_core (struct Simulation*const sim, const char*const ctrl_name)
{
/// \todo Remove commented variables below.

	set_ctrl_name_full(sim,ctrl_name);
	FILE *ctrl_file = fopen_checked(sim->ctrl_name_full);

	// Read information
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"pde_name"))  read_skip_const_c(line,sim->pde_name);
		if (strstr(line,"pde_spec"))  read_skip_const_c(line,sim->pde_spec);
		if (strstr(line,"geom_name")) read_skip_const_c(line,sim->geom_name);
		if (strstr(line,"geom_spec")) read_skip_const_c(line,sim->geom_spec);

		if (strstr(line,"dimension"))        read_skip_const_i(line,&sim->d);
		if (strstr(line,"mesh_level"))       read_skip_const_i(line,&sim->ml);
		if (strstr(line,"mesh_unrealistic")) read_skip_const_b(line,&sim->mesh_unrealistic);

//		if (strstr(line,"NodeType"))  read_skip_const_c(line,sim->node_type);
//		if (strstr(line,"BasisType")) read_skip_const_c(line,sim->basis_type);

//		if (strstr(line,"Vectorized")) read_skip_const_b(line,&sim->vectorized);
//		if (strstr(line,"Collocated")) read_skip_const_b(line,&sim->collocated);

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
	strcpy((char*)sim->ctrl_name_full,"../control_files/");
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

// Level 1 ********************************************************************************************************** //

/** \brief Extract the name of the mesh generation file (excluding the path and the extension).
 *	\return See brief. */
static const char* extract_mesh_gen_name
	(const char*const mesh_generator ///< Defined in \ref Mesh_Ctrl_Data.
	);

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

	strcpy(mesh_name_full,"");
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_path,NULL);
	strcat_path_c(mesh_name_full,sim->geom_name,"/");
	strcat_path_c(mesh_name_full,sim->pde_name,"/");
	strcat_path_c(mesh_name_full,sim->pde_spec,"/");
	strcat_path_c(mesh_name_full,sim->geom_spec,"/");
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_domain,"__");

	const char*const mesh_gen_name = extract_mesh_gen_name(mesh_ctrl_data->mesh_generator); // free
	strcat_path_c(mesh_name_full,mesh_gen_name,"__");
	free((void*)mesh_gen_name);

	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_elem_type,"_ml");
	strcat_path_i(mesh_name_full,sim->ml);
	strcat_path_c(mesh_name_full,mesh_ctrl_data->mesh_extension,NULL);
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

// Level 2 ********************************************************************************************************** //

static const char* extract_mesh_gen_name (const char*const mesh_generator)
{
	int len_gen_name = 0;

	const char* beg = NULL, // Pointer to the first element of the mesh_gen_name.
	          * end = NULL; // Pointer to one past the last element of the mesh_gen_name.

	bool found_extension = false;
	int ind = 0;
	for (int i = strlen(mesh_generator)-1; i >= 0; --i) {
		const char token = mesh_generator[i];
		if (token == '/') {
			beg = &mesh_generator[i+1];
			break;
		}

		if (token == '.') {
			end = &mesh_generator[i];
			found_extension = true;
		}

		if (!found_extension)
			continue;

		++len_gen_name; // Note: One longer than the number of characters.
	}

	char*const mesh_gen_name = calloc(len_gen_name , sizeof *mesh_gen_name); // returned

	ind = 0;
	while (end != beg) {
		mesh_gen_name[ind] = *beg++;
		++ind;
	}

	return mesh_gen_name;
}
