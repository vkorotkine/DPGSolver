// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "Simulation.h"

#include <stdlib.h>
#include <string.h>

#include "Macros.h"
#include "allocators.h"

static void set_ctrl_name_full      (struct Simulation*const sim, const char*const ctrl_name);
static void set_const_ui            (const char*const line, const unsigned int*const sim_var);
static void set_const_char_p        (const char*const line, const char*const sim_var);
static void set_const_bool          (const char*const line, const bool*const sim_var);
static void set_string_associations (struct Simulation*const sim);



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
/// \todo Potentially remove commented variables below.

	// Open control file
	set_ctrl_name_full(sim,ctrl_name);

	FILE *ctrl_file = fopen(sim->ctrl_name_full,"r");
	if (ctrl_file == NULL) {
		printf("Control file: '%s' is not present.\n",sim->ctrl_name_full);
		EXIT_UNSUPPORTED;
	}

	// Read information
	char* line = mallocator(CHAR_T,1,STRLEN_MAX); // free
	while(fscanf(ctrl_file,"%[^\n]\n",line) == 1) {
		if (strstr(line,"PDEName"))  set_const_char_p(line,sim->pde_name);
		if (strstr(line,"PDESpec"))  set_const_char_p(line,sim->pde_spec);
		if (strstr(line,"Geometry")) set_const_char_p(line,sim->geom_name);
		if (strstr(line,"GeomSpec")) set_const_char_p(line,sim->geom_spec);

		if (strstr(line,"MeshPath"))    set_const_char_p(line,sim->mesh_path);
		if (strstr(line,"MeshType"))    set_const_char_p(line,sim->mesh_elem_type);
		if (strstr(line,"MeshCurving")) set_const_char_p(line,sim->mesh_type);
		if (strstr(line,"MeshLevel"))   set_const_ui(line,&sim->ml);

//		if (strstr(line,"NodeType"))  set_const_char_p(line,sim->node_type);
//		if (strstr(line,"BasisType")) set_const_char_p(line,sim->basis_type);

		if (strstr(line,"Vectorized")) set_const_bool(line,&sim->vectorized);
		if (strstr(line,"Collocated")) set_const_bool(line,&sim->collocated);

//		if (strstr(line,"Dimension")) set_const_ui(line,&sim->d);
//		if (strstr(line,"Method"))    set_const_ui(line,&sim->method);
		if (strstr(line,"Adapt"))     set_const_ui(line,&sim->adapt_type);
		if (strstr(line,"PGlobal"))   set_const_ui(line,&sim->p);
		if (strstr(line,"PMax"))      set_const_ui(line,&sim->p_max);

//		if (strstr(line,"LevelsMax")) set_const_ui(line,&sim->ml_max);
	}
	deallocator(line,CHAR_T,1,STRLEN_MAX); // free

if (0)
	set_string_associations(sim);
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

/// \brief Set `const unsigned int` variable in \ref Simulation.
static void set_const_ui
	(const char*const line,           ///< Line from which to read data.
	 const unsigned int*const sim_var ///< Simulation variable.
	)
{
	sscanf(line,"%*s %u",(unsigned int*)sim_var);
}

/// \brief Set `const char*const` variable in \ref Simulation.
static void set_const_char_p
	(const char*const line,   ///< Line from which to read data.
	 const char*const sim_var ///< Simulation variable.
	)
{
	sscanf(line,"%*s %s",(char*)sim_var);
}

static void set_const_bool
	(const char*const line,   ///< Line from which to read data.
	 const bool*const sim_var ///< Simulation variable.
	)
{
	sscanf(line,"%*s %d",(int*)&sim_var);
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

