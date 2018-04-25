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

#include "simulation.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include "mpi.h"
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_core.h"
#include "definitions_mesh.h"
#include "definitions_intrusive.h"
#include "definitions_test_case.h"

#include "element.h"
#include "volume.h"
#include "face.h"

#include "const_cast.h"
#include "file_processing.h"
#include "geometry.h"
#include "intrusive.h"
#include "mesh.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

///\{ \name Invalid value for a polynomial order (which would be unlikely to be chosen inadvertently).
#define P_INVALID -99999
///\}

/// \brief Set the MPI related members of \ref Simulation.
static void set_simulation_mpi
	(struct Simulation*const sim ///< Standard.
	);

/** \brief Set core parameters for the simulation as specified in the control file.
 *
 *  Requires `sim` to be dynamically allocated. This allows for the definition of `const` members after the declaration
 *  which would otherwise be undefined behaviour.
 */
static void set_simulation_core
	(struct Simulation*const sim, ///< Standard.
	 const char*const ctrl_name   ///< Control file name (excluding the file extension).
	);

/// \brief Set additional parameters requiring that \ref set_up_fopen_input has been called.
static void set_simulation_additional
	(struct Simulation*const sim ///< Standard.
	);

/// \brief Set \ref Simulation parameters to invalid values so that it can be recognized if they were not read.
static void set_simulation_invalid
	(struct Simulation*const sim ///< Standard.
	);

/// \brief Set uninitialized \ref Simulation parameters to default values.
static void set_simulation_default
	(struct Simulation*const sim ///< Standard.
	);

/** \brief Check that necessary \ref Simulation parameters were set.
 *
 *  Mesh generation related parameters need not be checked as mesh generation would have failed if they were not
 *  present.
 */
static void check_necessary_simulation_parameters
	(struct Simulation*const sim ///< Standard.
	);

/** \brief Check whether h/p-adaption should be enabled based on the input array.
 *  \return `true` if yes. */
static bool is_adaptive
	(const int var[2] ///< The array of minimal and maximal orders/mesh levels.
	);

/** \brief Get the path to the input files.
 *  \return A statically allocated `char[]` holding the input path. */
static const char* get_input_path
	(struct Simulation*const sim ///< Standard.
	);

/** \brief Return a stack allocated `char*` holding the name (including full path) of the restart file.
 *  \return See brief. */
static const char* get_restart_name ( );

// Interface functions ********************************************************************************************** //

struct Simulation* constructor_Simulation__no_mesh (const char*const ctrl_name)
{
	struct Simulation* sim = calloc(1,sizeof *sim); // returned;

	set_simulation_invalid(sim);
	set_simulation_mpi(sim);
	set_simulation_core(sim,ctrl_name);
	set_up_fopen_input(sim->ctrl_name_full,get_input_path(sim));
	set_simulation_additional(sim);

	check_necessary_simulation_parameters(sim);
	set_simulation_default(sim);

	sim->elements = NULL;
	sim->volumes  = NULL;
	sim->faces    = NULL;

	sim->test_case_rc = constructor_Test_Case_rc_real(sim);

	return sim;
}

struct Simulation* constructor_Simulation (const char*const ctrl_name)
{
	struct Simulation* sim = constructor_Simulation__no_mesh(ctrl_name); // returned

	set_Simulation_elements(sim,constructor_Elements(DIM)); // destructed

	struct Mesh_Input mesh_input = set_Mesh_Input(sim);
	struct Mesh* mesh = constructor_Mesh(&mesh_input,sim->elements); // destructed
	remove_absent_Elements(sim->elements);

	sim->volumes = constructor_Volumes(sim,mesh); // destructed
	sim->faces   = constructor_Faces(sim,mesh);   // destructed

	destructor_Mesh(mesh);

	return sim;
}

struct Simulation* constructor_Simulation_restart (const struct Simulation*const sim_main)
{
	struct Simulation* sim = constructor_Simulation__no_mesh(sim_main->ctrl_name); // returned

	set_Simulation_elements(sim,constructor_Elements(DIM)); // destructed

	struct Mesh_Input mesh_input = set_Mesh_Input(sim);
	mesh_input.mesh_name_full = get_restart_name();
	struct Mesh* mesh = constructor_Mesh(&mesh_input,sim->elements); // destructed
	remove_absent_Elements(sim->elements);

	sim->volumes = constructor_Volumes(sim,mesh); // destructed
	sim->faces   = constructor_Faces(sim,mesh);   // destructed

	destructor_Mesh(mesh);

	return sim;
}

void destructor_Simulation (struct Simulation* sim)
{
	if (sim->elements) {
		assert(sim->elements->name == IL_ELEMENT);
		destructor_const_Elements(sim->elements);
	}

	if (sim->volumes) {
		assert(sim->volumes->name == IL_VOLUME);
		destructor_Volumes(sim->volumes);
	}

	if (sim->faces) {
		assert(sim->faces->name == IL_FACE);
		destructor_Faces(sim->faces);
	}

	destructor_Test_Case_rc_real(sim->test_case_rc);

	free(sim);
}

const char* set_ctrl_name_full (const char*const ctrl_name)
{
	static char ctrl_name_full[STRLEN_MAX] = { 0, };

	strcpy(ctrl_name_full,PROJECT_INPUT_DIR);
	if (strstr(ctrl_name,"TEST"))
		strcat(ctrl_name_full,"testing/");
	if (strstr(ctrl_name,"RESULTS"))
		strcat(ctrl_name_full,"results/");
	strcat(ctrl_name_full,"control_files/");
	strcat(ctrl_name_full,ctrl_name);
	strcat(ctrl_name_full,".ctrl");

	return ctrl_name_full;
}

struct Mesh_Input set_Mesh_Input (const struct Simulation*const sim)
{
	struct Mesh_Input mesh_input =
		{ .d                = DIM,
		  .domain_type      = sim->domain_type,
		  .mesh_unrealistic = sim->mesh_unrealistic,
		  .mesh_name_full   = sim->mesh_name_full,
		  .geom_name        = sim->geom_name,
		  .geom_spec        = sim->geom_spec,
		};

	return mesh_input;
}

int compute_adapt_type (const int p_ref[2], const int ml[2])
{
	const bool p_adapt = is_adaptive(p_ref),
	           h_adapt = is_adaptive(ml);

	if (!p_adapt && !h_adapt)
		return ADAPT_0;
	else if (p_adapt && !h_adapt)
		return ADAPT_P;
	else if (!p_adapt && h_adapt)
		return ADAPT_H;
	else if (p_adapt && h_adapt)
		return ADAPT_HP;

	EXIT_UNSUPPORTED;
	return -1;
}

void set_ml_p_curr (const int ml, const int p, struct Simulation* sim)
{
	const_cast_i_n(sim->ml_p_curr,(int[]){ml,p},2);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Holds data relating to the mesh as read from the control file.
struct Mesh_Ctrl_Data {
	char mesh_generator[STRLEN_MAX], ///< The name of the file used to generate the mesh.
	     mesh_format[STRLEN_MAX],    ///< Format of the input mesh. Options: gmsh.
	     mesh_domain[STRLEN_MIN],    ///< Type of the mesh. Options: straight, blended, parametric.
	     mesh_elem_type[STRLEN_MIN], ///< Type of elements present in the mesh; used to set the mesh file name.
	     mesh_path[STRLEN_MAX],      ///< Relative path to the meshes directory.

	     mesh_extension[STRLEN_MIN]; ///< File extension (set based on \ref mesh_format.
};

/// \brief Set the mesh parameters.
static void set_mesh_parameters
	(struct Simulation*const sim ///< Standard.
	);

/// \brief Set the orders of the various solution components based on the input data.
static void set_orders
	(struct Simulation*const sim ///< Standard.
	);

static void set_simulation_invalid (struct Simulation*const sim)
{
	const_cast_i(&sim->mpi_size,-1);
	const_cast_i(&sim->mpi_rank,-1);

	for (int i = 0; i < N_ST_STD; ++i) {
		const_cast_c(sim->nodes_interp[i],0);
		const_cast_c(sim->geom_blending[i],0);
	}

	const_cast_c(sim->test_case_extension,0);
	const_cast_c(sim->solution_extension,0);
	const_cast_c(sim->basis_geom,0);
	const_cast_c(sim->basis_sol,0);
	const_cast_c(sim->geom_rep,0);

	const_cast_i_n(sim->p_ref,(int[]){P_INVALID,P_INVALID},2);
	const_cast_i_n(sim->p_ig,(int[]){P_INVALID,P_INVALID},2);
	const_cast_i_n(sim->p_s_v,(int[]){P_INVALID,P_INVALID},2);
	const_cast_i_n(sim->p_s_f,(int[]){P_INVALID,P_INVALID},2);
	const_cast_i_n(sim->p_sg_v,(int[]){P_INVALID,P_INVALID},2);
	const_cast_i_n(sim->p_sg_f,(int[]){P_INVALID,P_INVALID},2);

	const_cast_i(&sim->p_s_v_p,P_INVALID);
	const_cast_i(&sim->p_s_f_p,P_INVALID);
	const_cast_i(&sim->p_sg_v_p,P_INVALID);
	const_cast_i(&sim->p_sg_f_p,P_INVALID);
	const_cast_i_n(sim->p_c_x,(int[]){P_INVALID,P_INVALID},2);
	const_cast_i_n(sim->p_c_p,(int[]){P_INVALID,P_INVALID},2);
	const_cast_i_n(sim->p_t_p,(int[]){P_INVALID,P_INVALID},2);

	set_ml_p_curr(-1,-1,sim);

	const_cast_i(&sim->method,-1);

	const_cast_b(&sim->collocated,false);
}

static void set_simulation_mpi (struct Simulation*const sim)
{
	int mpi_size = -1,
	    mpi_rank = -1;
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	const_cast_i(&sim->mpi_size,mpi_size);
	const_cast_i(&sim->mpi_rank,mpi_rank);
}

static void set_simulation_core (struct Simulation*const sim, const char*const ctrl_name)
{
	const_cast_c1(&sim->ctrl_name,ctrl_name);
	sim->ctrl_name_full = set_ctrl_name_full(ctrl_name);
	FILE *ctrl_file = fopen_checked(sim->ctrl_name_full);

	int d = -1;

	// Read information
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"pde_name"))  read_skip_const_c_1(line,sim->pde_name);
		if (strstr(line,"pde_spec"))  read_skip_const_c_1(line,sim->pde_spec);
		if (strstr(line,"geom_name")) read_skip_const_c_1(line,sim->geom_name);
		if (strstr(line,"geom_spec")) read_skip_const_c_1(line,sim->geom_spec);

		if (strstr(line,"dimension"))        read_skip_const_i_1(line,1,&d,1);
		if (strstr(line,"mesh_level"))       read_skip_const_i_1(line,1,sim->ml,2);
		if (strstr(line,"mesh_unrealistic")) read_skip_const_b(line,&sim->mesh_unrealistic);

		read_skip_name_const_c_1("test_case_extension",line,sim->test_case_extension);
		read_skip_name_const_c_1("solution_extension", line,sim->solution_extension);

		if (strstr(line,"interp_tp"))  read_skip_const_c_1(line,sim->nodes_interp[0]);
		if (strstr(line,"interp_si"))  read_skip_const_c_1(line,sim->nodes_interp[1]);
		if (strstr(line,"interp_pyr")) read_skip_const_c_1(line,sim->nodes_interp[2]);

		if (strstr(line,"basis_geom")) read_skip_const_c_1(line,sim->basis_geom);
		if (strstr(line,"basis_sol"))  read_skip_const_c_1(line,sim->basis_sol);
		if (strstr(line,"geom_rep"))   read_skip_const_c_1(line,sim->geom_rep);

		if (strstr(line,"geom_blending_tp"))  read_skip_const_c_1(line,sim->geom_blending[0]);
		if (strstr(line,"geom_blending_si"))  read_skip_const_c_1(line,sim->geom_blending[1]);
		if (strstr(line,"geom_blending_pyr")) read_skip_const_c_1(line,sim->geom_blending[2]);

		if (strstr(line,"p_ref"))    read_skip_const_i_1(line,1,sim->p_ref,2);
		if (strstr(line,"p_s_v_p"))  read_skip_const_i(line,&sim->p_s_v_p);
		if (strstr(line,"p_s_f_p"))  read_skip_const_i(line,&sim->p_s_f_p);
		if (strstr(line,"p_sg_v_p")) read_skip_const_i(line,&sim->p_sg_v_p);
		if (strstr(line,"p_sg_f_p")) read_skip_const_i(line,&sim->p_sg_f_p);
		if (strstr(line,"p_cub_x"))  read_skip_const_i_1(line,1,sim->p_c_x,2);
		if (strstr(line,"p_cub_p"))  read_skip_const_i_1(line,1,sim->p_c_p,2);
		if (strstr(line,"p_test_p")) read_skip_const_i_1(line,1,sim->p_t_p,2);

		if (strstr(line,"fe_method")) read_skip_const_i_1(line,1,&sim->method,1);

		if (strstr(line,"collocated")) read_skip_const_b(line,&sim->collocated);
	}
	fclose(ctrl_file);

	if (d != DIM)
		EXIT_ERROR("Using executable of incorrect dimension (%d (mesh) != %d (exec)).\n",d,DIM);

	set_mesh_parameters(sim);
	set_orders(sim);
}

static void set_simulation_additional (struct Simulation*const sim)
{
	if (is_internal_geom_straight()) {
		const_cast_i(&sim->p_ig[0],GSL_MIN(sim->p_ref[0],0));
		const_cast_i(&sim->p_ig[1],GSL_MAX(sim->p_ref[1],0));

		// In the current implementation, non-conforming meshes will have gaps if this condition is not satisfied.
		if (!((strcmp(sim->geom_rep,"isoparametric")          == 0) ||
	            (strcmp(sim->geom_rep,"superparametric")        == 0) ||
	            (strcmp(sim->geom_rep,"superparametric_p_le_1") == 0)))
			EXIT_ERROR("Cannot use internally straight geometry for geom_rep: %s",sim->geom_rep);
	} else {
		const_cast_i_n(sim->p_ig,sim->p_ref,2);
	}
}

static void check_necessary_simulation_parameters (struct Simulation*const sim)
{
	assert(sim->mpi_size >  0);
	assert(sim->mpi_rank >= 0);

	assert((sim->nodes_interp[0][0] != 0) ||
	       (sim->nodes_interp[1][0] != 0) ||
	       (sim->nodes_interp[2][0] != 0));

	assert((strcmp(sim->basis_geom,"lagrange") == 0) ||
	       (strcmp(sim->basis_geom,"bezier")   == 0) ||
	       (strcmp(sim->basis_geom,"nurbs")    == 0));
	assert((strcmp(sim->basis_sol,"orthonormal") == 0) ||
	       (strcmp(sim->basis_sol,"lagrange")    == 0) ||
	       (strcmp(sim->basis_sol,"bezier")      == 0));
	assert((strcmp(sim->geom_rep,"isoparametric")          == 0) ||
	       (strcmp(sim->geom_rep,"superparametric")        == 0) ||
	       (strcmp(sim->geom_rep,"superparametric2")       == 0) ||
	       (strcmp(sim->geom_rep,"superparametric_p_le_1") == 0) ||
	       (strstr(sim->geom_rep,"fixed")                  != 0));

	assert(sim->p_s_v[0] != P_INVALID);
	assert(sim->p_s_v[1] != P_INVALID);
	assert(sim->p_s_v[1] >= sim->p_s_v[0]);

	assert((sim->method == METHOD_DPG) || (sim->p_t_p[0] == P_INVALID && sim->p_t_p[1] == P_INVALID));

	assert((sim->method == METHOD_DG)   || (sim->method == METHOD_HDG) ||
	       (sim->method == METHOD_HDPG) || (sim->method == METHOD_DPG));
}

static void set_simulation_default (struct Simulation*const sim)
{
	for (int i = 0; i < 2; ++i) {
		if (sim->p_s_f[i] == P_INVALID)
			const_cast_i(&sim->p_s_f[i],sim->p_s_v[i]);

		if (sim->p_sg_v[i] == P_INVALID)
			const_cast_i(&sim->p_sg_v[i],sim->p_s_v[i]);

		if (sim->p_sg_f[i] == P_INVALID)
			const_cast_i(&sim->p_sg_f[i],sim->p_s_v[i]);

		if (sim->p_c_x[i] == P_INVALID)
			const_cast_i(&sim->p_c_x[i],2);
		else
			assert(sim->p_c_x[i] >= 2);

		if (sim->p_c_p[i] == P_INVALID)
			const_cast_i(&sim->p_c_p[i],0);
		else
			assert(sim->p_c_p[i] >= 0);

		if (sim->p_t_p[i] == P_INVALID)
			const_cast_i(&sim->p_t_p[i],0);
		else
			assert(sim->p_t_p[i] >= 0);
	}

	const_cast_i(&sim->n_hp,3);
	const_cast_i(&sim->adapt_type,compute_adapt_type(sim->p_ref,sim->ml));
}

void set_Simulation_elements (struct Simulation*const sim, struct const_Intrusive_List* elements)
{
	*(struct const_Intrusive_List**)&sim->elements = elements;
}

static bool is_adaptive (const int var[2])
{
	if (var[0] != var[1])
		return true;
	return false;
}

static const char* get_input_path (struct Simulation*const sim)
{
	static char input_path[4*STRLEN_MAX];
	snprintf(input_path,sizeof(input_path),"%s%s%s%s%s%s",
	        PROJECT_INPUT_DIR,"input_files/",sim->pde_name,"/",sim->pde_spec,"/");
	return input_path;
}

static const char* get_restart_name ( )
{
	static char restart_path[STRLEN_MAX];
	static bool needs_input = true;

	if (needs_input) {
		needs_input = false;
		FILE* input_file = fopen_input('c',NULL,NULL); // closed
		char line[STRLEN_MAX];
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"restart_path"))
				read_skip_c(line,restart_path);
		}
		fclose(input_file);
	}

	return restart_path;
}

// Level 1 ********************************************************************************************************** //

/// \brief Set the mesh extension based on the mesh generator used.
static void set_mesh_extension
	(const char*const mesh_generator, ///< The mesh generator.
	 char*const mesh_extension        ///< The mesh file extension.
	);

/** \brief Assemble the mesh name with full path and extension.
 *  If the standard mesh path (../meshes/) is used, assemble the full mesh name based on data included in the control
 *  file. Otherwise, it is assumed that the full mesh name (including the path) was provided in the control file.
 */
static void mesh_name_assemble
	(struct Simulation*const sim,                     ///< Standard.
	 const struct Mesh_Ctrl_Data*const mesh_ctrl_data ///< The \ref Mesh_Ctrl_Data.
	);

/** \brief Set the domain type based on the `MeshDomain` input.
 */
static void set_domain_type
	(struct Simulation*const sim,                     ///< Standard.
	 const struct Mesh_Ctrl_Data*const mesh_ctrl_data ///< The \ref Mesh_Ctrl_Data.
	);

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

static void set_orders (struct Simulation*const sim)
{
	if (sim->p_s_v_p == P_INVALID)
		const_cast_i(&sim->p_s_v_p,0);

	if (sim->p_s_f_p == P_INVALID)
		const_cast_i(&sim->p_s_f_p,0);

	if (sim->p_sg_v_p == P_INVALID)
		const_cast_i(&sim->p_sg_v_p,0);

	if (sim->p_sg_f_p == P_INVALID)
		const_cast_i(&sim->p_sg_f_p,0);

	for (int i = 0; i < 2; ++i) {
		int p = sim->p_ref[i]+sim->p_s_v_p;
		const_cast_i(&sim->p_s_v[i],p);

		p = sim->p_s_f_p + (sim->p_s_f_p == P_INVALID ? 0 : sim->p_ref[i]);
		const_cast_i(&sim->p_s_f[i],p);

		p = sim->p_sg_v_p + (sim->p_sg_v_p == P_INVALID ? 0 : sim->p_ref[i]);
		const_cast_i(&sim->p_sg_v[i],p);

		p = sim->p_sg_f_p + (sim->p_sg_f_p == P_INVALID ? 0 : sim->p_ref[i]);
		const_cast_i(&sim->p_sg_f[i],p);
	}
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

	int index = sprintf(mesh_name_full,"%s",mesh_ctrl_data->mesh_path);
	if (strstr(mesh_ctrl_data->mesh_path,"../meshes/")) {
		index += sprintf(mesh_name_full+index,"%s%s%s%s%s%s",sim->geom_name,"/",sim->pde_name,"/",sim->pde_spec,"/");
		if (strcmp(sim->geom_spec,"NONE") != 0)
			index += sprintf(mesh_name_full+index,"%s%s",sim->geom_spec,"/");
		index += sprintf(mesh_name_full+index,"%s%s",mesh_ctrl_data->mesh_domain,"__");

		const char*const mesh_gen_name = extract_name(mesh_ctrl_data->mesh_generator,true);
		index += sprintf(mesh_name_full+index,"%s%s",mesh_gen_name,"__");

		index += sprintf(mesh_name_full+index,"%s%s%d",mesh_ctrl_data->mesh_elem_type,"_ml",sim->ml[0]);
		index += sprintf(mesh_name_full+index,"%s",mesh_ctrl_data->mesh_extension);
	}
}

static void set_domain_type (struct Simulation*const sim, const struct Mesh_Ctrl_Data*const mesh_ctrl_data)
{
	if (strstr(mesh_ctrl_data->mesh_domain,"straight"))
		const_cast_i(&sim->domain_type,DOM_STRAIGHT);
	else if (strstr(mesh_ctrl_data->mesh_domain,"blended"))
		const_cast_i(&sim->domain_type,DOM_BLENDED);
	else if (strstr(mesh_ctrl_data->mesh_domain,"parametric"))
		const_cast_i(&sim->domain_type,DOM_PARAMETRIC);
	else
		EXIT_UNSUPPORTED;
}
