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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_dpg.h"
#include "definitions_test_case.h"


#include "def_templates_test_case.h"

#include "def_templates_flux.h"
#include "def_templates_solution_advection.h"
#include "def_templates_solution_euler.h"

// Static function declarations ************************************************************************************* //

/// \brief Set associations between `char*` and `int` variables.
static void set_string_associations
	(struct Test_Case_T* test_case,      ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Set pde related parameters.
static void set_pde_related
	(struct Test_Case_T* test_case,      ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Set the function pointer members of \ref Test_Case_T.
static void set_function_pointers
	(struct Test_Case_T* test_case,      ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Read members of \ref Test_Case_T from input file.
static void read_test_case_parameters
	(struct Test_Case_T* test_case,      ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

struct Test_Case_T* constructor_Test_Case_T (const struct Simulation* sim)
{
	struct Test_Case_T* test_case = calloc(1,sizeof *test_case); // returned

	set_string_associations(test_case,sim);
	set_pde_related(test_case,sim);

	read_test_case_parameters(test_case,sim);
	set_function_pointers(test_case,sim);

	test_case->solver_method_curr = 0;

	return test_case;
}

void destructor_Test_Case_T (const struct Test_Case_T* test_case)
{
	free((void*)test_case);
}

void increment_pointers_T (const int n_ptr, const Type**const ptrs)
{
	for (int i = 0; i < n_ptr; ++i)
		++ptrs[i];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// Container for input strings which are to be subsequently converted to integer parameters.
struct Test_Case_String_Inputs {
	char num_flux_1st[STRLEN_MIN]; ///< The name of the 1st order numerical flux scheme to be used.
	char num_flux_2nd[STRLEN_MIN]; ///< The name of the 2nd order numerical flux scheme to be used.

	char test_norm[STRLEN_MIN]; ///< The name of the norm to use for the optimal test function computation.
};

/** \brief Return a statically allocated \ref Test_Case_String_Inputs container which zero-initialized members.
 *  \return See brief. */
static struct Test_Case_String_Inputs set_Test_Case_String_Inputs
	();

/// \brief Set the string association relating to the \ref Test_Case_T input parameters.
static void set_string_associations_test_case
	(struct Test_Case_T* test_case,               ///< \ref Test_Case_T.
	 const struct Test_Case_String_Inputs* tcsi ///< \ref Test_Case_String_Inputs.
	);

/// \brief Correct invalid \ref Test_Case_T parameters if present.
static void correct_invalid_test_case_parameters
	(struct Test_Case_T* test_case, ///< \ref Test_Case_T.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Return a statically allocated array of \ref Flux_Input_T::compute_member flags depending on the pde and
 *         solver method.
 *  \return See brief. */
static const bool* get_compute_member_Flux_Input
	(const char type_ei,                 ///< 'e'xplicit/'i'mplicit type. Options: 'e', 'i'.
	 struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

/** \brief Return a statically allocated array of \ref Boundary_Value_Input_T::compute_member flags depending on the pde
 *         and solver method.
 *  \return See brief. */
static const bool* get_compute_member_Boundary_Value_Input
	(const char type_ei,                 ///< 'e'xplicit/'i'mplicit type. Options: 'e', 'i'.
	 struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

static void set_string_associations (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	// pde_index
	if (strstr(sim->pde_name,"advection"))
		const_cast_i(&test_case->pde_index,PDE_ADVECTION);
	else if (strstr(sim->pde_name,"poisson"))
		const_cast_i(&test_case->pde_index,PDE_POISSON);
	else if (strstr(sim->pde_name,"euler"))
		const_cast_i(&test_case->pde_index,PDE_EULER);
	else if (strstr(sim->pde_name,"navier_stokes"))
		const_cast_i(&test_case->pde_index,PDE_NAVIER_STOKES);
	else
		EXIT_ERROR("Unsupported: %s\n",sim->pde_name);
}

static void set_pde_related (struct Test_Case_T* test_case, const struct Simulation* sim)
{
	switch (test_case->pde_index) {
	case PDE_ADVECTION:
		const_cast_b(&test_case->is_linear,true);
		const_cast_i(&test_case->n_var,1);
		const_cast_i(&test_case->n_eq,1);
		const_cast_b(&test_case->has_1st_order,true);
		const_cast_b(&test_case->has_2nd_order,false);
		break;
	case PDE_POISSON:
		const_cast_b(&test_case->is_linear,true);
		const_cast_i(&test_case->n_var,1);
		const_cast_i(&test_case->n_eq,1);
		const_cast_b(&test_case->has_1st_order,false);
		const_cast_b(&test_case->has_2nd_order,true);
		break;
	case PDE_EULER:
		const_cast_b(&test_case->is_linear,false);
		const_cast_i(&test_case->n_var,DIM+2);
		const_cast_i(&test_case->n_eq,DIM+2);
		const_cast_b(&test_case->has_1st_order,true);
		const_cast_b(&test_case->has_2nd_order,false);
		break;
	case PDE_NAVIER_STOKES:
		const_cast_b(&test_case->is_linear,false);
		const_cast_i(&test_case->n_var,DIM+2);
		const_cast_i(&test_case->n_eq,DIM+2);
		const_cast_b(&test_case->has_1st_order,true);
		const_cast_b(&test_case->has_2nd_order,true);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}

	const bool* flux_comp_mem_e = get_compute_member_Flux_Input('e',test_case,sim),
	          * flux_comp_mem_i = get_compute_member_Flux_Input('i',test_case,sim);
	for (int i = 0; i < MAX_FLUX_OUT; ++i) {
		const_cast_b(&test_case->flux_comp_mem_e[i],flux_comp_mem_e[i]);
		const_cast_b(&test_case->flux_comp_mem_i[i],flux_comp_mem_i[i]);
	}

	const bool* boundary_value_comp_mem_e = get_compute_member_Boundary_Value_Input('e',test_case,sim),
	          * boundary_value_comp_mem_i = get_compute_member_Boundary_Value_Input('i',test_case,sim);
	for (int i = 0; i < MAX_BV_OUT; ++i) {
		const_cast_b(&test_case->boundary_value_comp_mem_e[i],boundary_value_comp_mem_e[i]);
		const_cast_b(&test_case->boundary_value_comp_mem_i[i],boundary_value_comp_mem_i[i]);
	}
}

static void set_function_pointers (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	switch (test_case->pde_index) {
		case PDE_ADVECTION:     set_function_pointers_solution_advection_T(test_case,sim);     break;
//		case PDE_POISSON:       set_function_pointers_solution_poisson(test_case,sim);       break;
		case PDE_EULER:         set_function_pointers_solution_euler_T(test_case,sim);         break;
//		case PDE_NAVIER_STOKES: set_function_pointers_solution_navier_stokes(test_case,sim); break;
		default: EXIT_ERROR("Unsupported: %d\n",test_case->pde_index); break;
	}
}

static void read_test_case_parameters (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input(sim->input_path,'t'); // closed

	struct Test_Case_String_Inputs tcsi = set_Test_Case_String_Inputs();

	int count_found = 0,
	    count_tmp = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"solver_proc")) {
			++count_found;
			read_skip_const_i(line,&test_case->solver_proc);
		}
		if (strstr(line,"solver_type_e")) read_skip_const_i(line,&test_case->solver_type_e);
		if (strstr(line,"solver_type_i")) read_skip_const_i(line,&test_case->solver_type_i);

		if (strstr(line,"num_flux_1st")) read_skip_c_1(line,tcsi.num_flux_1st);
		if (strstr(line,"num_flux_2nd")) read_skip_c_1(line,tcsi.num_flux_2nd);

		if (strstr(line,"test_norm")) read_skip_c_1(line,tcsi.test_norm);

		if (strstr(line,"time_final")) read_skip_const_d(line,&test_case->time_final,1,false);
		if (strstr(line,"time_step"))  read_skip_const_d(line,&test_case->dt,1,false);

		if (strstr(line,"use_schur_complement")) read_skip_const_b(line,&test_case->use_schur_complement);

		if (strstr(line,"display_progress"))    read_skip_const_b(line,&test_case->display_progress);
		if (strstr(line,"conv_order_discount")) read_skip_const_d(line,&test_case->conv_order_discount,1,false);

		read_skip_string_count_const_d("exit_tol_e",  &count_tmp,line,&test_case->exit_tol_e);
		read_skip_string_count_const_d("exit_ratio_e",&count_tmp,line,&test_case->exit_ratio_e);
		read_skip_string_count_const_d("exit_tol_i",  &count_tmp,line,&test_case->exit_tol_i);
		read_skip_string_count_const_d("exit_ratio_i",&count_tmp,line,&test_case->exit_ratio_i);
	}
	fclose(input_file);

	correct_invalid_test_case_parameters(test_case,sim);
	set_string_associations_test_case(test_case,&tcsi);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

// Level 1 ********************************************************************************************************** //

static struct Test_Case_String_Inputs set_Test_Case_String_Inputs ()
{
	assert(sizeof(struct Test_Case_String_Inputs) == 3*STRLEN_MIN*sizeof(char));

	struct Test_Case_String_Inputs tcsi;
	tcsi.num_flux_1st[0] = 0;
	tcsi.num_flux_2nd[0] = 0;
	tcsi.test_norm[0] = 0;
	return tcsi;
}

static void set_string_associations_test_case (struct Test_Case_T* test_case, const struct Test_Case_String_Inputs* tcsi)
{
/// \todo add string associations for solver types.
	// num_flux_1st
	if (strcmp(tcsi->num_flux_1st,"upwind") == 0)
		const_cast_i(&test_case->ind_num_flux[0],NUM_FLUX_UPWIND);
	else if (strcmp(tcsi->num_flux_1st,"Roe-Pike") == 0)
		const_cast_i(&test_case->ind_num_flux[0],NUM_FLUX_ROE_PIKE);
	else
		const_cast_i(&test_case->ind_num_flux[0],NUM_FLUX_INVALID);

	// num_flux_2nd
	if (strcmp(tcsi->num_flux_2nd,"BR2") == 0)
		const_cast_i(&test_case->ind_num_flux[1],NUM_FLUX_BR2);
	else
		const_cast_i(&test_case->ind_num_flux[1],NUM_FLUX_INVALID);

	// test_norm
	if (strcmp(tcsi->test_norm,"H0") == 0)
		const_cast_i(&test_case->ind_test_norm,TEST_NORM_H0);
	else if (strcmp(tcsi->test_norm,"H1") == 0)
		const_cast_i(&test_case->ind_test_norm,TEST_NORM_H1);
	else if (strcmp(tcsi->test_norm,"H1_upwind") == 0)
		const_cast_i(&test_case->ind_test_norm,TEST_NORM_H1_UPWIND);
	else
		const_cast_i(&test_case->ind_test_norm,TEST_NORM_INVALID);
}

static void correct_invalid_test_case_parameters (struct Test_Case_T* test_case, const struct Simulation* sim)
{
	switch (sim->method) {
	case METHOD_DG:
		const_cast_b(&test_case->use_schur_complement,false);
		break;
	case METHOD_DPG:
		break; // Do nothing.
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

static const bool* get_compute_member_Flux_Input
	(const char type_ei, struct Test_Case_T*const test_case, const struct Simulation*const sim)
{
	static const bool cm_100000[] = {1,0,0,0,0,0,},
	                  cm_110000[] = {1,1,0,0,0,0,},
	                  cm_110100[] = {1,1,0,1,0,0,};

	assert(type_ei == 'e' || type_ei == 'i');
	switch (test_case->pde_index) {
	case PDE_ADVECTION:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_110000;
		break;
	case PDE_EULER:
		switch (sim->method) {
		case METHOD_DG:
			if (type_ei == 'e')
				return cm_100000;
			else if (type_ei == 'i')
				return cm_110000;
			break;
		case METHOD_DPG:
			if (type_ei == 'e')
				return cm_110000;
			else if (type_ei == 'i')
				return cm_110100;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",sim->method);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
	EXIT_ERROR("Should not have reached this point.\n");
}

static const bool* get_compute_member_Boundary_Value_Input
	(const char type_ei, struct Test_Case_T*const test_case, const struct Simulation*const sim)
{
	UNUSED(sim);
	static const bool cm_100000[] = {1,0,0,0,0,0,},
	                  cm_110000[] = {1,1,0,0,0,0,};

	assert(type_ei == 'e' || type_ei == 'i');
	switch (test_case->pde_index) {
	case PDE_ADVECTION:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_110000;
		break;
	case PDE_EULER:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_110000;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
	EXIT_ERROR("Should not have reached this point.\n");
}
