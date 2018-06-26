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
#include "definitions_numerical_flux.h"
#include "definitions_test_case.h"


#include "def_templates_solution_advection.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_geometry.h"
#include "def_templates_math_functions.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the function pointers to the numerical flux computing functions.
static void set_function_pointers_num_flux_T
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_advection_T (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	const_cast_b(&test_case->has_analytical,true);
	const_cast_b(&test_case->copy_initial_rhs,true);
	test_case->set_grad = set_sg_do_nothing_T;
	if (strstr(sim->pde_spec,"peterson")) {
		test_case->constructor_xyz              = NULL;
		test_case->constructor_sol              = constructor_const_sol_peterson_T;
		test_case->set_sol                      = set_sol_peterson_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all;
	} else if ((strcmp(sim->pde_spec,"demkowicz_dpg_ii") == 0) ||
	           (strcmp(sim->pde_spec,"steady/default") == 0)) {
		test_case->constructor_xyz              = NULL;
		test_case->constructor_sol              = constructor_const_sol_advection_default_T;
		test_case->set_sol                      = set_sol_advection_default_T;
		test_case->compute_source_rhs           = compute_source_rhs_advection_default_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_advection_default_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all;
	} else if (strstr(sim->pde_spec,"steady/vortex")) {
		if (strstr(sim->geom_name,"n-cube")) {
			test_case->constructor_xyz = constructor_xyz_trigonometric_cube_parametric_xl_oct1_T;
		} else if (strstr(sim->geom_name,"n-cylinder_hollow_section")) {
			test_case->constructor_xyz = constructor_xyz_cylinder_parametric_T;
		} else {
			EXIT_ADD_SUPPORT;
		}
		test_case->constructor_sol              = constructor_const_sol_vortex_advection_T;
		test_case->set_sol                      = set_sol_vortex_advection_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;

		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all;
		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all_p_rhs;
	} else if (strstr(sim->pde_spec,"steady/free_stream")) {
		if (strstr(sim->geom_name,"n-cube")) {
			test_case->constructor_xyz = constructor_xyz_trigonometric_cube_parametric_xl_T;
		} else if (strstr(sim->geom_name,"n-cylinder_hollow_section")) {
			test_case->constructor_xyz = constructor_xyz_cylinder_parametric_T;
		} else {
			EXIT_ADD_SUPPORT;
		}
		test_case->constructor_sol              = constructor_const_sol_free_stream_advection_T;
		test_case->set_sol                      = set_sol_free_stream_advection_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all;
	} else if (strstr(sim->pde_spec,"unsteady/hyperbolic_tan")) {
		test_case->constructor_xyz              = constructor_xyz_fixed_cube_parametric_T;
		test_case->constructor_sol              = constructor_const_sol_hyperbolic_tan_T;
		test_case->set_sol                      = set_sol_hyperbolic_tan_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all;
		const_cast_b(&test_case->has_analytical,true);
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	if (sim->method == METHOD_DPG) {
		printf("*** Warning: Disabling computation of rhs_0 error for DPG. Add support.***\n");
		const_cast_b(&test_case->copy_initial_rhs,false);
		test_case->constructor_Error_CE = constructor_Error_CE_advection_all;
	}

	test_case->compute_Flux = compute_Flux_1_T;
	test_case->compute_Flux_iv[0] = compute_Flux_T_advection;

	set_function_pointers_num_flux_T(test_case,sim);

	test_case->constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_s_fcl_interp_T;
}

struct Sol_Data__Advection_T get_sol_data_advection_T ( )
{
	static bool need_input = true;

	static struct Sol_Data__Advection_T sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection_T(&sol_data);
	}

	return sol_data;
}

void read_data_advection_T (struct Sol_Data__Advection_T*const sol_data)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input('s',NULL,NULL); // closed

	int advection_type = 0;

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_convert_i(line,"advection_type",&advection_type,&count_found);
		if (strstr(line,"u_scale"))
			read_skip_d_1(line,1,&sol_data->u_scale,1);
		if (strstr(line,"use_constant_solution"))
			read_skip_b(line,&sol_data->use_constant_solution);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");

	switch (advection_type) {
		case ADVECTION_TYPE_CONST:  sol_data->compute_b_adv = compute_b_adv_constant_T; break;
		case ADVECTION_TYPE_VORTEX: sol_data->compute_b_adv = compute_b_adv_vortex_T;   break;
		default:                    EXIT_ERROR("Unsupported: %d\n",advection_type);     break;
	}
}

const Real* compute_b_adv_constant_T (const Type*const xyz)
{
	UNUSED(xyz);
	static bool need_input = true;
	static Real b_adv[DIM] = {0,};

	if (need_input) {
		need_input = false;

		const int count_to_find = 1;

		FILE* input_file = fopen_input('s',NULL,NULL); // closed

		int count_found = 0;
		char line[STRLEN_MAX];
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"b_adv")) {
				read_skip_d_1(line,1,b_adv,DIM);
				++count_found;
			}
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	return b_adv;
}

const Real* compute_b_adv_vortex_T (const Type*const xyz)
{
	static bool need_input = true;
	static Real b_mag = 0;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;

		FILE* input_file = fopen_input('s',NULL,NULL); // closed

		int count_found = 0;
		char line[STRLEN_MAX];
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"b_magnitude")) {
				read_skip_d_1(line,1,&b_mag,1);
				++count_found;
			}
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	static Real b_adv[DIM] = {0,};
	assert(DIM == 2);

	const Real x = real_T(xyz[0]),
	           y = real_T(xyz[1]);
MAYBE_UNUSED(x);

	// Denominator omitted to allow for exact GCL satisfaction.
//	const Real t = atan2(y,x);
//	IF_DIM_GE_1( b_adv[0] =  b_mag*sin(t); ) // Note:  sin(atan2(y,x)) ==  y/sqrt(x^2+y^2).
//	IF_DIM_GE_2( b_adv[1] = -b_mag*cos(t); ) //       -cos(atan2(y,x)) == -x/sqrt(x^2+y^2).
	IF_DIM_GE_1( b_adv[0] =  b_mag*y; )
	IF_DIM_GE_2( b_adv[1] = -b_mag*x; )

	return b_adv;
}

/* Additional 2D Divergence free advection fields:
 * - doublet: b_adv[0] =  (y*y-x*x)/pow(x*x+y*y,2); b_adv[1] = -(2*x*y)/pow(x*x+y*y,2);
 * - other: b_adv[0] = -x*x; b_adv[1] = 2*x*y;
 */

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_function_pointers_num_flux_T (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	switch (sim->method) {
	case METHOD_DG:  // fallthrough
	case METHOD_DPG: // fallthrough
	case METHOD_OPG: // fallthrough
	case METHOD_OPGC0:
		test_case->compute_Numerical_Flux = compute_Numerical_Flux_1_T;
		switch (test_case->ind_num_flux[0]) {
		case NUM_FLUX_UPWIND:
			test_case->compute_Numerical_Flux_e[0] = compute_Numerical_Flux_T_advection_upwind;
			test_case->compute_Numerical_Flux_i[0] = compute_Numerical_Flux_T_advection_upwind_jacobian;
			break;
		default:
			EXIT_ERROR("Unsupported: %d.\n",test_case->ind_num_flux[0]);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

#include "undef_templates_solution_advection.h"

#include "undef_templates_boundary.h"
#include "undef_templates_flux.h"
#include "undef_templates_geometry.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_solution.h"
#include "undef_templates_test_case.h"
