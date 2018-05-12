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
 *  \brief Provides the templated Euler solution functions.
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_numerical_flux.h"
#include "definitions_physics.h"


#include "def_templates_solution_euler.h"

#include "def_templates_multiarray.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_geometry.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_euler_T (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	test_case->set_grad = set_sg_do_nothing_T;
	if (strstr(sim->pde_spec,"periodic_vortex")) {
		test_case->constructor_sol              = constructor_const_sol_invalid_T;
		test_case->set_sol                      = set_sol_periodic_vortex_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_euler_all;
		const_cast_b(&test_case->has_analytical,true);
	} else if (strstr(sim->pde_spec,"supersonic_vortex")) {
		test_case->constructor_xyz              = constructor_xyz_cylinder_parametric_T;
		test_case->constructor_sol              = constructor_const_sol_supersonic_vortex_T;
		test_case->set_sol                      = set_sol_supersonic_vortex_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;

		if (sim->method == METHOD_DG) {
			const_cast_b(&test_case->copy_initial_rhs,true);
			test_case->constructor_Error_CE         = constructor_Error_CE_euler_all_p_rhs;
		} else { // not yet supported
			test_case->constructor_Error_CE         = constructor_Error_CE_euler_all;
		}
test_case->constructor_Error_CE = constructor_Error_CE_euler_all;
		test_case->constructor_Error_CE_restart_test = constructor_Error_CE_euler_all;
		const_cast_b(&test_case->has_analytical,true);
	} else if (strstr(sim->pde_spec,"free_stream")) {
		test_case->constructor_xyz              = constructor_xyz_trigonometric_cube_parametric_xl_T;
		test_case->constructor_sol              = constructor_const_sol_free_stream_T;
		test_case->set_sol                      = set_sol_free_stream_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_euler_all;
	} else if (strstr(sim->pde_spec,"joukowski")) {
		test_case->constructor_xyz              = constructor_xyz_joukowski_parametric_T;
		test_case->constructor_sol              = constructor_const_sol_free_stream_T;
		test_case->set_sol                      = set_sol_free_stream_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_euler_entropy;
		test_case->constructor_Error_CE_functionals = constructor_Error_CE_functionals__cl;
	} else if (strstr(sim->pde_spec,"gaussian_bump")) {
		test_case->constructor_xyz              = constructor_xyz_gaussian_bump_parametric_T;
		test_case->constructor_sol              = constructor_const_sol_free_stream_T;
		test_case->set_sol                      = set_sol_free_stream_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		if (sim->method == METHOD_DG) {
			const_cast_b(&test_case->copy_initial_rhs,true);
			test_case->constructor_Error_CE = constructor_Error_CE_euler_entropy_p_rhs;
		} else { // not yet supported
			test_case->constructor_Error_CE = constructor_Error_CE_euler_entropy;
		}
test_case->constructor_Error_CE = constructor_Error_CE_euler_entropy;
test_case->constructor_Error_CE = constructor_Error_CE_euler_all_p_rhs;
		test_case->constructor_Error_CE_restart_test = constructor_Error_CE_euler_entropy;
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	test_case->compute_Flux = compute_Flux_1_T;
	test_case->compute_Flux_iv[0] = compute_Flux_T_euler;
	test_case->compute_Flux_iv[1] = NULL;

	test_case->compute_Numerical_Flux = compute_Numerical_Flux_1_T;
	switch (test_case->ind_num_flux[0]) {
	case NUM_FLUX_ROE_PIKE:
		test_case->compute_Numerical_Flux_e[0] = compute_Numerical_Flux_T_euler_roe_pike;
		test_case->compute_Numerical_Flux_i[0] = compute_Numerical_Flux_T_euler_roe_pike_jacobian;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",test_case->ind_num_flux[0]);
		break;
	}

	test_case->constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_s_fcl_interp_T;
}

void convert_variables_T (struct Multiarray_T* vars, const char type_i, const char type_o)
{
	assert(vars->layout == 'C');
	assert(vars->extents[1] == NVAR);

	const ptrdiff_t ext_0 = vars->extents[0];


	switch (type_i) {
	case 'p': {
		Type* rho = get_col_Multiarray_T(0,vars),
		    * p   = get_col_Multiarray_T(NVAR-1,vars),
		    * E   = p;

/// \todo Use ARRAY_DIM here.
		Type* uvw[DMAX] = {            get_col_Multiarray_T(1,vars),
		                    (DIM > 1 ? get_col_Multiarray_T(2,vars) : NULL),
		                    (DIM > 2 ? get_col_Multiarray_T(3,vars) : NULL), };

		Type* rhouvw[DMAX] = { uvw[0], uvw[1], uvw[2], };
		switch (type_o) {
		case 'c':
			for (ptrdiff_t i = 0; i < ext_0; ++i) {
				Type V2 = 0.0;
				for (int d = 0; d < DIM; ++d) {
					V2 += uvw[d][i]*uvw[d][i];
					rhouvw[d][i] = rho[i]*uvw[d][i];
				}
				E[i] = p[i]/GM1 + 0.5*rho[i]*V2;
			}
			break;
		case 'p': // fallthrough
		default:
			EXIT_ERROR("Unsupported: %c\n",type_o);
			break;
		}
		break;
	} case 'c': {
		Type* rho = get_col_Multiarray_T(0,vars),
		    * p   = get_col_Multiarray_T(NVAR-1,vars),
		    * E   = p;

		Type* uvw[DMAX] = {            get_col_Multiarray_T(1,vars),
		                    (DIM > 1 ? get_col_Multiarray_T(2,vars) : NULL),
		                    (DIM > 2 ? get_col_Multiarray_T(3,vars) : NULL), };

		Type* rhouvw[DMAX] = { uvw[0], uvw[1], uvw[2], };
		switch (type_o) {
		case 'p':
			for (ptrdiff_t i = 0; i < ext_0; ++i) {
				const Type rho_inv = 1.0/rho[i];
				Type rho2V2 = 0.0;
				for (int d = 0; d < DIM; ++d) {
					rho2V2 += rhouvw[d][i]*rhouvw[d][i];
					uvw[d][i] = rhouvw[d][i]*rho_inv;
				}
				p[i] = GM1*(E[i]-0.5*rho2V2*rho_inv);
			}
			break;
		case 'c': // fallthrough
		default:
			EXIT_ERROR("Unsupported: %c\n",type_o);
			break;
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %c\n",type_i);
		break;
	}
}

const struct const_Multiarray_T* constructor_const_functionals_cd_cl_zero_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	UNUSED(sim);
	const ptrdiff_t n_n = xyz->extents[0];

	struct Multiarray_T* func = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,2}); // returned
	set_to_value_Multiarray_T(func,0.0);

	return (struct const_Multiarray_T*) func;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
