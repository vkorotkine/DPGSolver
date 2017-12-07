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

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

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

void set_function_pointers_solution_euler (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	test_case->set_grad = set_sg_do_nothing_T;
	if (strstr(sim->pde_spec,"periodic_vortex")) {
		test_case->constructor_sol      = constructor_const_sol_invalid_T;
		test_case->set_sol              = set_sol_periodic_vortex_T;
		test_case->compute_source_rhs   = compute_source_rhs_do_nothing_T;
		test_case->constructor_Error_CE = constructor_Error_CE_euler_all;
	} else if (strstr(sim->pde_spec,"supersonic_vortex")) {
		test_case->constructor_xyz      = constructor_xyz_cylinder_parametric_T;
		test_case->constructor_sol      = constructor_const_sol_supersonic_vortex_T;
		test_case->set_sol              = set_sol_supersonic_vortex_T;
		test_case->compute_source_rhs   = compute_source_rhs_do_nothing_T;
		test_case->constructor_Error_CE = constructor_Error_CE_euler_all;
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	const bool* flux_comp_mem_e = (bool[]){1,0,0},
	          * flux_comp_mem_i = (bool[]){1,1,0};
	for (int i = 0; i < MAX_NUM_FLUX_OUT; ++i) {
		const_cast_b(&test_case->flux_comp_mem_e[i],flux_comp_mem_e[i]);
		const_cast_b(&test_case->flux_comp_mem_i[i],flux_comp_mem_i[i]);
	}

	test_case->compute_Flux = compute_Flux_1_T;
	test_case->compute_Flux_e[0] = compute_Flux_T_euler;
	test_case->compute_Flux_e[1] = NULL;
	test_case->compute_Flux_i[0] = compute_Flux_T_euler_jacobian;
	test_case->compute_Flux_i[1] = NULL;

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
	assert(type_i != type_o);
	assert(vars->layout == 'C');

	const ptrdiff_t ext_0 = vars->extents[0];

	assert(vars->extents[1] == NVAR);

	switch (type_i) {
	case 'p': {
		Type* rho = get_col_Multiarray_T(0,vars),
		    * p   = get_col_Multiarray_T(NVAR-1,vars),
		    * E   = p;

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
		case 'p':
			return;
			break;
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
		case 'c':
			return;
			break;
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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
