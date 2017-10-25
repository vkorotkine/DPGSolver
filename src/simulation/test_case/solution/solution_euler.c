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

#include "solution_euler.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "compute_error.h"
#include "const_cast.h"
#include "flux.h"
#include "numerical_flux.h"
#include "simulation.h"
#include "solution.h"
#include "solution_periodic_vortex.h"
//#include "solution_supersonic_vortex.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_euler (struct Test_Case* test_case, const struct Simulation*const sim)
{
	test_case->set_grad = set_sg_do_nothing;
	if (strstr(sim->pde_spec,"periodic_vortex")) {
		test_case->set_sol = set_sol_periodic_vortex;
		test_case->compute_source = compute_source_do_nothing;
		test_case->constructor_Error_CE = constructor_Error_CE_euler_all;
	} else if (strstr(sim->pde_spec,"supersonic_vortex")) {
		EXIT_ADD_SUPPORT;
//		test_case->set_sol_v = set_sol_v_supersonic_vortex;
//		test_case->set_sol_f = set_sol_f_supersonic_vortex;
//		test_case->compute_source = compute_source_do_nothing;
//		test_case->compute_errors = compute_errors_euler_all;
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	const bool* flux_comp_mem_e = (bool[]){1,0,0},
	          * flux_comp_mem_i = (bool[]){1,1,0};
	for (int i = 0; i < MAX_NUM_FLUX_OUT; ++i) {
		const_cast_b(&test_case->flux_comp_mem_e[i],flux_comp_mem_e[i]);
		const_cast_b(&test_case->flux_comp_mem_i[i],flux_comp_mem_i[i]);
	}

	test_case->compute_Flux = compute_Flux_1;
	test_case->compute_Flux_e[0] = compute_Flux_euler;
	test_case->compute_Flux_e[1] = NULL;
	test_case->compute_Flux_i[0] = compute_Flux_euler_jacobian;
	test_case->compute_Flux_i[1] = NULL;

	test_case->compute_Numerical_Flux = compute_Numerical_Flux_1;
	switch (test_case->ind_num_flux[0]) {
	case NUM_FLUX_ROE_PIKE:
		test_case->compute_Numerical_Flux_e[0] = compute_Numerical_Flux_euler_roe_pike;
//		test_case->compute_Numerical_Flux_i[0] = compute_Numerical_Flux_euler_roe_pike_jacobian;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",test_case->ind_num_flux[0]);
		break;
	}

	test_case->constructor_s_l_fcl = constructor_s_l_fcl_interp;
	test_case->constructor_g_l_fcl = constructor_sg_fc_null;
}

void convert_variables (struct Multiarray_d* vars, const char type_i, const char type_o)
{
	assert(type_i != type_o);
	assert(vars->layout == 'C');

	const ptrdiff_t ext_0 = vars->extents[0],
	                n_var = vars->extents[1];

	const ptrdiff_t d = n_var-2;

	switch (type_i) {
	case 'p': {
		double* rho = get_col_Multiarray_d(0,vars),
		      * u   = get_col_Multiarray_d(1,vars),
		      * v   = (d > 1 ? get_col_Multiarray_d(2,vars) : NULL),
		      * w   = (d > 2 ? get_col_Multiarray_d(3,vars) : NULL),
		      * p   = get_col_Multiarray_d(n_var-1,vars);

		double* rhou = u,
		      * rhov = v,
		      * rhow = w,
		      * E    = p;
		switch (type_o) {
		case 'c':
			for (ptrdiff_t i = 0; i < ext_0; ++i) {
				double V2 = 0.0;
				switch (d) {
				case 3:
					V2 += w[i]*w[i];
					rhow[i] = rho[i]*w[i];
					// fallthrough
				case 2:
					V2 += v[i]*v[i];
					rhov[i] = rho[i]*v[i];
					// fallthrough
				case 1:
					V2 += u[i]*u[i];
					rhou[i] = rho[i]*u[i];
					E[i] = p[i]/GM1 + 0.5*rho[i]*V2;
					break;
				default:
					EXIT_ERROR("Unsupported: %td\n",d);
					break;
				}
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
		double* rho  = get_col_Multiarray_d(0,vars),
		      * rhou = get_col_Multiarray_d(1,vars),
		      * rhov = (d > 1 ? get_col_Multiarray_d(2,vars) : NULL),
		      * rhow = (d > 2 ? get_col_Multiarray_d(3,vars) : NULL),
		      * E    = get_col_Multiarray_d(n_var-1,vars);

		double* u = rhou,
		      * v = rhov,
		      * w = rhow,
		      * p = E;
		switch (type_o) {
		case 'p':
			for (ptrdiff_t i = 0; i < ext_0; ++i) {
				const double rho_inv = 1.0/rho[i];
				double rho2V2 = 0.0;
				switch (d) {
				case 3:
					rho2V2 += rhow[i]*rhow[i];
					w[i] = rhow[i]*rho_inv;
					// fallthrough
				case 2:
					rho2V2 += rhov[i]*rhov[i];
					v[i] = rhov[i]*rho_inv;
					// fallthrough
				case 1:
					rho2V2 += rhou[i]*rhou[i];
					u[i] = rhou[i]*rho_inv;
					p[i] = GM1*(E[i]-0.5*rho2V2*rho_inv);
					break;
				default:
					EXIT_ERROR("Unsupported: %td\n",d);
					break;
				}
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

void compute_entropy (struct Multiarray_d* s, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(s->extents[0] == vars->extents[0]);
	assert(s->extents[1] == 1);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const double* rho = get_col_const_Multiarray_d(0,vars),
	            * p   = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = s->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		s->data[i] = p[i]*pow(rho[i],-GAMMA);

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

void compute_mach (struct Multiarray_d* mach, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(mach->extents[0] == vars->extents[0]);
	assert(mach->extents[1] == 1);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const ptrdiff_t d = vars->extents[1]-2;

	const double* rho = get_col_const_Multiarray_d(0,vars),
		      * u   = get_col_const_Multiarray_d(1,vars),
		      * v   = (d > 1 ? get_col_const_Multiarray_d(2,vars) : NULL),
		      * w   = (d > 2 ? get_col_const_Multiarray_d(3,vars) : NULL),
	            * p   = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = mach->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		double V2 = 0.0;
		switch (d) {
			case 3: V2 += w[i]*w[i]; // fallthrough
			case 2: V2 += v[i]*v[i]; // fallthrough
			case 1: V2 += u[i]*u[i]; break;
			default: EXIT_ERROR("Unsupported: %td\n",d); break;
		}
		const double c2 = GAMMA*p[i]/rho[i];
		mach->data[i] = sqrt(V2/c2);
	}

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
