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
 *  \brief Provides the templated Navier-Stokes solution functions.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_numerical_flux.h"
#include "definitions_physics.h"


#include "def_templates_solution_navier_stokes.h"

#include "def_templates_multiarray.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_geometry.h"
#include "def_templates_math_functions.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

/// \brief Set the function pointers to the numerical flux computing functions.
static void set_function_pointers_num_flux
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_navier_stokes_T (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	if (strstr(sim->pde_spec,"taylor_couette")) {
		test_case->constructor_xyz              = constructor_xyz_cylinder_parametric_T;
		test_case->constructor_sol              = constructor_const_sol_taylor_couette_T;
		test_case->constructor_grad             = constructor_const_grad_taylor_couette_T;
		test_case->set_sol                      = set_sol_taylor_couette_T;
		test_case->set_grad                     = set_grad_taylor_couette_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_navier_stokes_uvwt;
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	test_case->compute_Flux = compute_Flux_12_T;
	test_case->compute_Flux_iv[0] = compute_Flux_T_euler;
	test_case->compute_Flux_iv[1] = compute_Flux_T_navier_stokes;

	test_case->compute_Numerical_Flux = compute_Numerical_Flux_12_T;
	set_function_pointers_num_flux(test_case,sim);

	test_case->constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_s_fcl_interp_T;
}

void convert_variables_gradients_T
	(struct Multiarray_T*const grad, const struct const_Multiarray_T*const sol, const char type_i, const char type_o)
{
	assert(type_i != type_o);
	assert(grad->layout == 'C');
	assert(sol->layout == 'C');
	assert(grad->extents[0] == sol->extents[0]);
	assert(grad->extents[1] == sol->extents[1]);
	assert(grad->extents[1] == NVAR);
	assert(grad->extents[2] == DIM);

	const ptrdiff_t ext_0 = grad->extents[0];

	switch (type_i) {
	case 'p': {
		const Type*const rho   = get_col_const_Multiarray_T(0,sol),
		          *const uvw[] = ARRAY_DIM( get_col_const_Multiarray_T(1,sol),
		                                    get_col_const_Multiarray_T(2,sol),
		                                    get_col_const_Multiarray_T(3,sol) );

		for (int d_g = 0; d_g < DIM; ++d_g) {
			Type*const g_rho   = get_col_Multiarray_T(d_g*NVAR+0,grad),
			    *const g_uvw[] = ARRAY_DIM( get_col_Multiarray_T(d_g*NVAR+1,grad),
			                                get_col_Multiarray_T(d_g*NVAR+2,grad),
			                                get_col_Multiarray_T(d_g*NVAR+3,grad) ),
			    *const g_p     = get_col_Multiarray_T(d_g*NVAR+NVAR-1,grad),

			    *const*const g_rhouvw = g_uvw,
			    *const g_E            = g_p;

			switch (type_o) {
			case 'c':
				for (ptrdiff_t i = 0; i < ext_0; ++i) {
					Type V2   = 0.0,
					     g_V2 = 0.0;
					for (int d = 0; d < DIM; ++d) {
						V2   += uvw[d][i]*uvw[d][i];
						g_V2 += 2.0*(uvw[d][i]*g_uvw[d][i]);

						g_rhouvw[d][i] = g_rho[i]*uvw[d][i] + rho[i]*g_uvw[d][i];
					}
					g_E[i] = g_p[i]/GM1 + 0.5*(g_rho[i]*V2 + rho[i]*g_V2);
				}
				break;
			case 'p': // fallthrough
			default:
				EXIT_ERROR("Unsupported: %c\n",type_o);
				break;
			}
		}
		break;
	} case 'c': {
		const Type*const rho      = get_col_const_Multiarray_T(0,sol),
		          *const rhouvw[] = ARRAY_DIM( get_col_const_Multiarray_T(1,sol),
		                                       get_col_const_Multiarray_T(2,sol),
		                                       get_col_const_Multiarray_T(3,sol) );

		for (int d_g = 0; d_g < DIM; ++d_g) {
			Type*const g_rho      = get_col_Multiarray_T(d_g*NVAR+0,grad),
			    *const g_rhouvw[] = ARRAY_DIM( get_col_Multiarray_T(d_g*NVAR+1,grad),
			                                   get_col_Multiarray_T(d_g*NVAR+2,grad),
			                                   get_col_Multiarray_T(d_g*NVAR+3,grad) ),
			    *const g_E        = get_col_Multiarray_T(d_g*NVAR+NVAR-1,grad),

			    *const*const g_uvw = g_rhouvw,
			    *const g_p         = g_E;

			switch (type_o) {
			case 'p':
				for (ptrdiff_t i = 0; i < ext_0; ++i) {
					Type rho2V2   = 0.0,
					     g_rho2V2 = 0.0;
					for (int d = 0; d < DIM; ++d) {
						rho2V2   += rhouvw[d][i]*rhouvw[d][i];
						g_rho2V2 += 2.0*(rhouvw[d][i]*g_rhouvw[d][i]);

						g_uvw[d][i] = 1.0/(rho[i]*rho[i])*(rho[i]*g_rhouvw[d][i] - g_rho[i]*rhouvw[d][i]);
					}
					g_p[i] = GM1*(g_E[i] - 0.5/(rho[i]*rho[i])*(rho[i]*g_rho2V2 - g_rho[i]*rho2V2));
				}
				break;
			case 'c': // fallthrough
			default:
				EXIT_ERROR("Unsupported: %c\n",type_o);
				break;
			}
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %c\n",type_i);
		break;
	}
}

compute_mu_fptr_T get_compute_mu_fptr_T (const char*const input_path)
{
	static int viscosity_type = VISCOSITY_INVALID;
	static bool need_input = true;
	set_viscosity_type_T(input_path,&viscosity_type,&need_input);

	switch (viscosity_type) {
		case VISCOSITY_CONSTANT:   return compute_mu_constant_T;                  break;
		case VISCOSITY_SUTHERLAND: return compute_mu_sutherland_T;                break;
		default:                   EXIT_ERROR("Unsupported: %d.",viscosity_type); break;
	};
}

void set_viscosity_type_T
	(const char*const input_path, int*const viscosity_type_ptr, bool*const need_input)
{
	if (*need_input) {
		*need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_convert_i(line,"viscosity_type",viscosity_type_ptr,&count_found);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}
}

Type compute_mu_constant_T (const char*const input_path, const Type rho, const Type*const rhouvw, const Type E)
{
	UNUSED(rho); UNUSED(rhouvw); UNUSED(E);
	static Real mu = 0.0;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_string_count_d("mu",&count_found,line,&mu);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}
	return mu;
}

Type compute_mu_sutherland_T (const char*const input_path, const Type rho, const Type*const rhouvw, const Type E)
{
	static Real r_s = 0.0;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_string_count_d("r_s",&count_found,line,&r_s);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	const Type V2 = compute_V2_from_rhouvw_T(rho,rhouvw),
	           p = GM1*(E-0.5*rho*V2),
		     T = p/(rho*r_s);

	static const Real c1 = 1.46e-6,
	                  c2 = 112;

	return c1/(1+c2/T)*sqrt_T(T);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_function_pointers_num_flux (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	// Inviscid flux
	switch (test_case->ind_num_flux[0]) {
	case NUM_FLUX_ROE_PIKE:
		test_case->compute_Numerical_Flux_e[0] = compute_Numerical_Flux_T_euler_roe_pike;
		test_case->compute_Numerical_Flux_i[0] = compute_Numerical_Flux_T_euler_roe_pike_jacobian;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",test_case->ind_num_flux[0]);
		break;
	}

	// Viscous flux
	switch (sim->method) {
	case METHOD_DG: // fallthrough
	case METHOD_DPG:
		switch (test_case->ind_num_flux[1]) {
		case NUM_FLUX_BR2_STABLE: // fallthrough
		case NUM_FLUX_CDG2:
			test_case->compute_Numerical_Flux_e[1] = compute_Numerical_Flux_T_navier_stokes_central;
			test_case->compute_Numerical_Flux_i[1] = compute_Numerical_Flux_T_navier_stokes_central_jacobian;
			break;
		default:
			EXIT_ERROR("Unsupported: %d.\n",test_case->ind_num_flux[1]);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}
