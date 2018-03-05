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

#include "solution_navier_stokes.h"

#include <float.h>

#include "multiarray.h"

#include "boundary.h"
#include "compute_error_euler.h"
#include "compute_error_navier_stokes.h"
#include "file_processing.h"
#include "flux_euler.h"
#include "flux_navier_stokes.h"
#include "geometry.h"
#include "geometry_parametric.h"
#include "numerical_flux_euler.h"
#include "numerical_flux_navier_stokes.h"
#include "simulation.h"
#include "solution.h"
#include "solution_euler.h"
#include "test_case.h"

#include "free_stream/solution_free_stream.h"
#include "taylor_couette/solution_taylor_couette.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_navier_stokes_T.c"

double compute_kappa_const_cp (const double mu, const double Cp, const double Pr)
{
	return mu*Cp/Pr;
}

double compute_cp_ideal_gas (const double r_s)
{
	return GAMMA/GM1*r_s;
}

void compute_viscosity (struct Multiarray_d* mu, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(mu->extents[0] == vars->extents[0]);
	assert(mu->extents[1] == 1);
	assert(vars->extents[1]-2 == DIM);
	assert(vars->layout == 'C');

	compute_mu_fptr_d compute_mu = get_compute_mu_fptr_d();

	if (var_type != 'c')
		convert_variables((struct Multiarray_d*)vars,var_type,'c');

	const double*const rho         = get_col_const_Multiarray_d(0,vars),
	            *const rhouvw[DIM] = ARRAY_DIM( get_col_const_Multiarray_d(1,vars),
		                                      get_col_const_Multiarray_d(2,vars),
		                                      get_col_const_Multiarray_d(3,vars) ),
	            *const E           = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = mu->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		const double rhouvw_i[DIM] = ARRAY_DIM( rhouvw[0][i], rhouvw[1][i], rhouvw[2][i] );
		mu->data[i] = compute_mu(rho[i],rhouvw_i,E[i]);
	}

	if (var_type != 'c')
		convert_variables((struct Multiarray_d*)vars,'c',var_type);
}

double get_normal_flux_Energy ( )
{
	static double nf_E = DBL_MIN;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		char line[STRLEN_MAX];
		int count_to_find = 0;
		int count_found = 0;
		FILE* input_file = NULL;

		int diabatic_flux_type = VISCOUS_BC_INVALID;

		count_to_find = 1;
		count_found = 0;
		input_file = fopen_input('s',NULL,NULL); // closed
		while (fgets(line,sizeof(line),input_file))
			read_skip_convert_i(line,"diabatic_flux_type",&diabatic_flux_type,&count_found);
		fclose(input_file);
		assert(count_found == count_to_find);

		if (diabatic_flux_type == DIABATIC_FLUX_CONSTANT_ZERO) {
			nf_E = 0.0;
		} else {
			assert(diabatic_flux_type == DIABATIC_FLUX_CONSTANT);

			int viscosity_type = VISCOSITY_INVALID;
			double Pr   = DBL_MAX,
			       r_s  = DBL_MAX,
			       mu   = DBL_MAX,
			       dTdn = DBL_MAX;

			count_to_find = 5;
			count_found = 0;
			input_file = fopen_input('s',NULL,NULL); // closed
			while (fgets(line,sizeof(line),input_file)) {
				read_skip_convert_i(line,"viscosity_type",&viscosity_type,&count_found);

				read_skip_string_count_d("Pr",  &count_found,line,&Pr);
				read_skip_string_count_d("r_s", &count_found,line,&r_s);
				read_skip_string_count_d("mu",  &count_found,line,&mu);
				read_skip_string_count_d("dTdn",&count_found,line,&dTdn);
			}
			fclose(input_file);
			assert(count_found == count_to_find);

			assert(viscosity_type == VISCOSITY_CONSTANT); // Otherwise nf_E still varies for constant dTdn.

			const double Cp = compute_cp_ideal_gas(r_s);

			/// See comments in \ref flux_navier_stokes_T.h regarding the negation used below.
			nf_E = (-1.0)*mu*Cp/Pr*dTdn;

			assert(nf_E != 0.0);
		}
	}
	return nf_E;
}

double get_r_s ( )
{
	static double r_s = 0.0;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input('s',NULL,NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_string_count_d("r_s",&count_found,line,&r_s);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	assert(r_s != 0.0);
	return r_s;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
