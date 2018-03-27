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

#include "solution_advection.h"

#include <math.h>

#include "macros.h"
#include "definitions_solution.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error.h"
#include "compute_error_advection.h"
#include "const_cast.h"
#include "file_processing.h"
#include "flux_advection.h"
#include "numerical_flux_advection.h"
#include "simulation.h"
#include "solution.h"
#include "test_case.h"

#include "solution_advection_default.h"
#include "free_stream_advection/solution_free_stream_advection.h"
#include "peterson/solution_peterson.h"
#include "vortex_advection/solution_vortex_advection.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_advection_T.c"

struct Sol_Data__Advection get_sol_data_advection ( )
{
	static bool need_input = true;

	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(&sol_data);
	}

	return sol_data;
}

void read_data_advection (struct Sol_Data__Advection*const sol_data)
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
		if (strstr(line,"u_coef_polynomial4"))
			read_skip_d_1(line,1,sol_data->u_coef_polynomial4,5);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");

	switch (advection_type) {
		case ADVECTION_TYPE_CONST:  sol_data->compute_b_adv = compute_b_adv_constant; break;
		case ADVECTION_TYPE_VORTEX: sol_data->compute_b_adv = compute_b_adv_vortex;   break;
		default:                    EXIT_ERROR("Unsupported: %d\n",advection_type);   break;
	}
}

const double* compute_b_adv_constant (const double*const xyz)
{
	UNUSED(xyz);
	static bool need_input = true;
	static double b_adv[DIM] = {0,};

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

const double* compute_b_adv_vortex (const double*const xyz)
{
	static bool need_input = true;
	static double b_mag = 0;
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

	static double b_adv[DIM] = {0,};
	assert(DIM == 2);

	const double x = xyz[0],
	             y = xyz[1],
	             t = atan2(y,x);

	IF_DIM_GE_1( b_adv[0] =  sin(t); )
	IF_DIM_GE_2( b_adv[1] = -cos(t); )

	return b_adv;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
