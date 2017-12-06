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
#include "solution_advection_default.h"
#include "peterson/solution_peterson.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_advection_T.c"

struct Sol_Data__Advection get_sol_data_advection (const struct Simulation* sim)
{
	static bool need_input = true;

	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(sim->input_path,&sol_data);
	}

	return sol_data;
}

void read_data_advection (const char*const input_path, struct Sol_Data__Advection*const sol_data)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input(input_path,'s'); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"b_adv")) {
			read_skip_d_1(line,1,sol_data->b_adv,DMAX);
			++count_found;
		}
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
