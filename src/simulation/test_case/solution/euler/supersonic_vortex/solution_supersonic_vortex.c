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

#include "solution_supersonic_vortex.h"

#include "multiarray.h"

#include "element_solution.h"
#include "volume_solver.h"

#include "file_processing.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solution data relating to 's'upersonic 'v'ortex.
struct Sol_Data__sv {
	// Read parameters
	double r_i,   ///< The 'r'adius of the 'i'nternal cylinder.
	       m_i,   ///< The 'm'ach number at the 'i'nternal cylinder face.
	       rho_i, ///< The density (\f$ rho \f$) at the 'i'nternal cylinder face.
	       V_i;   ///< The magnitude of the 'V'elocity at the 'i'nternal cylinder face.
};

/** \brief Return the statically allocated \ref Sol_Data__sv container.
 *  \return See brief. */
static struct Sol_Data__sv get_sol_data
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Read the required solution data into \ref Sol_Data__sv.
static void read_data_supersonic_vortex
	(const char*const input_path,       ///< Defined in \ref fopen_input.
	 struct Sol_Data__sv*const sol_data ///< \ref Sol_Data__sv.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_supersonic_vortex_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Sol_Data__sv get_sol_data (const struct Simulation* sim)
{
	static bool need_input = true;

	static struct Sol_Data__sv sol_data;
	if (need_input) {
		need_input = false;
		read_data_supersonic_vortex(sim->input_path,&sol_data);
	}

	return sol_data;
}

static void read_data_supersonic_vortex (const char*const input_path, struct Sol_Data__sv*const sol_data)
{
	const int count_to_find = 4;

	FILE* input_file = fopen_input(input_path,'s'); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("r_i",  &count_found,line,&sol_data->r_i);
		read_skip_string_count_d("m_i",  &count_found,line,&sol_data->m_i);
		read_skip_string_count_d("rho_i",&count_found,line,&sol_data->rho_i);
		read_skip_string_count_d("V_i",  &count_found,line,&sol_data->V_i);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
