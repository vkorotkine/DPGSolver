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

#include "test_integration_convergence.h"

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_intrusive.h"
#include "definitions_visualization.h"

#include "test_base.h"
#include "test_integration.h"

#include "computational_elements.h"

#include "compute_error.h"
#include "geometry.h"
#include "simulation.h"
#include "solve.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

/// \brief 'Con'/'De'strucor for a \ref Simulation depending on `adapt_type`.
static void structor_simulation
	(struct Simulation** sim,   ///< Pointer to the \ref Simulation.
	 const char mode,           ///< Mode of operation. Options: 'c'onstructor, 'd'estructor.
	 const int adapt_type,      ///< \ref Integration_Test_Info::adapt_type.
	 const int p,               ///< The order of the current simulation.
	 const int ml,              ///< The mesh level of the current simulation.
	 const int p_prev,          ///< The order of the previous simulation.
	 const int ml_prev,         ///< The mesh level of the previous simulation.
	 const char*const ctrl_name ///< The name of the control file (used for \ref constructor_Simulation).
	);

/** \brief Set the name of the current control file to be used.
 *  \return See brief. */
const char* set_ctrl_name_curr
	(const int adapt_type,      ///< \ref Integration_Test_Info::adapt_type.
	 const int p,               ///< The order of the current simulation.
	 const int ml,              ///< The mesh level of the current simulation.
	 const char*const ctrl_name ///< The name of the control file (used for \ref constructor_Simulation).
	);

// Interface functions ********************************************************************************************** //

void test_integration_convergence (struct Test_Info*const test_info, const char*const ctrl_name)
{

//	constructor_derived_computational_elements(sim,IL_SOLVER);

	const char* ctrl_name_full = set_ctrl_name_full(ctrl_name);
	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name_full);

	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	const int adapt_type = int_test_info->adapt_type;

	struct Simulation* sim = NULL;
	for (int p = p_ref[0], p_prev = p; p <= p_ref[1]; ++p) {
	for (int ml = ml_ref[0], ml_prev = ml; ml <= ml_ref[1]; ++ml) {
p = 2;
ml = 2;


		const char*const ctrl_name_curr = set_ctrl_name_curr(adapt_type,p,ml,ctrl_name);
		structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr);

		output_visualization(sim,VIS_GEOM_EDGES);
		output_visualization(sim,VIS_SOLUTION);

		output_error(sim);

		p_prev  = p;
		ml_prev = ml;
		structor_simulation(&sim,'d',adapt_type,p,ml,p_prev,ml_prev,NULL);
	}}

	destructor_Integration_Test_Info(int_test_info);

UNUSED(test_info);
EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void structor_simulation
	(struct Simulation** sim, const char mode, const int adapt_type, const int p, const int ml, const int p_prev,
	 const int ml_prev, const char*const ctrl_name)
{
	assert(mode == 'c' || mode == 'd');

	switch (adapt_type) {
	case ADAPT_0:
		if (mode == 'c') {
// Group these calls in a function when working.
			*sim = constructor_Simulation(ctrl_name); // destructed
			constructor_derived_computational_elements(*sim,IL_SOLVER); // destructed
			set_up_solver_geometry(*sim);
			solve_for_solution(*sim);
		} else if (mode == 'd') {
			destructor_derived_computational_elements(*sim,IL_BASE);
			destructor_Simulation(*sim);
		}
		break;
	case ADAPT_P:
		EXIT_ADD_SUPPORT;
		if (ml != ml_prev) {
			structor_simulation(sim,mode,ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name);
		} else {
			; // p-adapt
		}
		break;
	case ADAPT_H:
		EXIT_ADD_SUPPORT;
		if (p != p_prev) {
			structor_simulation(sim,mode,ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name);
		} else {
			; // h-adapt
		}
		break;
	case ADAPT_HP:
		EXIT_ADD_SUPPORT;
		if (ml != ml_prev)
			{ ; } // h-adapt
		if (p != p_prev)
			{ ; } // p-adapt
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}
}

const char* set_ctrl_name_curr (const int adapt_type, const int p, const int ml, const char*const ctrl_name)
{
	static char ctrl_name_curr[STRLEN_MAX] = { 0, };
	strcpy(ctrl_name_curr,ctrl_name);

	assert(p  >= 0 && p  <= 9); // Constraints are only imposed because of the currently limited string processing.
	assert(ml >= 0 && ml <= 9);

	char* index = NULL;
	switch (adapt_type) {
	case ADAPT_0:
		index = strstr(ctrl_name_curr,"__ml");
		index[4] = '0'+ml;
		assert(!isdigit(index[5]));
		index = strstr(ctrl_name_curr,"__p");
		index[3] = '0'+p;
		assert(!isdigit(index[4]));
		break;
	case ADAPT_P:
		assert(strstr(ctrl_name_curr,"__ml") == NULL);
		index = strstr(ctrl_name_curr,"__p");
		index[3] = '0'+p;
		assert(!isdigit(index[4]));
		break;
	case ADAPT_H:
		assert(strstr(ctrl_name_curr,"__p") == NULL);
		index = strstr(ctrl_name_curr,"__ml");
		index[4] = '0'+ml;
		assert(!isdigit(index[5]));
		break;
	case ADAPT_HP:
		assert(strstr(ctrl_name_curr,"__ml") == NULL);
		assert(strstr(ctrl_name_curr,"__p") == NULL);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}

	return ctrl_name_curr;
}
