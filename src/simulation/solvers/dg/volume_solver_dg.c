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

#include "volume_solver_dg.h"

#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Check if the memory allocation for \ref DG_Solver_Volume::sol_coef_p is needed.
 *  \return `true` if needed. */
static bool check_need_sol_coef_p
	(const struct Test_Case* test_case ///< \ref Test_Case.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Volume (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Solver_Volume* s_volume  = (struct Solver_Volume*) volume_ptr;
	struct DG_Solver_Volume* volume = (struct DG_Solver_Volume*) volume_ptr;

	const int order = s_volume->sol_coef->order;
	ptrdiff_t* extents = s_volume->sol_coef->extents;

	volume->rhs        = constructor_empty_Multiarray_d('C',order,extents); // destructed
	volume->sol_coef_p = // destructed
		( check_need_sol_coef_p(sim->test_case) ? constructor_empty_Multiarray_d('C',order,extents) : NULL );
}

void destructor_derived_DG_Solver_Volume (struct Volume* volume_ptr)
{
	struct DG_Solver_Volume* volume = (struct DG_Solver_Volume*) volume_ptr;

	destructor_Multiarray_d(volume->rhs);
	if (volume->sol_coef_p != NULL)
		destructor_Multiarray_d(volume->sol_coef_p);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static bool check_need_sol_coef_p (const struct Test_Case* test_case)
{
	bool allocate_sol_coef_p = false;

	switch (test_case->solver_proc) {
	case SOLVER_E: // fallthrough
	case SOLVER_EI:
		switch (test_case->solver_type_e) {
		case SOLVER_E_SSP_RK_33: // fallthrough
		case SOLVER_E_LS_RK_54:
			allocate_sol_coef_p = true;
			break;
		case SOLVER_E_EULER:
			// Do nothing
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",test_case->solver_type_e);
			break;
		}
		break;
	case SOLVER_I:
		// Do nothing
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}

	return allocate_sol_coef_p;
}
