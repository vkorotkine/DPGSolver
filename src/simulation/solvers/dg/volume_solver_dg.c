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

#include "element_solver_dg.h"
#include "volume.h"
#include "volume_solver.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// Container holding flags for which members of \ref DG_Solver_Volume are needed.
struct Needed_Members {
	bool sol_coef_p, ///< Flag for \ref DG_Solver_Volume::sol_coef_p.
	     m_inv;      ///< Flag for \ref DG_Solver_Volume::m_inv.
};

/** \brief Return a statically allocated \ref Needed_Members container with values set.
 *  \return See brief. */
static struct Needed_Members set_needed_members
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_matrix_d.h"
#include "def_templates_multiarray_d.h"
#include "def_templates_vector_d.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"
#include "volume_solver_dg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Needed_Members set_needed_members (const struct Simulation* sim)
{
	const struct Test_Case* test_case = sim->test_case;
	struct Needed_Members needed_members =
		{ .sol_coef_p = false,
		  .m_inv      = false, };

	switch (test_case->solver_proc) {
	case SOLVER_E: // fallthrough
	case SOLVER_EI:
		if (!sim->collocated)
			needed_members.m_inv = true;
		switch (test_case->solver_type_e) {
		case SOLVER_E_SSP_RK_33: // fallthrough
		case SOLVER_E_LS_RK_54:
			needed_members.sol_coef_p = true;
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

	return needed_members;
}
