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
 *  \brief Provides the static templated boundary condition functions.
 */

#include "def_templates_multiarray.h"
#include "def_templates_test_case.h"

/** \brief Constructor for the Multiarray_\* holding the solution at the given xyz coordinates.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_sol_bv
	(const struct const_Multiarray_d* xyz, ///< xyz coordinates.
	 const struct Simulation* sim          ///< \ref Simulation.
	)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	return test_case->constructor_sol(xyz,sim); // returned
}

#include "undef_templates_multiarray.h"
#include "undef_templates_test_case.h"
