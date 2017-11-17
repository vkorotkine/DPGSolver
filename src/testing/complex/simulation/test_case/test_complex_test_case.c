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

#include "test_complex_test_case.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the function pointer members of \ref Test_Case.
static void set_function_pointers
	(struct Complex_Test_Case* test_case,     ///< \ref Complex_Test_Case.
	 const struct Test_Case*const test_case_b ///< \ref Test_Case.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Complex_Test_Case (struct Simulation* sim)
{
	struct Complex_Test_Case* test_case = calloc(1,sizeof *test_case); // free

	memcpy(test_case,sim->test_case,sizeof(struct Test_Case)); // shallow copy of the base.

	set_function_pointers(test_case,sim->test_case);

	destructor_Test_Case(sim->test_case);
	sim->test_case = (struct Test_Case*) test_case;
}

void destructor_derived_Complex_Test_Case (struct Simulation* sim)
{
	struct Test_Case* test_case_b = calloc(1,sizeof *test_case_b); // moved
	memcpy(test_case_b,sim->test_case,sizeof(struct Test_Case)); // shallow copy of the base.

	free((void*)sim->test_case);
	sim->test_case = test_case_b;
}

bool has_complex_Jacobians (const int method)
{
	switch (method) {
	case METHOD_DG:
		return false;
		break;
	case METHOD_DPG:
		return true;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",method);
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_function_pointers (struct Complex_Test_Case* test_case, const struct Test_Case*const test_case_b)
{
	// Boundary_Value_Input
	if (test_case_b->constructor_Boundary_Value_Input_face_fcl == constructor_Boundary_Value_Input_face_s_fcl_interp)
		test_case->constructor_Boundary_Value_Input_c_face_fcl = constructor_Boundary_Value_Input_c_face_s_fcl_interp;
	else
		EXIT_UNSUPPORTED;
}
