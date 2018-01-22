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

#include "const_cast.h"
#include "file_processing.h"
#include "simulation.h"
#include "test_case.h"

#include "test_complex_solution.h"
#include "test_complex_solution_advection.h"
#include "test_complex_solution_diffusion.h"
#include "test_complex_solution_euler.h"
#include "test_complex_test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "test_case_T.c"

void convert_to_Test_Case_rc (struct Simulation* sim, const char type_rc_o)
{
	struct Test_Case_rc* test_case_rc = sim->test_case_rc;

	switch (type_rc_o) {
	case 'c':
		assert(test_case_rc->is_real == true);
		destructor_Test_Case(test_case_rc->tc);

		const_cast_b(&test_case_rc->is_real,false);
		test_case_rc->tc = (void*)constructor_Test_Case_c(sim); // keep
		break;
	case 'r':
		assert(test_case_rc->is_real == false);
		destructor_Test_Case_c(test_case_rc->tc);

		const_cast_b(&test_case_rc->is_real,true);
		test_case_rc->tc = (void*)constructor_Test_Case(sim); // keep
		break;
	default:
		EXIT_ERROR("Unsupported: %c.\n",type_rc_o);
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
