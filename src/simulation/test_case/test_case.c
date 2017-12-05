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

#include "test_case.h"

#include "const_cast.h"
#include "file_processing.h"
#include "simulation.h"
#include "solution_advection.h"
#include "solution_euler.h"

// Static function declarations ************************************************************************************* //

/// Container for input strings which are to be subsequently converted to integer parameters.
struct Test_Case_String_Inputs {
	char num_flux_1st[STRLEN_MIN]; ///< The name of the 1st order numerical flux scheme to be used.
	char num_flux_2nd[STRLEN_MIN]; ///< The name of the 2nd order numerical flux scheme to be used.

	char test_norm[STRLEN_MIN]; ///< The name of the norm to use for the optimal test function computation.
};

/** \brief Return a statically allocated \ref Test_Case_String_Inputs container which zero-initialized members.
 *  \return See brief. */
static struct Test_Case_String_Inputs set_Test_Case_String_Inputs
	();

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_flux.h"
#include "def_templates_test_case.h"
#include "test_case_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Test_Case_String_Inputs set_Test_Case_String_Inputs ()
{
	assert(sizeof(struct Test_Case_String_Inputs) == 3*STRLEN_MIN*sizeof(char));

	struct Test_Case_String_Inputs tcsi;
	tcsi.num_flux_1st[0] = 0;
	tcsi.num_flux_2nd[0] = 0;
	tcsi.test_norm[0] = 0;
	return tcsi;
}

