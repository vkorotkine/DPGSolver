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

#ifndef DPG__test_case_h__INCLUDED
#define DPG__test_case_h__INCLUDED
/** \file
 *  \brief Provides container(s) and functions relating to the test cases.
 */

#include "flux.h"
#include "geometry.h"
#include "numerical_flux.h"
#include "solution.h"
#include "compute_error.h"

#include "def_templates_type_d.h"
#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_geometry.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"
#include "test_case_T.h"
#include "undef_templates_type.h"
#include "undef_templates_boundary.h"
#include "undef_templates_flux.h"
#include "undef_templates_geometry.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_solution.h"
#include "undef_templates_test_case.h"

/** \brief Container providing an additional layer of dereferencing to either a 'r'eal or 'c'omplex \ref Test_Case_T.
 *  This container is provided so that \ref Simulation does not need to be templated.
 */
struct Test_Case_rc {
	const bool is_real; ///< Flag for whether \ref Test_Case_rc::test_case is real (as opposed to complex).

	void* tc; ///< Pointer to the \ref Test_Case_T.
};

/** \brief Constructor for a real \ref Test_Case_rc.
 *  \return See brief. */
struct Test_Case_rc* constructor_Test_Case_rc_real
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a real \ref Test_Case_rc.
void destructor_Test_Case_rc_real
	(struct Test_Case_rc* test_case_rc ///< \ref Test_Case_rc.
	);

/** \brief Check whether the test case under consideration requires positivity of certain variables.
 *  \return `true` if yes; `false` otherwise. */
bool test_case_requires_positivity
	(const struct Test_Case*const test_case ///< \ref Test_Case_T.
	);

#endif // DPG__test_case_h__INCLUDED
