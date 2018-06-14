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

#include "test_case_c.h"

#include "boundary.h"
#include "flux.h"
#include "geometry.h"
#include "numerical_flux.h"
#include "solution.h"
#include "compute_error.h"

#include "def_templates_type_d.h"
#include "test_case_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "test_case_T.h"
#include "undef_templates_type.h"

/** \brief Container providing an additional layer of dereferencing to either a 'r'eal or 'c'omplex \ref Test_Case_T.
 *  This container is provided so that \ref Simulation does not need to be templated.
 */
struct Test_Case_rc {
	const bool is_real; ///< Flag for whether \ref Test_Case_rc::tc is real (as opposed to complex).

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

/** \brief Check whether the test case explicitly enforces conservation.
 *  \return `true` if yes; `false` otherwise. */
bool test_case_explicitly_enforces_conservation
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Return whether the solution is being initialized from a restart file.
 *  \return See brief. */
bool using_restart ( );

/** \brief Return whether a restart file should be outputted.
 *  \return See brief. */
bool outputting_restart ( );

/** \brief Return a statically allocated array of `int*` holding the values of the number of variables/equations.
 *  \return See brief.
 *
 *  n_variables = out[0]
 *  n_equations = out[1]
 *
 *  Passing a non-NULL input for `new_vals` sets the statically allocated array to the contained values.
 */
const int* get_set_n_var_eq
	(const int*const new_vals ///< New values.
	);

/** \brief Return a statically allocated array of `bool*` holding the values of the flags for whether the PDE under
 *         consideration has 1st and/or 2nd order components.
 *  \return See brief.
 *
 *  Passing a non-NULL input for `new_vals` sets the statically allocated array to the contained values.
 */
const bool* get_set_has_1st_2nd_order
	(const bool*const new_vals ///< New values.
	);

/** \brief Return a statically allocated `int` holding the value of the pde index.
 *  \return See brief.
 *
 *  Passing a non-NULL input for `new_val` sets the statically allocated int to the value.
 */
int get_set_pde_index
	(const int*const new_val ///< New value.
		);

#endif // DPG__test_case_h__INCLUDED
