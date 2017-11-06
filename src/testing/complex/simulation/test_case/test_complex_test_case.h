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

#ifndef DPG__test_complex_test_case_h__INCLUDED
#define DPG__test_complex_test_case_h__INCLUDED
/** \file
 *  \brief Provides `complex` version of container(s) and functions defined in \ref test_case.h.
 */

#include <stdbool.h>

#include "test_case.h"
#include "test_complex_boundary.h"

struct Simulation;

/// \brief Derived `complex` version of \ref Test_Case.
struct Complex_Test_Case {
	struct Test_Case test_case; ///< Base \ref Test_Case.

	constructor_Boundary_Value_Input_c_face_fptr constructor_Boundary_Value_Input_c_face_fcl; ///< See brief.
};

/** \brief Constructor for a derived \ref Complex_Test_Case.
 *  \return See brief. */
void constructor_derived_Complex_Test_Case
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Complex_Test_Case.
void destructor_derived_Complex_Test_Case
	(struct Simulation* sim ///< \ref Simulation.
	);

#endif // DPG__test_complex_test_case_h__INCLUDED
