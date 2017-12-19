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

#include "compute_error.h"
#include "test_complex_flux.h"
#include "test_complex_geometry.h"
#include "test_complex_numerical_flux.h"
#include "test_complex_solution.h"

#include "def_templates_type_dc.h"
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

struct Simulation;

/// \brief Convert a \ref Test_Case_T from `real` to `complex` or vice versa.
void convert_to_Test_Case_rc
	(struct Simulation* sim, ///< \ref Simulation.
	 const char type_rc_o    ///< The output 'r'eal/'c'omplex type.
	);

#endif // DPG__test_complex_test_case_h__INCLUDED
