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

#ifndef DPG__solution_euler_h__INCLUDED
#define DPG__solution_euler_h__INCLUDED
/** \file
 *  \brief Provides real functions relating to the Euler solutions.
 */

#include "def_templates_type_d.h"
#include "def_templates_solution_euler.h"
#include "def_templates_multiarray.h"
#include "def_templates_test_case.h"
#include "solution_euler_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solution_euler.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_test_case.h"

struct Multiarray_d;
struct const_Multiarray_d;

/// \brief Compute the entropy measure: s = p/pow(rho,GAMMA).
void compute_entropy
	(struct Multiarray_d* s,                ///< The container to hold the entropy data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

/// \brief Compute the mach number: mach = V/c.
void compute_mach
	(struct Multiarray_d* mach,             ///< The container to hold the mach data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

#endif // DPG__solution_euler_h__INCLUDED
