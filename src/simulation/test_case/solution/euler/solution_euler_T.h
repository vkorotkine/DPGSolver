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
 *  \brief Provides templated functions relating to the Euler solutions.
 */

#include "def_templates_solution_euler.h"
#include "def_templates_multiarray.h"
#include "def_templates_test_case.h"

struct Test_Case_T;
struct Simulation;
struct Multiarray_T;
struct const_Multiarray_T;

/// \brief Set the solution function pointer members of an Euler \ref Test_Case_T.
void set_function_pointers_solution_euler_T
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Convert between supported variable types.
 *
 *  Supported types include:
 *  - 'p'rimitive:    [rho u v w p].
 *  - 'c'onservative: [rho*[1 u v w] E]; E = p/(GAMMA-1) + 0.5*rho*V^2.
 */
void convert_variables_T
	(struct Multiarray_T* vars, ///< The container holding the data.
	 const char type_i,         ///< The input variable type.
	 const char type_o          ///< The output variable type.
	);

#include "undef_templates_solution_euler.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_test_case.h"
