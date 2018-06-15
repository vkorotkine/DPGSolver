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
 *  \brief Provides the interface to templated functions used to set the solution for the taylor-couette test case
 *  (Section 7, \cite Illingworth1950).
 */

#include "def_templates_solution.h"
#include "def_templates_solution_navier_stokes.h"
#include "def_templates_multiarray.h"

struct const_Multiarray_T;
struct Simulation;
struct Solution_Container_T;

/// \brief Function to be used for \ref Test_Case_T::set_sol for the taylor-couette test case.
void set_sol_taylor_couette_T
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/// \brief Function to be used for \ref Test_Case_T::set_grad for the taylor-couette test case.
void set_grad_taylor_couette_T
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/** \brief Function to be used for \ref Test_Case_T::constructor_sol for the taylor-couette test case.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_sol_taylor_couette_T
	(const struct const_Multiarray_T* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
	);

/** \brief Function to be used for \ref Test_Case_T::constructor_grad for the taylor-couette test case.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_grad_taylor_couette_T
	(const struct const_Multiarray_T* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
	);

#include "undef_templates_solution.h"
#include "undef_templates_solution_navier_stokes.h"
#include "undef_templates_multiarray.h"
