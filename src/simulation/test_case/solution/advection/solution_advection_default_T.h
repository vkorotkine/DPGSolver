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
 *  \brief Provides the interface to functions used to set the default solution for linear advection test cases.
 */

struct Simulation;
struct Solution_Container_T;
struct Solver_Volume_T;
struct Multiarray_T;
struct const_Multiarray_R;

/// \brief Function to be used for \ref Test_Case_T::set_sol for the default linear advection solution.
void set_sol_advection_default_T
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_fptr.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr.
	);

/** \brief Function to be used for \ref Test_Case_T::constructor_sol for the default linear advection solution.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_sol_advection_default_T
	(const struct const_Multiarray_R* xyz, ///< Defined for \ref constructor_sol_fptr.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr.
	);

/// \brief Version of \ref compute_source_rhs_fptr for the default linear advection solution.
void compute_source_rhs_advection_default_T
	(const struct Simulation* sim,        ///< See brief.
	 const struct Solver_Volume_T* s_vol, ///< See brief.
	 struct Multiarray_T* rhs             ///< See brief.
	);
