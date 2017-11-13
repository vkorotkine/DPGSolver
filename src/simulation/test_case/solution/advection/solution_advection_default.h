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

#ifndef DPG__solution_advection_default_h__INCLUDED
#define DPG__solution_advection_default_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to set the default solution for linear advection test cases.
 */

struct Simulation;
struct Solution_Container;
struct Solver_Volume;
struct Multiarray_d;

/// \brief Function to be used for \ref Test_Case::set_sol for the default linear advection solution.
void set_sol_advection_default
	(const struct Simulation* sim,      ///< Defined for \ref set_sol_fptr.
	 struct Solution_Container sol_cont ///< Defined for \ref set_sol_fptr.
	);

/** \brief Function to be used for \ref Test_Case::constructor_sol for the default linear advection solution.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_const_sol_advection_default
	(const struct const_Multiarray_d* xyz, ///< Defined for \ref constructor_sol_fptr.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr.
	);

/// \brief Version of \ref compute_source_rhs_fptr for the default linear advection solution.
void compute_source_rhs_advection_default
	(const struct Simulation* sim,      ///< See brief.
	 const struct Solver_Volume* s_vol, ///< See brief.
	 struct Multiarray_d* rhs           ///< See brief.
	);

#endif // DPG__solution_advection_default_h__INCLUDED
