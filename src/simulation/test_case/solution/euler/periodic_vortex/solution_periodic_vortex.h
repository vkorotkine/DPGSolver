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

#ifndef DPG__solution_periodic_vortex_h__INCLUDED
#define DPG__solution_periodic_vortex_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to set the solution for the periodic vortex test case.
 */

struct Simulation;
struct Solution_Container;

/// \brief Function to be used for \ref Test_Case::set_sol for the periodic vortex test case.
void set_sol_periodic_vortex
	(const struct Simulation* sim,      ///< Defined for \ref set_sol_fptr.
	 struct Solution_Container sol_cont ///< Defined for \ref set_sol_fptr.
	);

#endif // DPG__solution_periodic_vortex_h__INCLUDED
