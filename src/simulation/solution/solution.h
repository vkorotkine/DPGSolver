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

#ifndef DPG__solution_h__INCLUDED
#define DPG__solution_h__INCLUDED
/**	\file
 *	\brief Provides the interface to functions used for solution specification.
 */

struct Simulation;
struct Intrusive_List;

/** \brief Set up the initial solution for the simulation. Computes:
 *	- \ref Solver_Volume::sol_coef;
 *	- \ref Solver_Volume::grad_coef (if applicable).
 */
void set_up_solution
	(struct Simulation* sim,               ///< \ref Simulation.
	 struct Intrusive_List* solver_volumes ///< The \ref Solver_Volume list for which to set up the solution.
	);

#endif // DPG__solution_h__INCLUDED
