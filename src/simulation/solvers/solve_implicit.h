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

#ifndef DPG__solve_implicit_h__INCLUDED
#define DPG__solve_implicit_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to solve for the solution using implicit procedures.
 *
 *  \todo Update this.
 */

struct Simulation;

/// \brief Solve for the solution using an implicit solver.
void solve_implicit
	(struct Simulation* sim ///< \ref Simulation.
	);

#endif // DPG__solve_implicit_h__INCLUDED
