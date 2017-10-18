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

#ifndef DPG__solve_h__INCLUDED
#define DPG__solve_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to solve for the solution.
 */

struct Simulation;

/// \brief Solve for the solution.
void solve_for_solution
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the volume and face rhs terms for the the given method.
 *  \return The maximum absolute value of the rhs.
 *
 *  \todo Make the math pretty.
 *  The rhs includes all terms of the discretization except for the time-varying term and **does not** include the
 *  inverse of the mass matrix (i.e. M_v d/dt sol_coef = rhs).
 */
double compute_rhs
	(const struct Simulation* sim ///< \ref Simulation.
	);

#endif // DPG__solve_h__INCLUDED
