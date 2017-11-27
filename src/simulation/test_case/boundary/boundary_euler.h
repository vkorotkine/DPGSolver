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

#ifndef DPG__boundary_euler_h__INCLUDED
#define DPG__boundary_euler_h__INCLUDED
/** \file
 *  \brief Provides containers and functions relating to boundary conditions for the Euler equation.
 */

struct Boundary_Value_Input;
struct Boundary_Value;
struct Solver_Face;
struct Simulation;

/** \brief Version of \ref constructor_Boundary_Value_fptr computing members using the Riemann invariant values.
 *  \return See brief.
 *
 *  Reference: (Section 2.2, \cite Carlson2011); note the typo in eq. (14).
 */
void constructor_Boundary_Value_euler_riemann
	(struct Boundary_Value* bv,               ///< See brief.
	 const struct Boundary_Value_Input* bv_i, ///< See brief.
	 const struct Solver_Face* face,          ///< See brief.
	 const struct Simulation* sim             ///< See brief.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr computing members using the slip wall values.
 *  \return See brief.
 *
 *  The slip wall boundary condition sets the ghost state density and total energy equal to the internal values and uses
 *  the opposite normal velocity.
 */
void constructor_Boundary_Value_euler_slipwall
	(struct Boundary_Value* bv,               ///< See brief.
	 const struct Boundary_Value_Input* bv_i, ///< See brief.
	 const struct Solver_Face* face,          ///< See brief.
	 const struct Simulation* sim             ///< See brief.
	);

#endif // DPG__boundary_euler_h__INCLUDED
