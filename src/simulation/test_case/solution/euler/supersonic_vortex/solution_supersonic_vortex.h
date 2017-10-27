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

#ifndef DPG__solution_supersonic_vortex_h__INCLUDED
#define DPG__solution_supersonic_vortex_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to set the solution for the supersonic vortex test case.
 */

struct Simulation;
struct Solver_Volume;
struct Solver_Face;

/// \brief \todo update this.
void set_sol_coef_v_supersonic_vortex
	(const struct Simulation* sim, ///< \todo update this.
	 struct Solver_Volume* volume  ///< \todo update this.
	);

/// \brief \todo update this.
void set_sol_coef_f_supersonic_vortex
	(const struct Simulation* sim, ///< \todo update this.
	 struct Solver_Face* face      ///< \todo update this.
	);

#endif // DPG__solution_supersonic_vortex_h__INCLUDED
