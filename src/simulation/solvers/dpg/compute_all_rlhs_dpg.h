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

#ifndef DPG__compute_all_rlhs_dpg_h__INCLUDED
#define DPG__compute_all_rlhs_dpg_h__INCLUDED
/** \file
 *  \brief Provides functions used for all contributions to the right and left-hand side (rlhs) terms of the DPG
 *         scheme.
 */

struct Simulation;
struct Solver_Storage_Implicit;
struct Volume;

/// \brief Compute all contributions to the rhs and lhs terms for the DPG scheme.
void compute_all_rlhs_dpg
	(const struct Simulation* sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi ///< \ref Solver_Storage_Implicit.
	);

#endif // DPG__compute_all_rlhs_dpg_h__INCLUDED
