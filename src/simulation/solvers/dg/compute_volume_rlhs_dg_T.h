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
 *  \brief Provides templated functions used for computing the volume contributions to the right and left-hand side
 *         (rlhs) terms of the DG scheme.
 */

struct Simulation;
struct Solver_Storage_Implicit;
struct Intrusive_List;

/// \brief Compute the volume contributions to the rhs (and optionally lhs) terms for the DG scheme.
void compute_volume_rlhs_dg_T
	(const struct Simulation* sim,        ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi, ///< \ref Solver_Storage_Implicit.
	 struct Intrusive_List* volumes       ///< The list of volumes.
	);
