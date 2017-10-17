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

#ifndef DPG__compute_volume_rlhs_dg_h__INCLUDED
#define DPG__compute_volume_rlhs_dg_h__INCLUDED
/** \file
 *  \brief Provides functions used for computing the volume contributions to the right and left-hand sides (rlhs) terms
 *         of the DG scheme.
 */

struct Simulation;

/** \brief Compute the volume contributions to the rhs terms for the DG scheme.
 *
 *  Computes:
 *  - \ref DG_Solver_Volume::rhs.
 */
void compute_volume_rhs_dg
	(const struct Simulation* sim ///< \ref Simulation.
	);

#endif // DPG__compute_volume_rlhs_dg_h__INCLUDED
