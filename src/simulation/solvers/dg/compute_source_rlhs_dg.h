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

#ifndef DPG__compute_source_rlhs_dg_h__INCLUDED
#define DPG__compute_source_rlhs_dg_h__INCLUDED
/** \file
 *  \brief Provides functions used for computing the source contributions to the right-hand side (rhs) term of the DG
 *         scheme.
 */

struct Simulation;

/** \brief Compute the source contributions to the rhs term for the DG scheme.
 *
 *  \warning Should not be used for non-linear sources.
 *
 *  Updates:
 *  - \ref DG_Solver_Volume_T::rhs.
 */
void compute_source_rhs_dg
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Compute the contribution of the source term integrals to the flux imbalances for the DG scheme.
void compute_flux_imbalances_source_dg
	(const struct Simulation*const sim ///< \ref Simulation.
	);

#endif // DPG__compute_source_rlhs_dg_h__INCLUDED
