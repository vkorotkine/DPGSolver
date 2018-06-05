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
 *  \brief Provides templated functions used for computing the face contributions to the right and left-hand side (rlhs)
 *         terms of the OPG scheme.
 *
 *  Please see the comments in \ref compute_face_rlhs_dg_T.h for the notation convention for the cross-terms in the
 *  linearized equations.
 */

struct Simulation;
struct Solver_Storage_Implicit;
struct Intrusive_List;

/// \brief Compute the face contributions to the rhs (and optionally lhs) terms for the DG scheme.
void compute_face_rlhs_opg_T
	(const struct Simulation*const sim,        ///< Standard.
	 struct Solver_Storage_Implicit*const ssi, ///< Standard.
	 struct Intrusive_List*const faces         ///< The list of faces.
	);
