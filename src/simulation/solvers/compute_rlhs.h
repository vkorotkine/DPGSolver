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

#ifndef DPG__compute_rlhs_h__INCLUDED
#define DPG__compute_rlhs_h__INCLUDED
/** \file
 *  \brief Provides functions used for computing the general contributions to the right and left-hand side (rlhs) terms
 *         of supported schemes.
 */

#include <stdbool.h>
#include "definitions_core.h"

struct Matrix_d;
struct const_Matrix_d;
struct Simulation;

/** \brief Compute the source contributions to the rhs term as required by the DG scheme.
 *
 *  The function is "dg-like" in the sense that it computes the rhs term as for the DG scheme, but can be used with test
 *  space not equal to trial space (Petrov-Galerkin).
 *
 *  Updates:
 *  - \ref Solver_Volume_T::rhs.
 */
void compute_source_rhs_dg_like
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Compute the maximum absolute value of the rhs term stored in \ref Solver_Volume_T::rhs.
 *  \return See brief. */
double compute_max_rhs_dg_like
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Add a contribution to the d_grad_coef__d_s_coef matrix, the 'r'ight matrix in the 'p'artial 'l'eft-'h'and
 *         's'ide entry for a volume/face.
 *
 *  This function is used to add any necessary terms to partially or fully corrected weak gradient linearization terms.
 */
void add_to_lhs_p_r
	(const double alpha,                             ///< Scaling constant.
	 const struct const_Matrix_d*const dgc_dsc[DIM], ///< The linearization term.
	 struct Matrix_d*const lhs_p_r,                  ///< The 'r'ight 'p'artial term for the lhs matrix.
	 const bool boundary_face_term                   ///< Flag for whether the entry is from a boundary face.
	);

#endif // DPG__compute_rlhs_h__INCLUDED
