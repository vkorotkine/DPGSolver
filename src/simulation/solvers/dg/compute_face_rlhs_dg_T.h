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
 *         terms of the DG scheme.
 *
 *  The convention for the notation for the cross-terms in the linearized equations is as follows. The contributions on
 *  either side of the face are named according to their interpretation as being on the 'l'eft or on the 'r'ight of the
 *  face. Let the solution variables in the global system be denoted by \f$ s \f$. Then the Jacobian terms in the
 *  appropriate sub-range of the linearized matrix will be
 *  \f[
 *  \begin{bmatrix}
 *  	L_{ll} & L_{lr} \\
 *  	L_{rl} & L_{rr} \\
 *  \end{bmatrix}
 *  \begin{bmatrix}
 *  	s_{l} \\
 *  	s_{r} \\
 *  \end{bmatrix}
 *  =
 *  -
 *  \begin{bmatrix}
 *  	R_{l} \\
 *  	R_{r} \\
 *  \end{bmatrix}
 *  \f]
 *
 *  where \f$ L_{xy} := \frac{\partial R_{x}}{\partial s_{y}} \f$ (i.e. \f$ L_{xy} \f$ accounts for the effect of
 *  the change of the \f$ s_y \f$ (variables) on \f$ R_x \f$ (residuals)). The choice of convention may be more
 *  intuitive when reading the left-hand-side of the above equation from right to left.
 *
 *  In terminology used elsewhere, \f$ y \f$ is the affector and \f$ x \f$ is the affectee.
 */

struct Simulation;
struct Solver_Storage_Implicit;
struct Face;

/** \brief Compute the face contributions to the rhs (and optionally lhs) terms for the DG scheme.
 *
 *  \todo update this list when complete.
 *  Computes:
 */
void compute_face_rlhs_dg_T
	(const struct Simulation* sim,             ///< \ref Simulation.
	 struct Solver_Storage_Implicit* s_store_i ///< \ref Solver_Storage_Implicit.
	);
