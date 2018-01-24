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
 *  \brief Provides the interface for the templated \ref DG_Solver_Face_T container and associated functions.
 *
 *  These faces are needed by the 'D'iscontinuous 'G'alerkin solver functions.
 */

#include "definitions_core.h"

/// \brief Container for data relating to the DG solver faces.
struct DG_Solver_Face_T {
	struct Solver_Face_T face; ///< The base \ref Solver_Face_T.

	// Members required for 2nd order PDE terms.

	/// \brief Container for information relating to volumes on either side of the \ref DG_Solver_Face_T.
	struct Neigh_Info_DG {
		struct Multiarray_T* grad_coef_f; ///< The face contributions to the solution gradient coefficients.

		/** Linearization of \ref DG_Solver_Face_T::grad_coef_f wrt \ref Solver_Volume_T::sol_coef from the 'L'eft
		 *  (index 0) and 'R'ight (index 1) for all dimensions.
		 *
		 *  \note For non-boundary faces, this term is independent of the equation and variable under consideration
		 *        and a single entry is stored for all combinations. However, for boundary faces, as the boundary
		 *        values of the ghost state generally depend on all solution variables, all entries are required.
		 */
		const struct const_Matrix_R* d_g_coef_f__d_s_coef[2][DIM];
	} neigh_info[2]; ///< \ref Neigh_Info_DG. Uses the same indexing convention as that of \ref Face::neigh_info.
};

/// \brief Constructor for a derived \ref DG_Solver_Face_T.
void constructor_derived_DG_Solver_Face_T
	(struct Face* face_ptr,       ///< Pointer to the face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref DG_Solver_Face_T.
void destructor_derived_DG_Solver_Face_T
	(struct Face* face_ptr ///< Pointer to the face.
	);
