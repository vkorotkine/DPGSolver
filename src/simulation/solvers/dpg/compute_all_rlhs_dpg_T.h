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
 *  \brief Provides templated functions used for all contributions to the right and left-hand side (rlhs) terms of
 *         the DPG scheme.
 */

#include <stddef.h>
#include <stdbool.h>

struct Simulation;
struct Solver_Storage_Implicit;
struct Intrusive_List;
struct DPG_Solver_Volume_T;
struct DPG_Solver_Face_T;
struct Solver_Volume_T;
struct Vector_T;
struct Matrix_T;

/// \brief Compute all contributions to the rhs and lhs terms for the DPG scheme.
void compute_all_rlhs_dpg_T
	(const struct Simulation* sim,        ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi, ///< \ref Solver_Storage_Implicit.
	 struct Intrusive_List* volumes       ///< The list of volumes.
	);

/** \brief Get the appropriate sub-range of the \ref DPG_Solver_Element::cvt1_vt_vc operators.
 *  \return See brief. */
struct Multiarray_Operator get_operator__cvt1_vt_vc__rlhs_T
	(const struct DPG_Solver_Volume_T* dpg_s_vol ///< The current volume.
	);

/** \brief Construct the lhs left operator for the internal face dof for the dpg scheme.
 *  \return See brief. */
const struct const_Matrix_R* constructor_lhs_l_internal_face_dpg_T
	(const struct DPG_Solver_Volume_T* dpg_s_vol, ///< Pointer to the current volume.
	 const struct DPG_Solver_Face_T* dpg_s_face   ///< Pointer to the current face.
	);

/** \brief Return the number of degrees of freedom for all \ref Solver_Face::nf_coef adjacent to the current volume.
 *  \return See brief. */
ptrdiff_t compute_n_dof_nf_T
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/** \brief Constructor for the \ref Vector_T\* of indices of the global matrix in which to insert values for the current
 *         \ref Solver_Volume.
 *  \return See brief. */
/// \todo see if this can be made static.
const struct const_Vector_i* constructor_petsc_idxm_dpg_T
	(const ptrdiff_t n_dof,            ///< The number of local degrees of freedom.
	 const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/** \brief Add the face contributions to the rhs and lhs for 1st order equations.
 *  \note Columns are added to the lhs matrix for each of the face coefficient degrees of freedom.
 */
void add_to_rlhs__face_T
	(struct Vector_T* rhs,                        ///< Holds the values of the rhs.
	 struct Matrix_T** lhs_ptr,                   ///< Pointer to the matrix holding the values of the lhs.
	 const struct DPG_Solver_Volume_T* dpg_s_vol, ///< \ref DPG_Solver_Volume_T.
	 const struct Simulation* sim,                ///< \ref Simulation.
	 const bool include_internal                  /**< Flag for whether the internal face contributions should be
	                                               *   included. */
	);
