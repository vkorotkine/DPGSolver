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

#include <stddef.h>

struct Simulation;
struct Solver_Storage_Implicit;
struct DPG_Solver_Volume;
struct DPG_Solver_Face;
struct Solver_Volume;

/// \brief Compute all contributions to the rhs and lhs terms for the DPG scheme.
void compute_all_rlhs_dpg
	(const struct Simulation* sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi ///< \ref Solver_Storage_Implicit.
	);

/** \brief Get the appropriate sub-range of the \ref DPG_Solver_Element::cvt1_vt_vc operators.
 *  \return See brief. */
struct Multiarray_Operator get_operator__cvt1_vt_vc__rlhs
	(const struct DPG_Solver_Volume* dpg_s_vol ///< The current volume.
	);

/** \brief Construct the lhs left operator for the internal face dof for the dpg scheme.
 *  \return See brief. */
const struct const_Matrix_d* constructor_lhs_l_internal_face_dpg
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< Pointer to the current volume.
	 const struct DPG_Solver_Face* dpg_s_face   ///< Pointer to the current face.
	);

/** \brief Return the number of degrees of freedom for the \ref Solver_Face::nf_coef adjacent to the current volume.
 *  \return See brief. */
ptrdiff_t compute_n_dof_nf
	(const struct Solver_Volume* s_vol ///< The current volume.
	);

/** \brief Constructor for the \ref Vector_i\* of indices of the global matrix in which to insert values for the current
 *         \ref Solver_Volume.
 *  \return See brief. */
const struct const_Vector_i* constructor_petsc_idxm_dpg
	(const ptrdiff_t n_dof,            ///< The number of local degrees of freedom.
	 const struct Solver_Volume* s_vol ///< The current volume.
	);

#endif // DPG__compute_all_rlhs_dpg_h__INCLUDED
