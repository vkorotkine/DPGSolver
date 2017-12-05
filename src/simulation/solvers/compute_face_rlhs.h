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

#ifndef DPG__compute_face_rlhs_h__INCLUDED
#define DPG__compute_face_rlhs_h__INCLUDED
/** \file
 *  \brief Provides functions used for computing the face contributions to the right and left-hand side (rlhs) terms
 *         of supported schemes.
 */

struct Matrix_d;
struct Solver_Face;
struct Numerical_Flux_Input;
struct Numerical_Flux;
struct Simulation;

/** \brief Get the pointer to the appropriate \ref Solver_Element::tw0_vt_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__tw0_vt_fc
	(const int side_index,            ///< The index of the side of the face under consideration.
	 const struct Solver_Face* s_face ///< The current \ref Face.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vs_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vs_fc
	(const int side_index,            ///< The index of the side of the face under consideration.
	 const struct Solver_Face* s_face ///< The current \ref Face.
	);

/** \brief Permute the input matrix such that its ordering is such that it is in the reference coordinates of the
 *         face cubature nodes of the opposite volume. */
void permute_Matrix_d_fc
	(struct Matrix_d* data,           ///< The data to be permuted.
	 const char perm_layout,          ///< The layout in which to permute.
	 const int side_index_dest,       ///< The side index of the destination.
	 const struct Solver_Face* s_face ///< \ref Solver_Face.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::nc_fc \ref const_Vector_T\*.
 *  \return See brief. */
const struct const_Vector_i* get_operator__nc_fc
	(const int side_index_dest,       ///< Defined for \ref permute_Multiarray_d_fc.
	 const struct Solver_Face* s_face ///< Defined for \ref permute_Multiarray_d_fc.
	);

/** \brief Construct the data members of the \ref Numerical_Flux_Input container which are specific to the face under
 *         consideration. */
void constructor_Numerical_Flux_Input_data
	(struct Numerical_Flux_Input* num_flux_i, ///< \ref Numerical_Flux_Input.
	 const struct Solver_Face* s_face,        ///< \ref Solver_Face.
	 const struct Simulation* sim             ///< \ref Simulation.
	);

/// \brief Destructor for the data members of the \ref Numerical_Flux_Input container.
void destructor_Numerical_Flux_Input_data
	(struct Numerical_Flux_Input* num_flux_i ///< \ref Numerical_Flux_Input.
	);

/** \brief Constructor for the lhs face term of 1st order equations only.
 *  \return See brief. */
struct Matrix_d* constructor_lhs_f_1
	(const int side_index[2],               /**< The indices of the affectee, affector, respectively. See the
	                                         *   comments in \ref compute_face_rlhs_dg.h for the convention. */
	 const struct Numerical_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Solver_Face* s_face       ///< Defined for \ref compute_rlhs_fptr.
	);

#endif // DPG__compute_face_rlhs_h__INCLUDED
