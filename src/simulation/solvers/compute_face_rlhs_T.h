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
 *  \brief Provides functions used for computing the face contributions to the right and left-hand side (rlhs) terms
 *         of supported schemes.
 */

struct const_Vector_R;
struct Matrix_R;
struct Matrix_T;
struct const_Matrix_T;
struct Multiarray_T;
struct Solver_Face_T;
struct Numerical_Flux_Input_T;
struct Numerical_Flux_T;
struct Simulation;
struct Solver_Storage_Implicit;

/** \brief Pointer to the function used to evaluate the rhs (and optionally lhs) face terms.
 *
 *  \param num_flux \ref Numerical_Flux_T.
 *  \param s_face   \ref Solver_Face_T.
 *  \param ssi      \ref Solver_Storage_Implicit.
 */
typedef void (*compute_rlhs_f_fptr_T)
	(const struct Numerical_Flux_T*const num_flux,
	 struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::tw0_vt_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__tw0_vt_fc_T
	(const int side_index,              ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T* s_face ///< The current \ref Face.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vs_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vs_fc_T
	(const int side_index,              ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T* s_face ///< The current \ref Face.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vr_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vr_fc_T
	(const int side_index,              ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T* s_face ///< The current \ref Face.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_ff_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_ff_fc_T
	(const int side_index,                   ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const s_face ///< The current face.
	);

/** \brief Permute the input matrix such that its ordering is such that it is in the reference coordinates of the
 *         face cubature nodes of the opposite volume. */
void permute_Matrix_T_fc
	(struct Matrix_T* data,             ///< The data to be permuted.
	 const char perm_layout,            ///< The layout in which to permute.
	 const int side_index_dest,         ///< The side index of the destination.
	 const struct Solver_Face_T* s_face ///< \ref Solver_Face_T.
	);
#if TYPE_RC == TYPE_COMPLEX
/// \brief Version of \ref permute_Matrix_T_fc taking a real input matrix.
void permute_Matrix_R_fc
	(struct Matrix_R* data,             ///< See brief.
	 const char perm_layout,            ///< See brief.
	 const int side_index_dest,         ///< See brief.
	 const struct Solver_Face_T* s_face ///< See brief.
	);
#endif
/** \brief Get the pointer to the appropriate \ref Solver_Element::nc_fc \ref const_Vector_T\*.
 *  \return See brief. */
const struct const_Vector_i* get_operator__nc_fc_T
	(const int side_index_dest,         ///< Defined for \ref permute_Multiarray_T_fc.
	 const struct Solver_Face_T* s_face ///< Defined for \ref permute_Multiarray_T_fc.
	);

/** \brief Construct the data members of the \ref Numerical_Flux_Input_T container which are specific to the face under
 *         consideration. */
void constructor_Numerical_Flux_Input_data_T
	(struct Numerical_Flux_Input_T* num_flux_i, ///< \ref Numerical_Flux_Input_T.
	 const struct Solver_Face_T* s_face,        ///< \ref Solver_Face_T.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

/// \brief Destructor for the data members of the \ref Numerical_Flux_Input_T container.
void destructor_Numerical_Flux_Input_data_T
	(struct Numerical_Flux_Input_T* num_flux_i ///< \ref Numerical_Flux_Input_T.
	);

/** \brief Constructor for the lhs face term of 1st order equations only.
 *  \return See brief. */
struct Matrix_T* constructor_lhs_f_1_T
	(const int side_index[2],                 /**< The indices of the affectee, affector, respectively. See the
	                                           *   comments in \ref compute_face_rlhs_dg.h for the convention. */
	 const struct Numerical_Flux_T* num_flux, ///< Defined for \ref compute_rlhs_f_fptr_T.
	 const struct Solver_Face_T* s_face       ///< Defined for \ref compute_rlhs_f_fptr_T.
	);

/** \brief Constructor for the partial lhs face term of 2nd order equations only (i.e. flux having dependence only on
 *         \ref Solver_Volume_T::grad_coef).
 *  \return See brief.
 *
 *  The lhs contribution computed here is termed 'p'artial as it must subsequently be multiplied by the
 *  d_g_coef__d_s_coef terms in order to obtain the final lhs terms.
 */
struct Matrix_T* constructor_lhs_p_f_2_T
	(const int side_index[2],                 ///< Defined for \ref constructor_lhs_f_1_T.
	 const struct Numerical_Flux_T* num_flux, ///< Defined for \ref compute_rlhs_f_fptr_T.
	 const struct Solver_Face_T* s_face       ///< Defined for \ref compute_rlhs_f_fptr_T.
	);

/** \brief Combine the input face cubature weights and normal flux values and add to the corresponding \ref
 *         Solver_Volume_T::flux_imbalance. */
void add_to_flux_imbalance_face_nf_w_T
	(const struct const_Matrix_T*const nf_fc, ///< The normal flux evaluated at the face cubature nodes.
	 const struct const_Vector_R*const w_fc,  ///< The face cubature weights.
	 const struct Solver_Face_T*const s_face  ///< The current face.
	);

/** \brief Version of \ref compute_rlhs_f_fptr_T computing only the rhs term.
 *
 *  The function is "dg-like" in the sense that it computes the rhs term as for the DG scheme, but can be used with test
 *  space not equal to trial space (Petrov-Galerkin).
 */
void compute_rhs_f_dg_like_T
	(const struct Numerical_Flux_T*const num_flux, ///< See brief.
	 struct Solver_Face_T*const s_face,            ///< See brief.
	 struct Solver_Storage_Implicit*const ssi      ///< See brief.
	);

/** \brief Permute the input multiarray such that its ordering is such that it is in the reference coordinates of the
 *         face cubature nodes of the opposite volume. */
void permute_Multiarray_T_fc
	(struct Multiarray_T* data,              ///< The data to be permuted.
	 const char perm_layout,                 ///< Defined for \ref permute_Multiarray_T_V.
	 const int side_index_dest,              ///< The side index of the destination.
	 const struct Solver_Face_T*const s_face ///< \ref Solver_Face_T.
	);
