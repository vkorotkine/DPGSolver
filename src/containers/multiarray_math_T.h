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
 *  \brief Provides Multiarray_\* math functions.
 */

#include <stdbool.h>
#include <stddef.h>

struct Multiarray_T;
struct const_Vector_i;
struct const_Vector_R;
struct const_Vector_T;
struct const_Matrix_R;
struct const_Matrix_T;
struct const_Multiarray_T;
struct Multiarray_Matrix_T;

/// \brief Transpose the \ref Multiarray_T\*'s memory by interpreting as a matrix with `ext_0 = extents[0]`.
void transpose_Multiarray_T
	(struct Multiarray_T* a, ///< Multiarray to be transposed.
	 const bool mem_only     ///< Flag for whether only the memory should be transposed (with extents unchanged).
	);

/// \brief Scale the \ref Multiarray_T\* by a constant value.
void scale_Multiarray_T
	(struct Multiarray_T* a, ///< The multiarray.
	 const Type val        ///< The value by which to scale.
	);

/// \brief Normalize the primary extents of a \ref Multiarray_T\*, optionally storing the value of the norms.
void normalize_Multiarray_T
	(struct Multiarray_T* a,      ///< Multiarray to be normalized.
	 const char*const norm_type,  ///< Defined for \ref norm_T.
	 const bool store_norms,      ///< Flag for whether norm values should be stored.
	 struct Multiarray_T* a_norms ///< Multiarray in which to store the norms if enabled.
	);

/// \brief Permute the the input multiarray according to the input permutation by reinterpretting as a matrix.
void permute_Multiarray_T
	(struct Multiarray_T* a, ///< Multiarray to be permuted.
	 const ptrdiff_t* p,     ///< Permutation indices.
	 const char perm_layout  ///< The layout in which to permute. Options: 'R'ow, 'C'olumn.
	);

/// \brief Call \ref permute_Multiarray_T using the indices of the vector as the permutation indices.
void permute_Multiarray_T_V
	(struct Multiarray_T* a,           ///< Defined for \ref permute_Multiarray_T.
	 const struct const_Vector_i* p_V, ///< Vector of permutation indices.
	 const char perm_layout            ///< Defined for \ref permute_Multiarray_T.
	);

/// \brief Call \ref scale_Matrix_T_by_Vector_R after interpreting the multiarray as a matrix.
void scale_Multiarray_T_by_Vector_R
	(const char side,                     ///< See brief.
	 const Real alpha,                    ///< See brief.
	 struct Multiarray_T*const a,         ///< See brief.
	 const struct const_Vector_R*const b, ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

/// \brief Sets `a = a + alpha*b`.
void add_in_place_Multiarray_T
	(const Real alpha,                  ///< Scaling constant.
	 struct Multiarray_T*const a,       ///< Multiarray to be modified.
	 const struct const_Multiarray_T* b ///< Multiarray to add to array to be modified.
	);

/// \brief Subtract the 2nd from the 1st multiarray in-place.
/// \deprecated Replace with calls to \ref add_in_place_Multiarray_T.
void subtract_in_place_Multiarray_T
	(struct Multiarray_T* a,            ///< 1st multiarray.
	 const struct const_Multiarray_T* b ///< 2nd multiarray.
	);

/** \brief Perform a matrix-matrix multiplication on a \ref const_Multiarray_T\*, interpreting the input multiarray as a
 *         a matrix with the appropriate extents, asserting that the input multiarrays have column-major layout.
 *
 *  The first extent **must** be equal to `ext_1` of the `a` matrix.
 *
 *  Defaults:
 *	- `trans_a_i = 'N'`;
 *	- `trans_b_i = 'N'`;
 *	- `layout = 'C'`.
 */
void mm_NNC_Multiarray_T
	(const Real alpha,                        ///< Defined for \ref mm_T.
	 const Real beta,                         ///< Defined for \ref mm_T.
	 const struct const_Matrix_R*const a,     ///< Defined for \ref mm_T.
	 const struct const_Multiarray_T*const b, ///< Input `b`.
	 struct Multiarray_T*const c              ///< Output `c`.
	);

/** \brief Perform a matrix-matrix multiplication on a \ref const_Multiarray_T\*, interpreting the input multiarray as a
 *         a matrix with the appropriate extents, asserting that the input multiarrays have column-major layout.
 *
 *  The first extent **must** be equal to `ext_1` of the `a` matrix.
 *  See comments in constructor_mm_NN1C_Matrix_T for the preset matrix-matrix multiplication parameters.
 */
void mm_NN1C_Multiarray_T
	(const struct const_Matrix_R*const a,     ///< Defined for \ref mm_T.
	 const struct const_Multiarray_T*const b, ///< Input `b`.
	 struct Multiarray_T*const c              ///< Output `c`.
	);

/// \brief Similar to \ref mm_NN1C_Multiarray_T but where the input `b` is overwritten by the result.
void mm_NN1C_overwrite_Multiarray_T
	(const struct const_Matrix_R*const a, ///< Defined for \ref mm_T.
	 struct Multiarray_T** b              ///< Input/Output `b`.
	);

/// \brief Reinterpret the \ref const_Multiarray_T\* as a \ref const_Matrix_T\* having the given input extents.
void reinterpret_const_Multiarray_as_Matrix_T
	(const struct const_Multiarray_T* a, ///< The multiarray.
	 const struct const_Matrix_T* a_M,   ///< The matrix.
	 const ptrdiff_t ext_0,              ///< The value of `ext_0` for the matrix.
	 const ptrdiff_t ext_1               ///< The value of `ext_1` for the matrix.
	);

/// \brief Reinterpret the \ref const_Matrix_T\* as a \ref const_Multiarray_T\* having the given input order/extents.
void reinterpret_const_Matrix_as_Multiarray_T
	(const struct const_Matrix_T* a_M,   ///< The matrix.
	 const struct const_Multiarray_T* a, ///< The multiarray.
	 const int order,                    ///< The value of `order` for the multiarray.
	 const ptrdiff_t* extents            ///< The value of `extents` for the multiarray.
	);

/** \brief Update the layout of each of the \ref Matrix_T\*s stored in the \ref Multiarray_Matrix_T\* to be that which
 *         is specified. */
void update_layout_Multiarray_Matrix_T
	(struct Multiarray_Matrix_T* a, ///< The input multiarray.
	 const char layout_o            ///< The desired output layout.
	);
