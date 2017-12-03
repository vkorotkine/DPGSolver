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

#ifndef DPG__multiarray_math_h__INCLUDED
#define DPG__multiarray_math_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* math functions.
 */

#include <stdbool.h>
#include <stddef.h>

struct Multiarray_d;
struct const_Vector_i;
struct const_Vector_d;
struct const_Matrix_d;
struct const_Multiarray_d;

/// \brief Transpose the \ref Multiarray_d\*'s memory by interpreting as a matrix with `ext_0 = extents[0]`.
void transpose_Multiarray_d
	(struct Multiarray_d* a, ///< Multiarray to be transposed.
	 const bool mem_only     ///< Flag for whether only the memory should be transposed (with extents unchanged).
	);

/// \brief Scale the \ref Multiarray_d\* by a constant value.
void scale_Multiarray_d
	(struct Multiarray_d* a, ///< The multiarray.
	 const double val        ///< The value by which to scale.
	);

/// \brief Normalize the primary extents of a \ref Multiarray_d\*, optionally storing the value of the norms.
void normalize_Multiarray_d
	(struct Multiarray_d* a,      ///< Multiarray to be normalized.
	 const char*const norm_type,  ///< Defined for \ref norm_d.
	 const bool store_norms,      ///< Flag for whether norm values should be stored.
	 struct Multiarray_d* a_norms ///< Multiarray in which to store the norms if enabled.
	);

/// \brief Permute the the input multiarray according to the input permutation by reinterpretting as a matrix.
void permute_Multiarray_d
	(struct Multiarray_d* a, ///< Multiarray to be permuted.
	 const ptrdiff_t* p,     ///< Permutation indices.
	 const char perm_layout  ///< The layout in which to permute. Options: 'R'ow, 'C'olumn.
	);

/// \brief Call \ref permute_Multiarray_d using the indices of the vector as the permutation indices.
void permute_Multiarray_d_V
	(struct Multiarray_d* a,           ///< Defined for \ref permute_Multiarray_d.
	 const struct const_Vector_i* p_V, ///< Vector of permutation indices.
	 const char perm_layout            ///< Defined for \ref permute_Multiarray_d.
	);

/// \brief Call \ref scale_Matrix_by_Vector_d after interpreting the multiarray as a matrix.
void scale_Multiarray_by_Vector_d
	(const char side,                     ///< Defined for \ref scale_Matrix_by_Vector_d.
	 const double alpha,                  ///< Defined for \ref scale_Matrix_by_Vector_d.
	 struct Multiarray_d*const a,         ///< Defined for \ref scale_Matrix_by_Vector_d.
	 const struct const_Vector_d*const b, ///< Defined for \ref scale_Matrix_by_Vector_d.
	 const bool invert_diag               ///< Defined for \ref scale_Matrix_by_Vector_d.
	);

/// \brief Subtract the 2nd from the 1st multiarray in-place.
void subtract_in_place_Multiarray_d
	(struct Multiarray_d* a,            ///< 1st multiarray.
	 const struct const_Multiarray_d* b ///< 2nd multiarray.
	);

/** \brief Perform a matrix-matrix multiplication on a \ref const_Multiarray_d\*, interpreting the input multiarray as a
 *         a matrix with the appropriate extents, asserting that the input multiarrays have column-major layout.
 *
 *  The first extent **must** be equal to `ext_1` of the `a` matrix.
 *
 *  Defaults:
 *	- `trans_a_i = 'N'`;
 *	- `trans_b_i = 'N'`;
 *	- `layout = 'C'`.
 */
void mm_NNC_Multiarray_d
	(const double alpha,                      ///< Defined for \ref mm_d.
	 const double beta,                       ///< Defined for \ref mm_d.
	 const struct const_Matrix_d*const a,     ///< Defined for \ref mm_d.
	 const struct const_Multiarray_d*const b, ///< Input `b`.
	 struct Multiarray_d*const c              ///< Output `c`.
	);

/** \brief Perform a matrix-matrix multiplication on a \ref const_Multiarray_d\*, interpreting the input multiarray as a
 *         a matrix with the appropriate extents, asserting that the input multiarrays have column-major layout.
 *
 *  The first extent **must** be equal to `ext_1` of the `a` matrix.
 *  See comments in \ref constructor_mm_NN1C_Matrix_d for the preset matrix-matrix multiplication parameters.
 */
void mm_NN1C_Multiarray_d
	(const struct const_Matrix_d*const a,     ///< Defined for \ref mm_d.
	 const struct const_Multiarray_d*const b, ///< Input `b`.
	 struct Multiarray_d*const c              ///< Output `c`.
	);

/** \brief Similar to \ref mm_NN1C_Multiarray_d but where the input `b` is overwritten by the result.
 *  See comments in \ref constructor_mm_NN1C_Matrix_d for the preset matrix-matrix multiplication parameters. */
void mm_NN1C_overwrite_Multiarray_d
	(const struct const_Matrix_d*const a, ///< Defined for \ref mm_d.
	 struct Multiarray_d** b              ///< Input/Output `b`.
	);

/** \brief Compute the extents of the output multiarray from a matrix-multiarray multiplication.
 *  \return Dynamically allocated extents. */
ptrdiff_t* compute_extents_mm_MMa
	(const ptrdiff_t ext_0,     ///< The value of `extents[0]`.
	 const int order,           ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t* extents_i ///< The input extents. Used to set all but the first entry.
	);

/// \brief Reinterpret the \ref const_Multiarray_d\* as a \ref const_Matrix_T\* having the given input extents.
void reinterpret_const_Multiarray_as_Matrix_d
	(const struct const_Multiarray_d* a, ///< The multiarray.
	 const struct const_Matrix_d* a_M,   ///< The matrix.
	 const ptrdiff_t ext_0,              ///< The value of `ext_0` for the matrix.
	 const ptrdiff_t ext_1               ///< The value of `ext_1` for the matrix.
	);

/// \brief Reinterpret the \ref const_Matrix_T\* as a \ref const_Multiarray_d\* having the given input order/extents.
void reinterpret_const_Matrix_as_Multiarray_d
	(const struct const_Matrix_d* a_M,   ///< The matrix.
	 const struct const_Multiarray_d* a, ///< The multiarray.
	 const int order,                    ///< The value of `order` for the multiarray.
	 const ptrdiff_t* extents            ///< The value of `extents` for the multiarray.
	);

#endif // DPG__multiarray_math_h__INCLUDED
