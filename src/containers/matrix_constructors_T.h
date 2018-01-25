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
 *  \brief Provides Matrix_\* constructors and destructors.
 *
 *  Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

struct Matrix_R;
struct Multiarray_T;
struct Multiarray_Matrix_T;
struct const_Vector_R;
struct const_Vector_T;
struct const_Vector_i;
struct const_Matrix_T;
struct const_Matrix_R;
struct const_Multiarray_T;
struct const_Multiarray_Matrix_T;

// Default constructors ********************************************************************************************* //

/** \brief Constructor for a default \ref Matrix_T\*.
 *  \return Standard. */
struct Matrix_T* constructor_default_Matrix_T ();

/** \brief Constructor for a default \ref const_Matrix_T\*.
 *  \return Standard. */
const struct const_Matrix_T* constructor_default_const_Matrix_T ();

// Empty constructors *********************************************************************************************** //

/** \brief Constructs an empty \ref Matrix_T\*.
 *  \return Standard. */
struct Matrix_T* constructor_empty_Matrix_T
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

/** \brief `const` version of \ref constructor_empty_Matrix_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_empty_const_Matrix_T
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

// Zero constructors ************************************************************************************************ //

/** \brief Same as \ref constructor_empty_Matrix_T but with data calloc'ed.
 *  \return Standard. */
struct Matrix_T* constructor_zero_Matrix_T
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Copy constructor for a \ref Matrix_T\* from a \ref Matrix_T\*.
 *  \return Standard. */
struct Matrix_T* constructor_copy_Matrix_T
	(const struct Matrix_T* src ///< The source matrix.
	);

/** \brief Copy constructor for a \ref const_Matrix_T\* from a \ref const_Matrix_T\*.
 *  \return Standard. */
const struct const_Matrix_T* constructor_copy_const_Matrix_T
	(const struct const_Matrix_T* src ///< The source matrix.
	);

/** \brief Copy constructor for a \ref Matrix_T\* from a `const Type*`.
 *  \return Standard. */
struct Matrix_T* constructor_copy_Matrix_T_T
	(const char layout,          ///< Standard.
	 const ptrdiff_t ext_0,      ///< Standard.
	 const ptrdiff_t ext_1,      ///< Standard.
	 const Type*const data_src ///< The source data.
	);

/** \brief Copy constructor for a \ref Matrix_T\* from a \ref Matrix_T\*.
 *  \return See brief. */
struct Matrix_T* constructor_copy_Matrix_T_Matrix_R
	(struct Matrix_R* src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_copy_Matrix_T_Matrix_R.
 *  \return See brief. */
const struct const_Matrix_T* constructor_copy_const_Matrix_T_Matrix_R
	(const struct const_Matrix_R* src ///< See brief.
	);

/** \brief `const` version of \ref constructor_copy_Matrix_T_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_copy_const_Matrix_T_T
	(const char layout,          ///< Defined for \ref constructor_copy_Matrix_T_T.
	 const ptrdiff_t ext_0,      ///< Defined for \ref constructor_copy_Matrix_T_T.
	 const ptrdiff_t ext_1,      ///< Defined for \ref constructor_copy_Matrix_T_T.
	 const Type*const data_src ///< Defined for \ref constructor_copy_Matrix_T_T.
	);

/** \brief Copy constructor for a \ref const_Matrix_T\* from a partial number of rows/columns of another.
 *  \return Standard. */
const struct const_Matrix_T* constructor_copy_extract_const_Matrix_T
	(const struct const_Matrix_T*const src,    ///< The source Matrix.
	 const struct const_Vector_i*const indices ///< The indices of the rows/columns to copy.
	);

/** \brief Copy constructor for a `const` \ref const_Matrix_T\* from a `const` \ref const_Matrix_T\*.
 *  \return Standard. */
void const_constructor_copy_Matrix_T
	(const struct const_Matrix_T*const* dest, ///< Destination.
	 const struct const_Matrix_T*const src    ///< Source.
	);

// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a \ref Matrix_T\* from a `Type*`.
 *  \return Standard. */
struct Matrix_T* constructor_move_Matrix_T_T
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 Type*const data      ///< Standard.
	);

/** \brief `const` version of constructor_move_Matrix_T_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_move_const_Matrix_T_T
	(const char layout,      ///< Defined for \ref constructor_move_Matrix_T_T.
	 const ptrdiff_t ext_0,  ///< Defined for \ref constructor_move_Matrix_T_T.
	 const ptrdiff_t ext_1,  ///< Defined for \ref constructor_move_Matrix_T_T.
	 const bool owns_data,   ///< Defined for \ref constructor_move_Matrix_T_T.
	 const Type*const data   ///< Defined for \ref constructor_move_Matrix_T_T.
	);

/// \brief Move Constructor for a `const` \ref const_Matrix_T `*const` from a \ref Matrix_T\*.
void const_constructor_move_Matrix_T
	(const struct const_Matrix_T*const* dest, ///< Destination.
	 struct Matrix_T* src                     ///< Source.
	);

/// \brief Move Constructor for a `const` \ref const_Matrix_T `*const` from a \ref const_Matrix_T\*.
void const_constructor_move_const_Matrix_T
	(const struct const_Matrix_T*const* dest, ///< Destination.
	 const struct const_Matrix_T* src         ///< Source.
	);

// Special constructors (only available for real/complex types) ***************************************************** //
#ifdef TYPE_RC

/** \brief Copy constructor for a \ref Matrix_T\* with permutation.
 *  \return See brief. */
struct Matrix_T* constructor_copy_permute_Matrix_T
	(const struct Matrix_T*const src,       ///< The source matrix.
	 const struct const_Vector_i*const p_V, ///< Defined for \ref permute_Matrix_T_V.
	 const char perm_layout                 ///< The layout in which to permute. Options: 'R' (permute rows).
	);

/** \brief `const` version of \ref constructor_copy_permute_Matrix_T.
 *  \return See brief. */
const struct const_Matrix_T* constructor_copy_permute_const_Matrix_T
	(const struct const_Matrix_T*const src, ///< See brief.
	 const struct const_Vector_i*const p_V, ///< See brief.
	 const char perm_layout                 ///< See brief.
	);

/** \brief Constructor for a \ref Matrix_T\* which is a sub-block of the input matrix.
 *  \return See brief. */
struct Matrix_T* constructor_sub_block_Matrix_T
	(const ptrdiff_t row0,      ///< Index of the first row in the source matrix.
	 const ptrdiff_t col0,      ///< Index of the first col in the source matrix.
	 const ptrdiff_t n_row,     ///< Number of rows of the sub-matrix.
	 const ptrdiff_t n_col,     ///< Number of cols of the sub-matrix.
	 const struct Matrix_T* src ///< The source matrix.
	);

/** \brief Constructor for a \ref const_Matrix_T\* from a subset of the input matrix.
 *  \return Standard. */
const struct const_Matrix_T* constructor_subset_const_Matrix_T
	(const struct const_Matrix_T* src,       ///< The source matrix.
	 const struct const_Vector_i* ind_subset ///< The indices of the subset matrix.
	);

/** \brief Constructor for a \ref Matrix_T\* as a copy of the transpose of the input matrix.
 *  \return Standard. */
struct Matrix_T* constructor_copy_transpose_Matrix_T
	(struct Matrix_T* a, ///< The input matrix.
	 const bool mem_only ///< Defined for \ref transpose_Matrix_T.
	);

/** \brief Constructor for a block diagonal \ref const_Matrix_T\* with blocks set to the input matrix.
 *  \return Standard.
 *
 *  \note The input matrix need not be square.
 */
const struct const_Matrix_T* constructor_block_diagonal_const_Matrix_T
	(const struct const_Matrix_T* src_b, ///< The source block.
	 const ptrdiff_t n_blocks            ///< The number of blocks in the output.
	);

/** \brief Constructor for a diagonal \ref Matrix_T\* with entries set to the input value.
 *  \return Standard. */
struct Matrix_T* constructor_diagonal_Matrix_T_T
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< The dimensions of the square matrix.
	 const Type val       ///< The value.
	);

/** \brief Constructor for an identity \ref Matrix_T\*.
 *  \return Standard. */
struct Matrix_T* constructor_identity_Matrix_T
	(const char layout,    ///< Standard.
	 const ptrdiff_t ext_0 ///< The dimensions of the square matrix.
	);

/** \brief `const` version of \ref constructor_identity_Matrix_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_identity_const_Matrix_T
	(const char layout,    ///< Standard.
	 const ptrdiff_t ext_0 ///< The dimensions of the square matrix.
	);

/** \brief Constructor for the inverse of the input \ref Matrix_T\*.
 *  \return Standard. */
struct Matrix_T* constructor_inverse_Matrix_T
	(struct Matrix_T* src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_inverse_Matrix_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_inverse_const_Matrix_T
	(const struct const_Matrix_T* src ///< Defined for \ref constructor_inverse_Matrix_T.
	);

/** \brief Constructor for a \ref Matrix_T\* from the solution of a linear system using LAPACKE_dsgesv.
 *  \return Standard.
 *
 *  Reference: [LAPACKE_dsgesv][dsgesv].
 *
 *  <!-- References: -->
 *  [dsgesv]: https://software.intel.com/en-us/mkl-developer-reference-c-gesv
 */
struct Matrix_T* constructor_sgesv_Matrix_T
	(struct Matrix_T* A_i, ///< The LHS input matrix A.
	 struct Matrix_T* B_i  ///< The RHS input matrix B.
	);

/** \brief `const` version of constructor_sgesv_Matrix_T.
 *  \return See brief. */
const struct const_Matrix_T* constructor_sgesv_const_Matrix_T
	(const struct const_Matrix_T* A_i, ///< Defined for \ref constructor_sgesv_Matrix_T.
	 const struct const_Matrix_T* B_i  ///< Defined for \ref constructor_sgesv_Matrix_T.
	);

/** \brief Constructor for a \ref Matrix_T\* from the solution of a linear system using LAPACKE_dsysv.
 *  \return Standard.
 *
 *  Reference: [LAPACKE_dsysv][dsysv].
 *
 *  <!-- References: -->
 *  [dsysv]: https://software.intel.com/en-us/mkl-developer-reference-c-sysv
 */
struct Matrix_T* constructor_sysv_Matrix_T
	(struct Matrix_T* A_i, ///< The LHS input matrix A.
	 struct Matrix_T* B_i  ///< The RHS input matrix B.
	);

/** \brief `const` version of constructor_sysv_Matrix_T.
 *  \return See brief. */
const struct const_Matrix_T* constructor_sysv_const_Matrix_T
	(const struct const_Matrix_T* A_i, ///< Defined for \ref constructor_sysv_Matrix_T.
	 const struct const_Matrix_T* B_i  ///< Defined for \ref constructor_sysv_Matrix_T.
	);

/** \brief Constructor for a \ref Matrix_T\* from a matrix-matrix multiplication.
 *  \return Standard. */
struct Matrix_T* constructor_mm_Matrix_T
	(const char trans_a_i,                ///< Defined for \ref mm_T.
	 const char trans_b_i,                ///< Defined for \ref mm_T.
	 const Real alpha,                    ///< Defined for \ref mm_T.
	 const struct const_Matrix_T*const a, ///< Defined for \ref mm_T.
	 const struct const_Matrix_T*const b, ///< Defined for \ref mm_T.
	 const char layout                    ///< The `layout` of the constructed \ref Matrix_T.
	);

/** \brief `const` version of \ref constructor_mm_Matrix_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_mm_const_Matrix_T
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const Real alpha,                    ///< See brief.
	 const struct const_Matrix_T*const a, ///< See brief.
	 const struct const_Matrix_T*const b, ///< See brief.
	 const char layout                    ///< See brief.
	);

/** \brief Version of \ref constructor_mm_Matrix_T with 'R'eal and 'T'ype inputs.
 *  \return See brief. */
struct Matrix_T* constructor_mm_RT_Matrix_T
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const Real alpha,                    ///< See brief.
	 const struct const_Matrix_R*const a, ///< See brief.
	 const struct const_Matrix_T*const b, ///< See brief.
	 const char layout                    ///< See brief.
	);

/** \brief `const` version of \ref constructor_mm_RT_Matrix_T.
 *  \return See brief. */
const struct const_Matrix_T* constructor_mm_RT_const_Matrix_T
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const Real alpha,                    ///< See brief.
	 const struct const_Matrix_R*const a, ///< See brief.
	 const struct const_Matrix_T*const b, ///< See brief.
	 const char layout                    ///< See brief.
	);

/** \brief Constructor for a \ref Matrix_T\* from a matrix-matrix multiplication using default values.
 *
 *  Defaults:
 *	- `trans_a_i = 'N'`;
 *	- `trans_b_i = 'N'`;
 *	- `alpha = 1.0`;
 *	- `layout = 'R'`.
 *
 *  \return Standard. */
struct Matrix_T* constructor_mm_NN1R_Matrix_T
	(const struct const_Matrix_T*const a, ///< Defined for \ref mm_T.
	 const struct const_Matrix_T*const b  ///< Defined for \ref mm_T.
	);

/** \brief `const` version of \ref constructor_mm_NN1R_Matrix_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_mm_NN1R_const_Matrix_T
	(const struct const_Matrix_T*const a, ///< Defined for \ref constructor_mm_NN1R_Matrix_T.
	 const struct const_Matrix_T*const b  ///< Defined for \ref constructor_mm_NN1R_Matrix_T.
	);

/** \brief Constructor for a \ref Matrix_T\* from a matrix-matrix multiplication using default values.
 *
 *  Defaults:
 *	- `trans_a_i = 'N'`;
 *	- `trans_b_i = 'N'`;
 *	- `alpha = 1.0`;
 *	- `layout = 'C'`.
 *
 *  \return Standard. */
struct Matrix_T* constructor_mm_NN1C_Matrix_T
	(const struct const_Matrix_T*const a, ///< Defined for \ref mm_T.
	 const struct const_Matrix_T*const b  ///< Defined for \ref mm_T.
	);

/** \brief `const` version of \ref constructor_mm_NN1C_Matrix_T.
 *  \return Standard. */
const struct const_Matrix_T* constructor_mm_NN1C_const_Matrix_T
	(const struct const_Matrix_T*const a, ///< Defined for \ref constructor_mm_NN1C_Matrix_T.
	 const struct const_Matrix_T*const b  ///< Defined for \ref constructor_mm_NN1C_Matrix_T.
	);

/** \brief Constructor for a \ref Matrix_T\* from a matrix-diagonal matrix multiplication taking a \ref Vector_T input.
 *  \return Standard.
 *
 *  The diagonal matrix is input as a vector and may be applied either from the left or the right.
 */
struct Matrix_T* constructor_mm_diag_Matrix_T_R
	(const Real alpha,                    ///< Defined for \ref mm_T.
	 const struct const_Matrix_T*const a, ///< Input matrix to be multiplied by the diagonal.
	 const struct const_Vector_R*const b, ///< Vector storing the entries of the diagonal matrix.
	 const char side,                     ///< The side from which to apply the diagonal matrix.
	 const bool invert_diag               ///< Defined for \ref scale_Matrix_T_by_Vector_R.
	);

/** \brief `const` version of \ref constructor_mm_diag_Matrix_T_R.
 *  \return See brief. */
const struct const_Matrix_T* constructor_mm_diag_const_Matrix_T_R
	(const Real alpha,                    ///< See brief.
	 const struct const_Matrix_T*const a, ///< See brief.
	 const struct const_Vector_R*const b, ///< See brief.
	 const char side,                     ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

/** \brief 'T'ype-'T'ype version of \ref constructor_mm_diag_Matrix_T_R.
 *  \return Standard. */
struct Matrix_T* constructor_mm_diag_Matrix_T
	(const Real alpha,                    ///< See brief.
	 const struct const_Matrix_T*const a, ///< See brief.
	 const struct const_Vector_T*const b, ///< See brief.
	 const char side,                     ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

/** \brief `const` version of \ref constructor_mm_diag_Matrix_T.
 *  \return See brief. */
const struct const_Matrix_T* constructor_mm_diag_const_Matrix_T
	(const Real alpha,                    ///< See brief.
	 const struct const_Matrix_T*const a, ///< See brief.
	 const struct const_Vector_T*const b, ///< See brief.
	 const char side,                     ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

/// \brief Set a \ref Matrix_T\* from a sub range of a \ref Multiarray_T\*.
void set_Matrix_from_Multiarray_T
	(struct Matrix_T* dest,            ///< The destination.
	 struct Multiarray_T* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Matrix_from_Multiarray_T.
void set_const_Matrix_from_Multiarray_T
	(const struct const_Matrix_T* dest,    ///< Defined for \ref set_Matrix_from_Multiarray_T.
	 const struct const_Multiarray_T* src, ///< Defined for \ref set_Matrix_from_Multiarray_T.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref set_Matrix_from_Multiarray_T.
	);

/// \brief Set a \ref Matrix_T\* from an entry of a \ref Multiarray_Matrix_T\*.
void set_Matrix_from_Multiarray_Matrix_T
	(struct Matrix_T* dest,            ///< The destination.
	 struct Multiarray_Matrix_T* src,  ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Matrix_from_Multiarray_Matrix_T.
void set_const_Matrix_from_Multiarray_Matrix_T
	(const struct const_Matrix_T* dest,           ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_T.
	 const struct const_Multiarray_Matrix_T* src, ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_T.
	 const ptrdiff_t*const sub_indices            ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_T.
	);
#endif
// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Matrix_T\*.
void destructor_Matrix_T
	(struct Matrix_T* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_Matrix_T.
void destructor_const_Matrix_T
	(const struct const_Matrix_T* a ///< Standard.
	);

/// \brief Destructs a \ref Matrix_T\* if it is not NULL.
void destructor_conditional_Matrix_T
	(struct Matrix_T* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_conditional_Matrix_T.
void destructor_conditional_const_Matrix_T
	(const struct const_Matrix_T* a ///< Standard.
	);
