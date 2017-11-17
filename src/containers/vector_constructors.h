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

#ifndef DPG__vector_constructors_h__INCLUDED
#define DPG__vector_constructors_h__INCLUDED
/** \file
 *  \brief Provides Vector_\* constructors and destructors.
 */

#include <stddef.h>
#include <stdbool.h>

struct Matrix_d;
struct Multiarray_d;
struct const_Vector_d;
struct const_Matrix_d;
struct const_Matrix_i;
struct const_Multiarray_d;

// Default constructors ********************************************************************************************* //

/** \brief Constructs a default \ref Vector_d\*.
 *  \return Standard. */
struct Vector_d* constructor_default_Vector_d ();

/** \brief Constructor for a default \ref const_Vector_d\*.
 *  \return Standard. */
const struct const_Vector_d* constructor_default_const_Vector_d ();

/** \brief Constructs a default \ref Vector_i\*.
 *  \return Standard. */
struct Vector_i* constructor_default_Vector_i ();

/** \brief Constructs a default \ref Vector_i\*\*.
 *  \return Standard. */
struct Vector_i** constructor_default_Vector_i_2
	(const ptrdiff_t n_dest ///< The number of \ref Vector_i\* components.
	);

// Empty constructors *********************************************************************************************** //

/** \brief Constructs an empty \ref Vector_d\*.
 *  \return Standard. */
struct Vector_d* constructor_empty_Vector_d
	(const ptrdiff_t ext_0 ///< Defined in \ref Vector_d.
	);

/** \brief Constructs an empty \ref Vector_i\*.
 *  \return Standard. */
struct Vector_i* constructor_empty_Vector_i
	(const ptrdiff_t ext_0 ///< The value of ext_0.
	);

// Zero constructors ************************************************************************************************ //

/** \brief The same as \ref constructor_empty_Vector_i but allocating using calloc.
 *  \return Standard. */
struct Vector_i* constructor_zero_Vector_i
	(const ptrdiff_t ext_0 ///< Defined for \ref constructor_empty_Vector_i.
	);

/** \brief The same as \ref constructor_empty_Vector_d but allocating using calloc.
 *  \return Standard. */
struct Vector_d* constructor_zero_Vector_d
	(const ptrdiff_t ext_0 ///< Defined for \ref constructor_empty_Vector_d.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Copy constructor for a \ref Vector_i\* from a `Vector_i*`.
 *  \return Standard. */
struct Vector_i* constructor_copy_Vector_i
	(const struct Vector_i*const src ///< The source data.
	);

/** \brief Copy constructor for a \ref Vector_d\* from a `const double*`.
 *  \return Standard. */
struct Vector_d* constructor_copy_Vector_d_d
	(const ptrdiff_t ext_0,      ///< The value of ext_0.
	 const double*const data_src ///< The source data.
	);

/** \brief `const` version of \ref constructor_copy_Vector_d_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_copy_const_Vector_d_d
	(const ptrdiff_t ext_0,      ///< See brief.
	 const double*const data_src ///< See brief.
	);

/** \brief Copy constructor for a \ref Vector_i\* from a `const int*`.
 *  \return Standard. */
struct Vector_i* constructor_copy_Vector_i_i
	(const ptrdiff_t ext_0,   ///< The value of ext_0.
	 const int*const data_src ///< The source data.
	);

/** \brief `const` version of \ref constructor_copy_Vector_i_i.
 *  \return Standard. */
const struct const_Vector_i* constructor_copy_const_Vector_i_i
	(const ptrdiff_t ext_0,   ///< Defined for \ref constructor_copy_Vector_i_i.
	 const int*const data_src ///< Defined for \ref constructor_copy_Vector_i_i.
	);

// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a \ref Vector_i\* from a `int*`.
 *  \return Standard. */
struct Vector_i* constructor_move_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of ext_0.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

/** \brief Move constructor for a \ref Vector_d\* from a `double*`.
 *  \return Standard. */
struct Vector_d* constructor_move_Vector_d_d
	(const ptrdiff_t ext_0, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 double*const data      ///< Standard.
	);

/** \brief Move constructor for a \ref const_Vector_i\* from a `const int*`.
 *  \return Standard. */
struct const_Vector_i* constructor_move_const_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of ext_0.
	 const bool owns_data,  ///< Standard.
	 const int*const data   ///< Standard.
	);

/** \brief Move constructor for a \ref const_Vector_i\* from a row of a \ref const_Matrix_i\*.
 *  \return Standard. */
const struct const_Vector_i* constructor_move_const_Vector_Matrix_row_i
	(const int row,                    ///< The row index.
	 const struct const_Matrix_i* src, ///< The source matrix.
	 const int owns_data               ///< Standard.
	);

/** \brief Move constructor for a \ref Vector_d\* from a \ref Matrix_d\*.
 *  \return Standard. */
struct Vector_d* constructor_move_Vector_d_Matrix_d
	(struct Matrix_d* src ///< The source matrix.
	);

/// \brief Move constructor for a `const` \ref const_Vector_d `*const`.
void const_constructor_move_Vector_d
	(const struct const_Vector_d*const* dest, ///< Destination.
	 struct Vector_d* src                     ///< Source.
	);

/// \brief Move constructor for a `const` \ref const_Vector_i `*const`.
void const_constructor_move_Vector_i
	(const struct const_Vector_i*const* dest, ///< Destination.
	 struct Vector_i* src                     ///< Source.
	);

/// \brief Move constructor for a `const` \ref const_Vector_i `*const` from a `const` \ref const_Vector_i `*const`.
void const_constructor_move_const_Vector_i
	(const struct const_Vector_i*const* dest, ///< Destination.
	 const struct const_Vector_i* src         ///< Source.
	);

/// \brief Move constructor for a `const` \ref const_Vector_d `*const` from a `const` \ref const_Vector_d `*const`.
void const_constructor_move_const_Vector_d
	(const struct const_Vector_d*const* dest, ///< Destination.
	 const struct const_Vector_d* src         ///< Source.
	);

// Set constructors ************************************************************************************************* //

/** \brief Constructor for a \ref Vector_d\* from a sub range of a \ref Multiarray_d\*.
 *  \return Standard. */
struct Vector_d* constructor_set_Vector_d_Multiarray_d
	(struct Multiarray_d* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/** \brief `const` version of \ref constructor_set_Vector_d_Multiarray_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_set_const_Vector_d_Multiarray_d
	(const struct const_Multiarray_d* src, ///< Defined for \ref constructor_set_Vector_d_Multiarray_d.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref constructor_set_Vector_d_Multiarray_d.
	);

// Special constructors ********************************************************************************************* //

/** \brief Constructor for a \ref Vector_d\* as the inverse of the source vector.
 *  \return See brief. */
struct Vector_d* constructor_inverse_Vector_d
	(const struct const_Vector_d* src ///< The source vector.
	);

/** \brief `const` version of \ref constructor_inverse_Vector_d.
 *  \return See brief. */
const struct const_Vector_d* constructor_inverse_const_Vector_d
	(const struct const_Vector_d* src ///< See brief.
	);

/** \brief Constructor for a \ref const_Vector_d\* using an element-wise multiplication of each entry of the inputs.
 *  \return See brief. */
const struct const_Vector_d* constructor_dot_mult_const_Vector_d
	(const struct const_Vector_d* a, ///< The 1st input.
	 const struct const_Vector_d* b, ///< The 2nd input.
	 const int n_repeated            ///< The number of times the sub-vector should be repeated.
	);

/** \brief Constructs a \ref Vector_d\* as the sum of the rows/columns of the input Matrix in the specified direction.
 *  \return Standard. */
struct Vector_d* constructor_sum_Vector_d_const_Matrix_d
	(const char sum_dir,                   ///< The direction in which to sum the entries. Options: 'R'ow, 'C'olumn.
	 const struct const_Matrix_d*const src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_sum_Vector_d_const_Matrix_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_sum_const_Vector_d_const_Matrix_d
	(const char sum_dir,                   ///< Defined for \ref constructor_sum_Vector_d_const_Matrix_d.
	 const struct const_Matrix_d*const src ///< Defined for \ref constructor_sum_Vector_d_const_Matrix_d.
	);

/** \brief Constructor for a \ref Vector_d\* from a matrix-vector multiplication.
 *  \return Standard. */
struct Vector_d* constructor_mv_Vector_d
	(const char trans_a_i,                ///< Defined for \ref mv_d.
	 const double alpha,                  ///< Defined for \ref mv_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref mv_d.
	 const struct const_Vector_d*const b  ///< Defined for \ref mv_d.
	);

/** \brief `const` version of \ref constructor_mv_Vector_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_mv_const_Vector_d
	(const char trans_a_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const struct const_Matrix_d*const a, ///< See brief.
	 const struct const_Vector_d*const b  ///< See brief.
	);

/** \brief Constructor for a \ref Vector_d\* from the solution of a linear system using [LAPACKE_dsgesv][dsgesv].
 *  \return Standard.
 *
 *  <!-- References: -->
 *  [dsgesv]: https://software.intel.com/en-us/mkl-developer-reference-c-gesv
 */
struct Vector_d* constructor_sgesv_Vector_d
	(struct Matrix_d* A_i, ///< The LHS input matrix A.
	 struct Vector_d* B_i  ///< The RHS input vector B.
	);

/** \brief `const` version of \ref constructor_sgesv_Vector_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_sgesv_const_Vector_d
	(const struct const_Matrix_d* A_i, ///< The LHS input matrix A.
	 const struct const_Vector_d* B_i  ///< The RHS input vector B.
	);

/// \brief Set a \ref Vector_d\* from a sub range of a \ref Matrix_d\*.
void set_Vector_from_Matrix_d
	(struct Vector_d* dest,            ///< The destination.
	 struct Matrix_d* src,             ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Vector_from_Matrix_d.
void set_const_Vector_from_Matrix_d
	(const struct const_Vector_d* dest, ///< Defined for \ref set_Vector_from_Matrix_d.
	 const struct const_Matrix_d* src,  ///< Defined for \ref set_Vector_from_Matrix_d.
	 const ptrdiff_t*const sub_indices  ///< Defined for \ref set_Vector_from_Matrix_d.
	);

/// \brief Set a \ref Vector_d\* from a sub range of a \ref Multiarray_d\*.
void set_Vector_from_Multiarray_d
	(struct Vector_d* dest,            ///< The destination.
	 struct Multiarray_d* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Vector_from_Multiarray_d.
void set_const_Vector_from_Multiarray_d
	(const struct const_Vector_d* dest,    ///< Defined for \ref set_Vector_from_Multiarray_d.
	 const struct const_Multiarray_d* src, ///< Defined for \ref set_Vector_from_Multiarray_d.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref set_Vector_from_Multiarray_d.
	);

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Vector_d\*.
void destructor_Vector_d
	(struct Vector_d* a ///< Standard.
	);

/// \brief Destructs a \ref const_Vector_d\*.
void destructor_const_Vector_d
	(const struct const_Vector_d* a ///< Standard.
	);

/// \brief Destructs a \ref Vector_i\*.
void destructor_Vector_i
	(struct Vector_i* a ///< Standard.
	);

/// \brief Destructs a \ref const_Vector_i\*.
void destructor_const_Vector_i
	(const struct const_Vector_i* a ///< Standard.
	);

#endif // DPG__vector_constructors_h__INCLUDED
