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

#ifndef DPG__Matrix_h__INCLUDED
#define DPG__Matrix_h__INCLUDED
/** \file
 *  \brief Provides Matrix_\* constructors and destructors.
 *
 *  Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

struct Multiarray_d;
struct Multiarray_Matrix_d;
struct const_Vector_i;
struct const_Matrix_i;
struct const_Multiarray_d;
struct const_Multiarray_Matrix_d;

// Default constructors ********************************************************************************************* //

/** \brief Constructor for a default \ref Matrix_d\*.
 *  \return Standard. */
struct Matrix_d* constructor_default_Matrix_d ();

/** \brief Constructor for a default \ref const_Matrix_d\*.
 *  \return Standard. */
const struct const_Matrix_d* constructor_default_const_Matrix_d ();

// Empty constructors *********************************************************************************************** //

/** \brief Constructs an empty \ref Matrix_d\*.
 *  \return Standard. */
struct Matrix_d* constructor_empty_Matrix_d
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

/** \brief Constructs an empty \ref Matrix_i\*.
 *  \return Standard. */
struct Matrix_i* constructor_empty_Matrix_i
	(const char layout,  ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Copy constructor for a \ref Matrix_d\* from a \ref Matrix_d\*.
 *  \return Standard. */
struct Matrix_d* constructor_copy_Matrix_d
	(struct Matrix_d* src /// The source matrix.
	);

/** \brief Copy constructor for a \ref const_Matrix_d\* from a \ref const_Matrix_d\*.
 *  \return Standard. */
const struct const_Matrix_d* constructor_copy_const_Matrix_d
	(const struct const_Matrix_d* src /// The source matrix.
	);

/** \brief Copy constructor for a \ref Matrix_i\* from a `const int*`.
 *  \return Standard. */
struct Matrix_i* constructor_copy_Matrix_i_i
	(const char layout,       ///< Standard.
	 const ptrdiff_t ext_0,   ///< Standard.
	 const ptrdiff_t ext_1,   ///< Standard.
	 const int*const data_src ///< The source data.
	);

/** \brief Copy constructor for a \ref Matrix_d\* from a `const double*`.
 *  \return Standard. */
struct Matrix_d* constructor_copy_Matrix_d_d
	(const char layout,          ///< Standard.
	 const ptrdiff_t ext_0,      ///< Standard.
	 const ptrdiff_t ext_1,      ///< Standard.
	 const double*const data_src ///< The source data.
	);

/** \brief `const` version of \ref constructor_copy_Matrix_d_d.
 *  \return Standard. */
const struct const_Matrix_d* constructor_copy_const_Matrix_d_d
	(const char layout,          ///< Defined for \ref constructor_copy_Matrix_d_d.
	 const ptrdiff_t ext_0,      ///< Defined for \ref constructor_copy_Matrix_d_d.
	 const ptrdiff_t ext_1,      ///< Defined for \ref constructor_copy_Matrix_d_d.
	 const double*const data_src ///< Defined for \ref constructor_copy_Matrix_d_d.
	);

/** \brief Copy constructor for a \ref const_Matrix_d\* from a partial number of rows/columns of another.
 *  \return Standard. */
const struct const_Matrix_d* constructor_copy_extract_const_Matrix_d
	(const struct const_Matrix_d*const src,    ///< The source Matrix.
	 const struct const_Vector_i*const indices ///< The indices of the rows/columns to copy.
	);

/** \brief Copy constructor for a `const` \ref const_Matrix_d\* from a `const` \ref const_Matrix_d\*.
 *  \return Standard. */
void const_constructor_copy_Matrix_d
	(const struct const_Matrix_d*const* dest, ///< Destination.
	 const struct const_Matrix_d*const src    ///< Source.
	);

// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a \ref Matrix_d\* from a `double*`.
 *  \return Standard. */
struct Matrix_d* constructor_move_Matrix_d_d
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 double*const data      ///< Standard.
	);

/** \brief `const` version of constructor_move_Matrix_d_d.
 *  \return Standard. */
const struct const_Matrix_d* constructor_move_const_Matrix_d_d
	(const char layout,      ///< Defined for \ref constructor_move_Matrix_d_d.
	 const ptrdiff_t ext_0,  ///< Defined for \ref constructor_move_Matrix_d_d.
	 const ptrdiff_t ext_1,  ///< Defined for \ref constructor_move_Matrix_d_d.
	 const bool owns_data,   ///< Defined for \ref constructor_move_Matrix_d_d.
	 const double*const data ///< Defined for \ref constructor_move_Matrix_d_d.
	);

/** \brief Move constructor for a \ref Matrix_i\* from a `int*`.
 *  \return Standard. */
struct Matrix_i* constructor_move_Matrix_i_i
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

/// \brief Move Constructor for a `const` \ref const_Matrix_d `*const`.
void const_constructor_move_Matrix_d
	(const struct const_Matrix_d*const* dest, ///< Destination.
	 struct Matrix_d* src                     ///< Source.
	);

/** \brief Move Constructor for a `const` \ref const_Matrix_i `*const`.
 *  \return Standard. */
void const_constructor_move_Matrix_i
	(const struct const_Matrix_i*const* dest, ///< Destination.
	 struct Matrix_i* src                     ///< Source.
	);

// Special constructors ********************************************************************************************* //

/** \brief Constructor for a \ref const_Matrix_d\* from a subset of the input matrix.
 *  \return Standard. */
const struct const_Matrix_d* constructor_subset_const_Matrix_d
	(const struct const_Matrix_d* src,       ///< The source matrix.
	 const struct const_Vector_i* ind_subset ///< The indices of the subset matrix.
	);

/** \brief Constructor for a \ref Matrix_d\* as a copy of the transpose of the input matrix.
 *  \return Standard. */
struct Matrix_d* constructor_copy_transpose_Matrix_d
	(struct Matrix_d* a, ///< The input matrix.
	 const bool mem_only ///< Defined for \ref transpose_Matrix_d.
	);

/** \brief Constructor for a diagonal \ref Matrix_d\* with entries set to the input value.
 *  \return Standard. */
struct Matrix_d* constructor_diagonal_Matrix_d_d
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< The dimensions of the square matrix.
	 const double val       ///< The value.
	);

/** \brief Constructor for an identity \ref Matrix_d\*.
 *  \return Standard. */
struct Matrix_d* constructor_identity_Matrix_d
	(const char layout,    ///< Standard.
	 const ptrdiff_t ext_0 ///< The dimensions of the square matrix.
	);

/** \brief `const` version of \ref constructor_identity_Matrix_d.
 *  \return Standard. */
const struct const_Matrix_d* constructor_identity_const_Matrix_d
	(const char layout,    ///< Standard.
	 const ptrdiff_t ext_0 ///< The dimensions of the square matrix.
	);

/** \brief Constructor for the inverse of the input \ref Matrix_d\*.
 *  \return Standard. */
struct Matrix_d* constructor_inverse_Matrix_d
	(struct Matrix_d* src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_inverse_Matrix_d.
 *  \return Standard. */
const struct const_Matrix_d* constructor_inverse_const_Matrix_d
	(const struct const_Matrix_d* src ///< Defined for \ref constructor_inverse_Matrix_d.
	);

/** \brief Constructor for a \ref Matrix_d\* from a matrix-matrix multiplication.
 *  \return Standard. */
struct Matrix_d* constructor_mm_Matrix_d
	(const char trans_a_i,                ///< Defined for \ref mm_d.
	 const char trans_b_i,                ///< Defined for \ref mm_d.
	 const double alpha,                  ///< Defined for \ref mm_d.
	 const double beta,                   ///< Defined for \ref mm_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref mm_d.
	 const struct const_Matrix_d*const b, ///< Defined for \ref mm_d.
	 const char layout                    ///< The `layout` of the constructed \ref Matrix_d.
	);

/** \brief `const` version of \ref constructor_mm_Matrix_d.
 *  \return Standard. */
const struct const_Matrix_d* constructor_mm_const_Matrix_d
	(const char trans_a_i,                ///< Defined for \ref constructor_mm_Matrix_d.
	 const char trans_b_i,                ///< Defined for \ref constructor_mm_Matrix_d.
	 const double alpha,                  ///< Defined for \ref constructor_mm_Matrix_d.
	 const double beta,                   ///< Defined for \ref constructor_mm_Matrix_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref constructor_mm_Matrix_d.
	 const struct const_Matrix_d*const b, ///< Defined for \ref constructor_mm_Matrix_d.
	 const char layout                    ///< Defined for \ref constructor_mm_Matrix_d.
	);

/// \brief Set a \ref Matrix_d\* from a sub range of a \ref Multiarray_d\*.
void set_Matrix_from_Multiarray_d
	(struct Matrix_d* dest,            ///< The destination.
	 struct Multiarray_d* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Matrix_from_Multiarray_d.
void set_const_Matrix_from_Multiarray_d
	(const struct const_Matrix_d* dest,    ///< Defined for \ref set_Matrix_from_Multiarray_d.
	 const struct const_Multiarray_d* src, ///< Defined for \ref set_Matrix_from_Multiarray_d.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref set_Matrix_from_Multiarray_d.
	);

/// \brief Set a \ref Matrix_d\* from an entry of a \ref Multiarray_Matrix_d\*.
void set_Matrix_from_Multiarray_Matrix_d
	(struct Matrix_d* dest,            ///< The destination.
	 struct Multiarray_Matrix_d* src,  ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Matrix_from_Multiarray_Matrix_d.
void set_const_Matrix_from_Multiarray_Matrix_d
	(const struct const_Matrix_d* dest,           ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_d.
	 const struct const_Multiarray_Matrix_d* src, ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_d.
	 const ptrdiff_t*const sub_indices            ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_d.
	);

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Matrix_d\*.
void destructor_Matrix_d
	(struct Matrix_d* a ///< Standard.
	);

/// \brief Destructs a \ref const_Matrix_d\*.
void destructor_const_Matrix_d
	(const struct const_Matrix_d* a ///< Standard.
	);

/// \brief Destructs a \ref Matrix_i\*.
void destructor_Matrix_i
	(struct Matrix_i* a ///< Standard.
	);

/// \brief Destructs a \ref const_Matrix_i\*.
void destructor_const_Matrix_i
	(const struct const_Matrix_i* a ///< Standard.
	);

#endif // DPG__Matrix_h__INCLUDED
