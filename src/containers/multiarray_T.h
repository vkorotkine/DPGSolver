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
 *  \brief Provides Multiarray_\* containers and related functions.
 */

#include <stddef.h>
#include <stdbool.h>

struct const_Vector_R;

/// \brief Templated Multiarray.
struct Multiarray_T {
	char layout; ///< The layout may be 'R'ow or 'C'olumn major.

	int order;          ///< Number of dimensions.
	ptrdiff_t* extents; ///< Size of array in each dimension.

	bool owns_data; /**< Flag for whether the data should be freed in the destructor. This would be false if a move
	                     constructor was used. */
	Type* data; ///< The data.
};

/// \brief `const` version of \ref Multiarray_T.
struct const_Multiarray_T { ///\{
	const char layout;

	const int order;
	const ptrdiff_t*const extents;

	const bool owns_data;
	const Type*const data;
}; ///\}

/// \brief Templated Multiarray of \ref Vector_T\*s.
struct Multiarray_Vector_T {
	int order;          ///< Defined in \ref Multiarray_T.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_T.

	bool owns_data;         ///< Defined in \ref Multiarray_T.
	struct Vector_T** data; ///< Defined in \ref Multiarray_T.
};

/// \brief `const` version of \ref Multiarray_Vector_T.
struct const_Multiarray_Vector_T { ///\{
	const int order;
	const ptrdiff_t*const extents;

	const bool owns_data;
	const struct const_Vector_T*const*const data;
}; ///\}

/// \brief Templated Multiarray of \ref Matrix_T\*s.
struct Multiarray_Matrix_T {
	int order;          ///< Defined in \ref Multiarray_T.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_T.

	bool owns_data;         ///< Defined in \ref Multiarray_T.
	struct Matrix_T** data; ///< Defined in \ref Multiarray_T.
};

/// \brief `const` version of \ref Multiarray_Matrix_T.
struct const_Multiarray_Matrix_T { ///\{
	const int order;
	const ptrdiff_t*const extents;

	const bool owns_data;
	const struct const_Matrix_T*const*const data;
}; ///\}

// Interface functions ********************************************************************************************** //

/** \brief Get pointer to row of row-major \ref Multiarray_T\* of order 2.
 *  \return Pointer to the first entry of the row. */
Type* get_row_Multiarray_T
	(const ptrdiff_t row,         ///< Desired row.
	 const struct Multiarray_T* a ///< Multiarray.
	);

/** \brief `const` version of \ref get_row_Multiarray_T.
 *  \return See brief. */
const Type* get_row_const_Multiarray_T
	(const ptrdiff_t row,               ///< Defined for \ref get_row_Multiarray_T.
	 const struct const_Multiarray_T* a ///< Defined for \ref get_row_Multiarray_T.
	);

/** \brief Get pointer to col of col-major \ref Multiarray_T\*.
 *  \return See brief. */
Type* get_col_Multiarray_T
	(const ptrdiff_t col,   ///< Desired column.
	 struct Multiarray_T* a ///< Multiarray.
	);

/** \brief `const` version of \ref get_col_Multiarray_T.
 *  \return See brief. */
const Type* get_col_const_Multiarray_T
	(const ptrdiff_t col,               ///< Defined for \ref get_col_Multiarray_T.
	 const struct const_Multiarray_T* a ///< Defined for \ref get_col_Multiarray_T.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Multiarray_T
	(struct Multiarray_T*const a, ///< Standard.
	 const Type val               ///< The value.
	);

/// \brief Set the values of the \ref Multiarray_Vector_T based on the input `Type*` data.
void set_Multiarray_Vector_T_T
	(struct Multiarray_Vector_T* a, ///< Standard.
	 const Type* data_V,             ///< Input data for the Vectors.
	 const int*const ext_V          ///< Defined in \ref constructor_copy_Multiarray_Vector_T_T.
	);

/// \brief Set the values of the output \ref Multiarray_T to those of the input \ref Multiarray_T.
void set_Multiarray_T
	(struct Multiarray_T* a_o,            ///< Output multiarray.
	 const struct const_Multiarray_T* a_i ///< Input multiarray.
	);

/// \brief Set the data of the \ref Multiarray_T\* container from that of the \ref Multiarray_R\* container.
void set_Multiarray_T_Multiarray_R
	(struct Multiarray_T* a,            ///< Multiarray with data to be set.
	 const struct const_Multiarray_R* b ///< Multiarray from which to take data.
	);

/** \brief Sort the data of the \ref Multiarray_Vector_T\*.
 *  \return Optionally return indices or `NULL`.
 */
struct Vector_i* sort_Multiarray_Vector_T
	(struct Multiarray_Vector_T* a, ///< Standard.
	 const bool return_indices      ///< Flag for whether the indices should also be returned.
	);

/** \brief Collapse a \ref Multiarray_Vector_T\* into a \ref Vector_T\* with copied data.
 *  \return The \ref Vector_T\*. */
struct Vector_T* collapse_Multiarray_Vector_T
	(const struct Multiarray_Vector_T*const src ///< The source.
	);

/// \brief Resize a \ref Multiarray_T\*, clearing any existing data entries.
void resize_Multiarray_T
	(struct Multiarray_T* a,  ///< The multiarray.
	 const int order,         ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t* extents ///< Defined in \ref Multiarray_T.
	);

/** \brief Get a pointer to a \ref const_Vector_T\* from a sub range of a \ref const_Multiarray_Vector_T\*.
 *  \return See brief. */
const struct const_Vector_T* get_const_Multiarray_Vector_T
	(const struct const_Multiarray_Vector_T* src, ///< The source.
	 const ptrdiff_t*const sub_indices            ///< The sub-indices specifying which part of the source to extract.
	);

/** \brief Return a copy of a stack allocated \ref const_Vector_T holding the data of the input multiarray of order 1.
 *  \return See brief. */
struct const_Vector_T interpret_const_Multiarray_as_Vector_T
	(const struct const_Multiarray_T* a_Ma ///< The input multiarray.
	);

/** \brief Return a copy of a stack allocated \ref Matrix_T holding the data of the input multiarray of order 2.
 *  \return See brief. */
struct Matrix_T interpret_Multiarray_as_Matrix_T
	(const struct Multiarray_T* a_Ma ///< The input multiarray.
	);
