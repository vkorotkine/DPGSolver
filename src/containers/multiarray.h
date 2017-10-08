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

#ifndef DPG__multiarray_h__INCLUDED
#define DPG__multiarray_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* containers and related functions.
 *
 *  \section s1_Multi General
 *
 *  \subsection s11_Multi Standard Data Types
 *
 *  For standard data types, the Multiarray container is intended to be used as a higher-dimensional matrix where
 *  move constructors are used to form matrix/vector containers for appropriate sub-blocks. As the data is stored
 *  contiguously in memory, the Multiarray may also be acted on over multiple "dimensions" at once.
 *
 *  Note that the `layout` parameter denotes the layout of the first two extents of the Multiarray (i.e. the layout of
 *  the Matrix interpretation of the 2d sub-blocks of the Multiarray). To act over multiple "dimensions" (more than 2)
 *  at once, the layout **must** be column-major.
 *
 *  \subsubsection s111_Multi Supported Types
 *
 *  Following the recommendation of many c++ experts, containers of unsigned integer types are **not** supported
 *  even if the variable to be represented is always positive. For more details, see the Eigen FAQ section on
 *  [Why Eigen's API is using signed integers for sizes, indices, etc.?][eigen_signed] and the linked video and
 *  posts.  [This][SO_signed] SO discussion provides a summary of some of the points made in the video referred to
 *  above if it is no longer accessible.
 *
 *  For similar reasons, the type used for array indexing and storing array sizes is `ptrdiff_t` and not `size_t`.
 *  It was mentioned that the use of `size_t` for container sizes (and indexing by association) is now seen as a
 *  mistake in the c++ standard.
 *
 *  \subsection s12_Multi Defined Types
 *
 *  \subsubsection s122_Multi Multiarray Specializations
 *
 *  Two specializations of the Multiarray exist: Matrix (2D Multiarray) and Vector (1D Multiarray). These containers
 *  are used when the data is most intuitively considered to be of the given form (e.g. mathematical operators are
 *  Matrices).
 *
 *  Further, functionality is provided for storing Multiarrays of these specialized types (i.e. Multiarrays of
 *  Matrices and Vectors). When comparing with multiply dereferenced specialized containers (e.g. struct Matrix**),
 *  this results in the fundamental advantage of containers carrying around their size information, easing the
 *  burden on the developer. This implies that containers with multiple levels of dereferencing should never be
 *  used.
 *
 *  As the intended functionality of the specialized Multiarrays does not include usage as higher-dimensional
 *  matrices for contiguous data storage, the `layout` parameter is omitted from these containers.
 *
 *  <!-- References: -->
 *  [eigen_signed]: http://eigen.tuxfamily.org/index.php?title=FAQ#Why_Eigen.27s_API_is_using_signed_integers_for_sizes.2C_indices.2C_etc..3F
 *  [SO_signed]: https://stackoverflow.com/a/18796234/5983549
 */

#include <stddef.h>
#include <stdbool.h>

#include "multiarray_constructors.h"
#include "multiarray_math.h"
#include "multiarray_print.h"

/// \brief Multiarray (`double`).
struct Multiarray_d {
	char layout; ///< The layout may be 'R'ow or 'C'olumn major.

	int order;          ///< Number of dimensions.
	ptrdiff_t* extents; ///< Size of array in each dimension.

	bool owns_data; /**< Flag for whether the data should be freed in the destructor. This would be false if a move
	                     constructor was used. */
	double* data; ///< The data.
};

/// \brief Multiarray (`const double`).
struct const_Multiarray_d {
	const char layout; ///< Defined in \ref Multiarray_d.

	const int order;               ///< Defined in \ref Multiarray_d.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_d.

	const bool owns_data;    ///< Defined in \ref Multiarray_d.
	const double*const data; ///< Defined in \ref Multiarray_d.
};

/// \brief Multiarray (`Vector_i*`).
struct Multiarray_Vector_i {
	int order;          ///< Defined in \ref Multiarray_d.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_d.

	bool owns_data;         ///< Defined in \ref Multiarray_d.
	struct Vector_i** data; ///< Defined in \ref Multiarray_d.
};

/// \brief Multiarray (`const Vector_i*`).
struct const_Multiarray_Vector_i {
	const int order;               ///< Defined in \ref Multiarray_d.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_d.

	const bool owns_data;                         ///< Defined in \ref Multiarray_d.
	const struct const_Vector_i*const*const data; ///< Defined in \ref Multiarray_d.
};

/// \brief Multiarray (`Matrix_d*`).
struct Multiarray_Matrix_d {
	int order;          ///< Defined in \ref Multiarray_d.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_d.

	bool owns_data;         ///< Defined in \ref Multiarray_d.
	struct Matrix_d** data; ///< Defined in \ref Multiarray_d.
};

/// \brief Multiarray (`const_Matrix_d*`).
struct const_Multiarray_Matrix_d {
	const int order;               ///< Defined in \ref Multiarray_d.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_d.

	const bool owns_data;                         ///< Defined in \ref Multiarray_d.
	const struct const_Matrix_d*const*const data; ///< Defined in \ref Multiarray_d.
};

// Interface functions ********************************************************************************************** //

/** \brief Computes the `size`, which is the product of the `extents`.
 *  \return See brief. */
ptrdiff_t compute_size
	(const int order,              ///< \ref Multiarray_d::order.
	 const ptrdiff_t*const extents ///< \ref Multiarray_d::extents.
	);

/** \brief Get pointer to row of row-major \ref Multiarray_d\* of order 2.
 *  \return Pointer to the first entry of the row. */
double* get_row_Multiarray_d
	(const ptrdiff_t row,         ///< Desired row.
	 const struct Multiarray_d* a ///< Multiarray.
	);

/** \brief `const` version of \ref get_row_Multiarray_d.
 *  \return See brief. */
const double* get_row_const_Multiarray_d
	(const ptrdiff_t row,               ///< Defined for \ref get_row_Multiarray_d.
	 const struct const_Multiarray_d* a ///< Defined for \ref get_row_Multiarray_d.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Multiarray_d
	(struct Multiarray_d*const a, ///< Standard.
	 const double val             ///< The value.
	);

/// \brief Set the values of the \ref Multiarray_Vector_i based on the input `int*` data.
void set_Multiarray_Vector_i_i
	(struct Multiarray_Vector_i* a, ///< Standard.
	 const int* data_V,             ///< Input data for the Vectors.
	 const int*const ext_V          ///< Defined in \ref constructor_copy_Multiarray_Vector_i_i.
	);

/** \brief Sort the data of the \ref Multiarray_Vector_i\*.
 *  \return Optionally return indices or `NULL`.
 */
struct Vector_i* sort_Multiarray_Vector_i
	(struct Multiarray_Vector_i* a, ///< Standard.
	 const bool return_indices      ///< Flag for whether the indices should also be returned.
	);

/** \brief Collapse a \ref Multiarray_Vector_i\* into a \ref Vector_i\* with copied data.
 *  \return The \ref Vector_i\*. */
struct Vector_i* collapse_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const src ///< The source.
	);

/** \brief Compute the index of the data of the sub-Vector based on the sub-indices.
 *  \deprecated Replace with \ref compute_index_sub_container.
 *  \return See brief. */
ptrdiff_t compute_index_sub_vector
	(const int order,                  ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents,    ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const sub_indices ///< The sub indices for the entries of order > 1.
	);

/** \brief Compute the index of the data of the sub-Matrix in a Multiarray based on the sub-indices.
 *  \deprecated Replace with \ref compute_index_sub_container.
 *  \return See brief. */
ptrdiff_t compute_index_sub_matrix
	(const int order,                  ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents,    ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const sub_indices ///< The sub indices for the entries of order > 2.
	);

/** \brief Compute the index of the data of the sub-container based on the sub-indices.
 *  \return See brief. */
ptrdiff_t compute_index_sub_container
	(const int order_i,                ///< The input order.
	 const int order_o,                ///< The output order.
	 const ptrdiff_t*const extents,    ///< The input extents.
	 const ptrdiff_t*const sub_indices ///< The sub indices for the highest `order_o` entries.
	);

/** \brief Version of of \ref compute_index_sub_container accepting `ptrdiff_t*` and `int*` inputs.
 *  \return See brief. */
ptrdiff_t compute_index_sub_container_pi
	(const int order_i,             ///< Defined for \ref compute_index_sub_container.
	 const int order_o,             ///< Defined for \ref compute_index_sub_container.
	 const ptrdiff_t*const extents, ///< Defined for \ref compute_index_sub_container.
	 const int*const sub_indices    ///< Defined for \ref compute_index_sub_container.
	);

/// \brief Resize a \ref Multiarray_d\*, clearing any existing data entries.
void resize_Multiarray_d
	(struct Multiarray_d* a,  ///< The multiarray.
	 const int order,         ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t* extents ///< Defined in \ref Multiarray_d.
	);

#endif // DPG__multiarray_h__INCLUDED
