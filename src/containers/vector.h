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

#ifndef DPG__vector_h__INCLUDED
#define DPG__vector_h__INCLUDED
/** \file
 *  \brief Provides Vector_\* containers and related functions.
 *
 *  Potentially relevant comments may be found in \ref multiarray.h.
 *
 *  Vectors are 1D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

#include "vector_constructors.h"
#include "vector_math.h"
#include "vector_print.h"

/// \brief Vector (`int`).
struct Vector_i {
	ptrdiff_t ext_0; ///< Defined in \ref Matrix_d.

	bool owns_data; ///< Defined in \ref Multiarray_d.
	int* data;      ///< Defined in \ref Multiarray_d.
};

/// \brief `const` version of \ref Vector_i.
struct const_Vector_i {
	const ptrdiff_t ext_0; ///< See brief.

	const bool owns_data; ///< See brief.
	const int*const data; ///< See brief.
};

/// \brief Vector (`double`).
struct Vector_d {
	ptrdiff_t ext_0; ///< Defined in \ref Matrix_d.

	bool owns_data; ///< Defined in \ref Multiarray_d.
	double* data;   ///< Defined in \ref Multiarray_d.
};

/// \brief `const` version of \ref Vector_d.
struct const_Vector_d {
	const ptrdiff_t ext_0; ///< See brief.

	const bool owns_data;    ///< See brief.
	const double*const data; ///< See brief.
};

// Interface functions ********************************************************************************************** //

/** \brief Reorder a \ref Vector_i based on the provided ordering.
 *	\warning This is not currently done in place.
 */
void reorder_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const int*const ordering ///< The ordering.
	);

/// \brief Resize a \ref Vector_i\*.
void resize_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const ptrdiff_t ext_0    ///< New value for ext_0.
	);

/// \brief Set all data entries to zero.
void set_to_zero_Vector_i
	(struct Vector_i*const a ///< Standard.
	);

/// \brief Set data entries to those of the src data.
void set_to_data_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const int*const data_src ///< The source data.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const int val            ///< The value.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Vector_d
	(struct Vector_d*const a, ///< Standard.
	 const double val         ///< The value.
	);

/// \brief Sort the data of the \ref Vector_i\*.
void sort_Vector_i
	(struct Vector_i* a ///< Standard.
	);

/** \brief See return.
 *  \return The sum of the components of the \ref Vector_i\*. */
int sum_Vector_i
	(const struct Vector_i* a ///< Standard.
	);

/** \brief See return.
 *  \return The product of the components of the \ref Vector_i\*. */
ptrdiff_t prod_Vector_i
	(const struct Vector_i* a ///< Standard.
	);

/** \brief `const` version of \ref prod_Vector_i.
 *  \return See brief. */
ptrdiff_t prod_const_Vector_i
	(const struct const_Vector_i* a ///< Standard.
	);

/**	\brief See return.
 *	\return `bool` indicating whether the \ref Vector_i\* inputs are equal. */
bool check_equal_Vector_i
	(const struct Vector_i*const a, ///< Input 1.
	 const struct Vector_i*const b  ///< Input 2.
	);

/**	\brief See return.
 *	\return `bool` indicating whether the \ref Vector_i\* and `int*` inputs are equal. */
bool check_equal_Vector_i_i
	(const struct Vector_i*const a, ///< The vector input.
	 const int* data_b              ///< The data input.
	);

/** \brief Comparison function for std::qsort between \ref Vector_i\*\* `a` and `b`.
 *	\return The lexicographical comparison of `a` and `b`.
 *
 *	\note Input Vectors must be have sorted data.
 */
int cmp_Vector_i
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/// \brief Copy the data of the `src` to that of the `dest` \ref Vector_i\*.
void copy_data_Vector_i_Vector_i
	(const struct Vector_i*const src, ///< Standard.
	 struct Vector_i*const dest       ///< Standard.
	);

/** \brief Push a value to the back of a \ref Vector_i\*.
 *	The vector may optionally be sorted and/or accept only unique values. */
void push_back_Vector_i
	(struct Vector_i*const src, ///< The source vector.
	 const int val,             ///< The value to be added.
	 const bool sorted,         ///< Flag for whether the vector should be sorted.
	 const bool unique          ///< Flag for whether the value should only be added if it is unique.
	);

/** \brief Check if the input value is present in the optionally sorted src \ref Vector_i\*.
 *	\return `true` if found. */
bool find_val_Vector_i
	(const struct const_Vector_i*const src, ///< The input vector.
	 const int val,                         ///< The value to find.
	 const bool sorted                      ///< Flag indicating whether the vector is sorted.
	);

#endif // DPG__vector_h__INCLUDED
