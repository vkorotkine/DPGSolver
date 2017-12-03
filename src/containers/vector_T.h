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

/// \brief Templated Vector.
struct Vector_T {
	ptrdiff_t ext_0; ///< Defined in \ref Matrix_T.

	bool owns_data; ///< Defined in \ref Matrix_T.
	Type* data;     ///< Defined in \ref Matrix_T.
};

/// \brief `const` version of \ref Vector_T.
struct const_Vector_T { ///\{
	const ptrdiff_t ext_0;

	const bool owns_data;
	const Type*const data;
}; ///\}

// Interface functions ********************************************************************************************** //

/** \brief Reorder a \ref Vector_T based on the provided ordering.
 *  \warning This is not currently done in place.
 */
void reorder_Vector_T
	(struct Vector_T*const a, ///< Standard.
	 const int*const ordering ///< The ordering.
	);

/// \brief Resize a \ref Vector_T\*.
void resize_Vector_T
	(struct Vector_T*const a, ///< Standard.
	 const ptrdiff_t ext_0    ///< New value for ext_0.
	);

/// \brief Set all data entries to zero.
void set_to_zero_Vector_T
	(struct Vector_T*const a ///< Standard.
	);

/// \brief Set data entries to those of the src data.
void set_to_data_Vector_T
	(struct Vector_T*const a, ///< Standard.
	 const Type*const data_src ///< The source data.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Vector_T
	(struct Vector_T*const a, ///< Standard.
	 const Type val            ///< The value.
	);

/// \brief Sort the data of the \ref Vector_T\*.
void sort_Vector_T
	(struct Vector_T* a ///< Standard.
	);

/** \brief See return.
 *  \return The sum of the components of the \ref Vector_T\*. */
Type sum_Vector_T
	(const struct Vector_T* a ///< Standard.
	);

/** \brief See return.
 *  \return The product of the components of the \ref Vector_T\*. */
ptrdiff_t prod_Vector_T
	(const struct Vector_T* a ///< Standard.
	);

/** \brief `const` version of \ref prod_Vector_T.
 *  \return See brief. */
ptrdiff_t prod_const_Vector_T
	(const struct const_Vector_T* a ///< Standard.
	);

/** \brief See return.
 *  \return `bool` indicating whether the \ref Vector_T\* inputs are equal. */
bool check_equal_Vector_T
	(const struct Vector_T*const a, ///< Input 1.
	 const struct Vector_T*const b  ///< Input 2.
	);

/** \brief See return.
 *  \return `bool` indicating whether the \ref Vector_T\* and `Type*` inputs are equal. */
bool check_equal_Vector_T_T
	(const struct Vector_T*const a, ///< The vector input.
	 const Type* data_b              ///< The data input.
	);

/** \brief Comparison function for std::qsort between \ref Vector_T\*\* `a` and `b`.
 *  \return The lexicographical comparison of `a` and `b`.
 *
 *  \note Input Vectors must be have sorted data.
 */
int cmp_Vector_T
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/// \brief Copy the data of the `src` to that of the `dest` \ref Vector_T\*.
void copy_data_Vector_T_Vector_T
	(const struct Vector_T*const src, ///< Standard.
	 struct Vector_T*const dest       ///< Standard.
	);

/** \brief Push a value to the back of a \ref Vector_T\*.
 *  The vector may optionally be sorted and/or accept only unique values. */
void push_back_Vector_T
	(struct Vector_T*const src, ///< The source vector.
	 const Type val,             ///< The value to be added.
	 const bool sorted,         ///< Flag for whether the vector should be sorted.
	 const bool unique          ///< Flag for whether the value should only be added if it is unique.
	);

/** \brief Check if the input value is present in the optionally sorted src \ref Vector_T\*.
 *  \return `true` if found. */
bool find_val_Vector_T
	(const struct const_Vector_T*const src, ///< The input vector.
	 const Type val,                         ///< The value to find.
	 const bool sorted                      ///< Flag indicating whether the vector is sorted.
	);
