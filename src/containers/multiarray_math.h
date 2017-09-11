// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__multiarray_math_h__INCLUDED
#define DPG__multiarray_math_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* math functions.
 */

#include <stdbool.h>

struct Multiarray_d;

/** \brief Transpose the \ref Multiarray_d\* optionally leaving the values of the extents unchanged if
 *         `mem_only = true`.
 *
 *  \note The multiarray must have `order = 2`.
 */
void transpose_Multiarray_d
	(struct Multiarray_d* a, ///< Multiarray to be transposed.
	 const bool mem_only     ///< Flag for whether only the memory should be transposed (with extents unchanged).
	);

#endif // DPG__multiarray_math_h__INCLUDED
