// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__matrix_math_h__INCLUDED
#define DPG__matrix_math_h__INCLUDED
/**	\file
 *	\brief Provides Matrix_\* math functions.
 */

#include <stddef.h>
#include <stdbool.h>

struct Matrix_d;

/** \brief Compute the norm of the specified row of the input \ref Matrix_d.
 *	\return See brief. */
double compute_norm_Matrix_d_row
	(const ptrdiff_t row,           ///< The row.
	 const struct Matrix_d*const a, ///< The input matrix.
	 const char*const norm_type     ///< The norm type.
	);

/// \brief Transpose the \ref Matrix_d\* optionally leaving the values of the extents unchanged in `mem_only = true`.
void transpose_Matrix_d
	(struct Matrix_d* a, ///< Matrix to be transposed.
	 const bool mem_only ///< Flag for whether only the memory should be transposed (with ext_0/ext_1 unchanged).
	);

#endif // DPG__matrix_math_h__INCLUDED
