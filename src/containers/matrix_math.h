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
struct const_Matrix_d;

/**	\brief Compute the norm of the specified row of the input \ref Matrix_d.
 *	\return See brief. */
double compute_norm_Matrix_d_row
	(const ptrdiff_t row,           ///< The row.
	 const struct Matrix_d*const a, ///< The input matrix.
	 const char*const norm_type     ///< The norm type.
	);

///	\brief Transpose the \ref Matrix_d\* optionally leaving the values of the extents unchanged in `mem_only = true`.
void transpose_Matrix_d
	(struct Matrix_d* a, ///< Matrix to be transposed.
	 const bool mem_only ///< Flag for whether only the memory should be transposed (with ext_0/ext_1 unchanged).
	);

/**	\brief Compute the (m)atrix-(m)atrix multiplication of input `double` matrices.
 *
 *	Computes: c = alpha*op(a)*op(b)+beta*c using the [cblas_dgemm function] is used for the computation.
 *
 *	op(): 'N'o transpose, 'T'ranspose.
 *
 *	<!-- References: -->
 *	[cblas_dgemm]: https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemm
 */
void mm_d
	(const char trans_a_i,                ///< Operator for input `a`. Options: 'N'o transpose, 'T'ranspose.
	 const char trans_b_i,                ///< Operator for input `b`. Options: 'N'o transpose, 'T'ranspose.
	 const double alpha,                  ///< Multiplicative constant.
	 const double beta,                   ///< Multiplicative constant.
	 const struct const_Matrix_d*const a, ///< Input \ref const_Matrix_d\* `a`.
	 const struct const_Matrix_d*const b, ///< Input \ref const_Matrix_d\* `b`.
	 struct Matrix_d*const c              ///< Input \ref Matrix_d\* `c`.
	);

///	\brief Specialization of \ref mm_d for `trans_a_i = 'N'`, `trans_b_i = 'N'`.
void mm_NN_d
	(const double alpha,                  ///< Defined for \ref mm_d.
	 const double beta,                   ///< Defined for \ref mm_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref mm_d.
	 const struct const_Matrix_d*const b, ///< Defined for \ref mm_d.
	 struct Matrix_d*const c              ///< Defined for \ref mm_d.
	);

#endif // DPG__matrix_math_h__INCLUDED
