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
 *  \brief Provides Matrix_\* math functions.
 */

#include <stddef.h>
#include <stdbool.h>

struct Matrix_T;
struct Vector_d;
struct const_Matrix_T;
struct const_Vector_d;
struct const_Vector_i;

/** \brief Compute the norm of the specified row of the input \ref Matrix_T.
 *  \return See brief. */
double compute_norm_Matrix_T_row
	(const ptrdiff_t row,           ///< The row.
	 const struct Matrix_T*const a, ///< The input matrix.
	 const char*const norm_type     ///< The norm type.
	);

/// \brief Transpose the \ref Matrix_T\* optionally leaving the values of the extents unchanged if `mem_only = true`.
void transpose_Matrix_T
	(struct Matrix_T* a, ///< Matrix to be transposed.
	 const bool mem_only ///< Flag for whether only the memory should be transposed (with ext_0/ext_1 unchanged).
	);

/// \brief Invert a sub-block of the input \ref Matrix_T\*.
void invert_sub_block_Matrix_T
	(struct Matrix_T* a,   ///< The input matrix.
	 const ptrdiff_t row0, ///< Index of the first row in the input matrix.
	 const ptrdiff_t col0, ///< Index of the first col in the input matrix.
	 const ptrdiff_t ext   ///< The number of rows/cols of the sub-matrix.
	);

/// \brief Scale the \ref Matrix_T\* by a constant value.
void scale_Matrix_T
	(struct Matrix_T* a, ///< The matrix.
	 const double val    ///< The value by which to scale.
	);

/** \brief Permute the the input matrix according to the input permutation.
 *
 *  The permutation is performed using the [gsl_permute_matrix] function and is from the left for a row-major input and
 *  from the right for a column major input.
 *
 *  <!-- References: -->
 *  [gsl_permuate_matrix]: https://www.gnu.org/software/gsl/doc/html/permutation.html#c.gsl_permute_matrix
 */
void permute_Matrix_T
	(struct Matrix_T* a, ///< Matrix to be permuted.
	 const ptrdiff_t* p  ///< Permutation indices.
	);

/// \brief Call \ref permute_Matrix_T using the indices of the vector as the permutation indices.
void permute_Matrix_T_V
	(struct Matrix_T* a,              ///< Defined for \ref permute_Matrix_T.
	 const struct const_Vector_i* p_V ///< Vector of permutation indices.
	);

/** \brief Compute the (m)atrix-(m)atrix multiplication of input `double` matrices.
 *
 *  Computes: c = alpha*op(a)*op(b)+beta*c using the [cblas_dgemm function][cblas_dgemm] for the computation.
 *
 *  op(): 'N'o transpose, 'T'ranspose.
 *
 *  <!-- References: -->
 *  [cblas_dgemm]: https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemm
 */
void mm_T
	(const char trans_a_i,                ///< Operator for input `a`. Options: 'N'o transpose, 'T'ranspose.
	 const char trans_b_i,                ///< Operator for input `b`. Options: 'N'o transpose, 'T'ranspose.
	 const Real alpha,                  ///< Multiplicative constant.
	 const Real beta,                   ///< Multiplicative constant.
	 const struct const_Matrix_T*const a, ///< Input \ref const_Matrix_T\* `a`.
	 const struct const_Matrix_T*const b, ///< Input \ref const_Matrix_T\* `b`.
	 struct Matrix_T*const c              ///< Input \ref Matrix_T\* `c`.
	);

/** \brief Compute the (m)atrix-(v)ector multiplication of input `double` containers.
 *
 *  Computes: c = alpha*op(a)*b+beta*c using the [cblas_dgemv function][cblas_dgemv] for the computation.
 *
 *  op(): 'N'o transpose, 'T'ranspose.
 *
 *  As the vectors do not have a layout, the `layout` parameter for the dgemv call is chosen as that of the input
 *  matrix.
 *
 *  <!-- References: -->
 *  [cblas_dgemv]: https://software.intel.com/en-us/mkl-developer-reference-c-cblas-gemv
 */
void mv_T
	(const char trans_a_i,                ///< Operator for input `a`. Options: 'N'o transpose, 'T'ranspose.
	 const Real alpha,                  ///< Multiplicative constant.
	 const Real beta,                   ///< Multiplicative constant.
	 const struct const_Matrix_T*const a, ///< Input \ref const_Matrix_T\* `a`.
	 const struct const_Vector_d*const b, ///< Input \ref const_Vector_d\* `b`.
	 struct Vector_d*const c              ///< Input \ref Vector_d\* `c`.
	);

/** \brief Computes the matrix-"matrix" multiplication of a matrix with a vector interpreted as a diagonal matrix
 *         in-place.
 *
 *  if (side == 'L')
 *  	computes: A = alpha*diag(b)*A.
 *  else if (side == 'R')
 *  	computes: A = alpha*A*diag(b).
 */
void scale_Matrix_by_Vector_T
	(const char side,                     ///< The side from which to apply the vector as a diagonal matrix.
	 const Real alpha,                  ///< Multiplicative constant.
	 struct Matrix_T*const a,             ///< Input \ref Matrix_T\*.
	 const struct const_Vector_d*const b, ///< Input \ref const_Vector_d\*.
	 const bool invert_diag               /**< Flag for whether the diagonal entries should be inverted before
	                                       *   application. */
	);

/** \brief Computes the matrix-"matrix" multiplication of a matrix with a vector interpreted as a diagonal matrix.
 *
 *  if (side == 'L')
 *  	computes: C = alpha*diag(b)*A + beta*C.
 *  else if (side == 'R')
 *  	computes: C = alpha*A*diag(b) + beta*C.
 *
 */
void mm_diag_T
	(const char side,                     ///< The side from which to apply the vector as a diagonal matrix.
	 const Real alpha,                  ///< Multiplicative constant.
	 const Real beta,                   ///< Multiplicative constant.
	 const struct const_Matrix_T*const a, ///< Input \ref const_Matrix_T\*.
	 const struct const_Vector_d*const b, ///< Input \ref const_Vector_d\*.
	 struct Matrix_T* c,                  ///< Input \ref Matrix_T\*.
	 const bool invert_diag               /**< Flag for whether the diagonal entries should be inverted before
	                                       *   application. */
	);

/// \brief Reinterpret the input \ref const_Matrix_T\* as having the input extents.
void reinterpret_const_Matrix_T
	(const struct const_Matrix_T* a, ///< The input matrix.
	 const ptrdiff_t ext_0,          ///< The new value for `ext_0`.
	 const ptrdiff_t ext_1           ///< The new value for `ext_1`.
	);
