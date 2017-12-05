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
 *  \brief Provides the templated \ref Operator related functions.
 *
 *  In most cases, only the standard, dense matrix operator will be available, however, for tensor-product elements, the
 *  sub-operators are also always available. The CSR sparse matrix format can be made available for operators known to
 *  have significant sparsity.
 *
 *  When used for matrix multiplication, the supported operator formats are denoted by the following `char`s:
 *  - 'd'efault -> standard.
 *  - 's'tandard;
 *  - 't'ensor-product;
 *  - 'c'sr.
 */

#include <stddef.h>
#include "definitions_core.h"

struct const_Multiarray_T;
struct Multiarray_T;
struct const_Multiarray_Matrix_T;

// Interface functions ********************************************************************************************** //
// Constructors ***************************************************************************************************** //

/** \brief Constructor for a \ref Multiarray_T\* from the operator-multiarray multiplication using the input operator
 *         format.
 *  \return See brief. */
struct Multiarray_T* constructor_mm_NN1_Operator_Multiarray_T
	(const struct Operator* op,    ///< \ref Operator.
	 const struct Multiarray_T* b, ///< The input multiarray.
	 const char layout_c,          ///< The desired layout of the output `c`.
	 const char op_format,         ///< The operator format to be used.
	 const int order_sub_ma,       ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b   ///< The sub-indices of the `b` multiarray if required.
	);

/** \brief `const` version of \ref  constructor_mm_NN1_Operator_Multiarray_T.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_mm_NN1_Operator_const_Multiarray_T
	(const struct Operator* op,          ///< See brief.
	 const struct const_Multiarray_T* b, ///< See brief.
	 const char layout_c,                ///< See brief.
	 const char op_format,               ///< See brief.
	 const int order_sub_ma,             ///< See brief.
	 const ptrdiff_t* sub_inds_b         ///< See brief.
	);

// General functions ************************************************************************************************ //

/// \brief Set the input arrays to hold the number of rows and columns of the tensor-product operators.
void set_ops_tp_n_rows_cols_T
	(int n_rows_sub[DMAX],                          ///< The number of rows of each sub-operator.
	 int n_cols_sub[DMAX],                          ///< The number of columns of each sub-operator.
	 const struct const_Multiarray_Matrix_T* ops_tp ///< The tensor-product operators.
	);

// Math functions *************************************************************************************************** //

/** \brief Apply a matrix-matrix multiplication of an operator with a \ref const_Multiarray_T\* using the input operator
 *         format, asserting that the two input multiarrays have column-major layout.
 *
 *  See comments in \ref mm_NNC_Multiarray_T for the preset matrix-matrix multiplication parameters.
 */
void mm_NNC_Operator_Multiarray_T
	(const Real alpha,                   ///< Defined for \ref mm_NNC_Multiarray_T.
	 const Real beta,                    ///< Defined for \ref mm_NNC_Multiarray_T.
	 const struct Operator* op,          ///< \ref Operator.
	 const struct const_Multiarray_T* b, ///< The input multiarray.
	 struct Multiarray_T* c,             ///< The output multiarray.
	 const char op_format,               ///< The operator format to be used.
	 const int order_sub_ma,             ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b,        ///< The sub-indices of the `b` multiarray if required.
	 const ptrdiff_t* sub_inds_c         ///< The sub-indices of the `c` multiarray if required.
	);

/** \brief Apply a matrix-matrix multiplication of an operator with a \ref const_Multiarray_T\* using the input operator
 *         format, asserting that the two input multiarrays have column-major layout.
 *
 *  See comments in constructor_mm_NN1C_Matrix_T for the preset matrix-matrix multiplication parameters.
 */
void mm_NN1C_Operator_Multiarray_T
	(const struct Operator* op,          ///< \ref Operator.
	 const struct const_Multiarray_T* b, ///< The input multiarray.
	 struct Multiarray_T* c,             ///< The output multiarray.
	 const char op_format,               ///< The operator format to be used.
	 const int order_sub_ma,             ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b,        ///< The sub-indices of the `b` multiarray if required.
	 const ptrdiff_t* sub_inds_c         ///< The sub-indices of the `c` multiarray if required.
	);

/** \brief Apply a matrix-matrix multiplication of an operator with a \ref const_Multiarray_T\* using the input operator
 *         format, transposing the multiarrays for the operation if they do not have a column-major layout.
 *
 *  See comments in constructor_mm_NN1C_Matrix_T for the preset matrix-matrix multiplication parameters, excluding
 *  the layout.
 */
void mm_NN1_Operator_Multiarray_T
	(const struct Operator* op,          ///< \ref Operator.
	 const struct const_Multiarray_T* b, ///< The input multiarray.
	 struct Multiarray_T* c,             ///< The output multiarray.
	 const char op_format,               ///< The operator format to be used.
	 const int order_sub_ma,             ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b,        ///< The sub-indices of the `b` multiarray if required.
	 const ptrdiff_t* sub_inds_c         ///< The sub-indices of the `c` multiarray if required.
	);

