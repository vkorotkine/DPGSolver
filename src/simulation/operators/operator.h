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

#ifndef DPG__operator_h__INCLUDED
#define DPG__operator_h__INCLUDED
/** \file
 *  \brief Provides the Operator\* container and related functions.
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

struct const_Multiarray_d;
struct Multiarray_d;

/// Container holding the operator matrices in various formats.
struct Operator {
	const struct const_Matrix_d*const op_std;             ///< The standard dense matrix operator.
	const struct const_Multiarray_Matrix_d*const  ops_tp; ///< The multiarray of tensor-product sub-operators.
	const struct const_Matrix_CSR_d*const op_csr;         ///< The sparse matrix operator in CSR format.
};

/// `mutable` version of \ref Operator
struct mutable_Operator {
	struct Matrix_d* op_std;            ///< The standard dense matrix operator.
	struct Multiarray_Matrix_d* ops_tp; ///< The multiarray of tensor-product sub-operators.
	struct Matrix_CSR_d* op_csr;        ///< The sparse matrix operator in CSR format.
};

// Interface functions ********************************************************************************************** //
// Constructors ***************************************************************************************************** //

/** \brief Constructor for a \ref const_Multiarray_d\* from the operator-multiarray multiplication using the input
 *         operator format.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_mm_NN1_Operator_const_Multiarray_d
	(const struct Operator* op,          ///< \ref Operator.
	 const struct const_Multiarray_d* b, ///< The input multiarray.
	 const char layout_c,                ///< The desired layout of the output `c`.
	 const char op_format,               ///< The operator format to be used.
	 const int order_sub_ma,             ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b         ///< The sub-indices of the `b` multiarray if required.
	);

// Destructors ****************************************************************************************************** //

/// \brief Destructor for a \ref Operator.
void destructor_Operator
	(const struct Operator* op ///< Standard.
	);

/// \brief Destructor for a \ref mutable_Operator.
void destructor_mutable_Operator
	(struct mutable_Operator* op ///< Standard.
	);

// General functions ************************************************************************************************ //

/// \brief Set the input arrays to hold the number of rows and columns of the tensor-product operators.
void set_ops_tp_n_rows_cols
	(int n_rows_sub[DMAX],                          ///< The number of rows of each sub-operator.
	 int n_cols_sub[DMAX],                          ///< The number of columns of each sub-operator.
	 const struct const_Multiarray_Matrix_d* ops_tp ///< The tensor-product operators.
	);

// Math functions *************************************************************************************************** //

/** \brief Apply a matrix-matrix multiplication of an operator with a \ref const_Multiarray_d\* using the input operator
 *         format, asserting that the two input multiarrays have column-major layout.
 *
 *  See comments in \ref mm_NNC_Multiarray_d for the preset matrix-matrix multiplication parameters.
 */
void mm_NNC_Operator_Multiarray_d
	(const double alpha,                 ///< Defined for \ref mm_NNC_Multiarray_d.
	 const double beta,                  ///< Defined for \ref mm_NNC_Multiarray_d.
	 const struct Operator* op,          ///< \ref Operator.
	 const struct const_Multiarray_d* b, ///< The input multiarray.
	 struct Multiarray_d* c,             ///< The output multiarray.
	 const char op_format,               ///< The operator format to be used.
	 const int order_sub_ma,             ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b,        ///< The sub-indices of the `b` multiarray if required.
	 const ptrdiff_t* sub_inds_c         ///< The sub-indices of the `c` multiarray if required.
	);

/** \brief Apply a matrix-matrix multiplication of an operator with a \ref const_Multiarray_d\* using the input operator
 *         format, asserting that the two input multiarrays have column-major layout.
 *
 *  See comments in \ref constructor_mm_NN1C_Matrix_d for the preset matrix-matrix multiplication parameters.
 */
void mm_NN1C_Operator_Multiarray_d
	(const struct Operator* op,          ///< \ref Operator.
	 const struct const_Multiarray_d* b, ///< The input multiarray.
	 struct Multiarray_d* c,             ///< The output multiarray.
	 const char op_format,               ///< The operator format to be used.
	 const int order_sub_ma,             ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b,        ///< The sub-indices of the `b` multiarray if required.
	 const ptrdiff_t* sub_inds_c         ///< The sub-indices of the `c` multiarray if required.
	);

/** \brief Apply a matrix-matrix multiplication of an operator with a \ref const_Multiarray_d\* using the input operator
 *         format, transposing the multiarrays for the operation if they do not have a column-major layout.
 *
 *  See comments in \ref constructor_mm_NN1C_Matrix_d for the preset matrix-matrix multiplication parameters, excluding
 *  the layout.
 */
void mm_NN1_Operator_Multiarray_d
	(const struct Operator* op,          ///< \ref Operator.
	 const struct const_Multiarray_d* b, ///< The input multiarray.
	 struct Multiarray_d* c,             ///< The output multiarray.
	 const char op_format,               ///< The operator format to be used.
	 const int order_sub_ma,             ///< The order of each of the sub-multiarrays to be operated on.
	 const ptrdiff_t* sub_inds_b,        ///< The sub-indices of the `b` multiarray if required.
	 const ptrdiff_t* sub_inds_c         ///< The sub-indices of the `c` multiarray if required.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Operator\* to the terminal displaying entries below the default tolerance as 0.0.
void print_Operator
	(const struct Operator*const a ///< Standard.
	);

/// \brief Print a \ref Operator\* to the terminal displaying entries below the tolerance as 0.0.
void print_Operator_tol
	(const struct Operator*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

#endif // DPG__operator_h__INCLUDED
