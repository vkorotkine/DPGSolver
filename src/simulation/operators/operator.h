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
 */

#include "definitions_core.h"

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

/// \brief Destructor for a \ref Operator.
void destructor_Operator
	(const struct Operator* op ///< Standard.
	);

/// \brief Destructor for a \ref mutable_Operator.
void destructor_mutable_Operator
	(struct mutable_Operator* op ///< Standard.
	);

/// \brief Set the input arrays to hold the number of rows and columns of the tensor-product operators.
void set_ops_tp_n_rows_cols
	(int n_rows_sub[DMAX],                          ///< The number of rows of each sub-operator.
	 int n_cols_sub[DMAX],                          ///< The number of columns of each sub-operator.
	 const struct const_Multiarray_Matrix_d* ops_tp ///< The tensor-product operators.
	);

#endif // DPG__operator_h__INCLUDED
