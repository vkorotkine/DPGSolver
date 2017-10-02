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
 */

#include "operator.h"
#include <assert.h>

#include "macros.h"
#include "definitions_core.h"

#include "multiarray.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Destructors ****************************************************************************************************** //

void destructor_Operator (const struct Operator* op)
{
	destructor_mutable_Operator((struct mutable_Operator*)op);
}

void destructor_mutable_Operator (struct mutable_Operator* op)
{
	if (op->op_std)
		destructor_Matrix_d(op->op_std);
	if (op->ops_tp) {
		assert(op->ops_tp->owns_data == false);
		destructor_Multiarray_Matrix_d(op->ops_tp);
	}
	if (op->op_csr)
		EXIT_ADD_SUPPORT;
//		destructor_Matrix_CSR_d(op->op_csr);
	free(op);
}

void set_ops_tp_n_rows_cols
	(int n_rows_sub[DMAX], int n_cols_sub[DMAX], const struct const_Multiarray_Matrix_d* ops_tp)
{
	const ptrdiff_t size = compute_size(ops_tp->order,ops_tp->extents);
	assert(2 <= size);
	assert(size <= DMAX);

	for (int i = 0; i < DMAX; ++i) {
		if (i < size && ops_tp->data[i]) {
			n_rows_sub[i] = ops_tp->data[i]->ext_0;
			n_cols_sub[i] = ops_tp->data[i]->ext_1;
		} else {
			n_rows_sub[i] = 1;
			n_cols_sub[i] = 1;
		}
	}
}


// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
