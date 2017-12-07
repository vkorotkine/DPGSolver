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
#include "definitions_tol.h"


#include "def_templates_operators.h"

#include "def_templates_multiarray.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Constructors ***************************************************************************************************** //

struct Multiarray_T* constructor_mm_NN1_Operator_Multiarray_T
	(const struct Operator* op, const struct Multiarray_T* b, const char layout_c, const char op_format,
	 const int order_sub_ma, const ptrdiff_t* sub_inds_b)
{
	ptrdiff_t* extents = compute_extents_mm_MMa(op->op_std->ext_0,b->order,b->extents); // keep

	struct Multiarray_T* c = constructor_empty_Multiarray_T_dyn_extents(layout_c,order_sub_ma,extents); // returned

	mm_NN1_Operator_Multiarray_T(op,(struct const_Multiarray_T*)b,c,op_format,order_sub_ma,sub_inds_b,NULL);

	return c;
}

const struct const_Multiarray_T* constructor_mm_NN1_Operator_const_Multiarray_T
	(const struct Operator* op, const struct const_Multiarray_T* b, const char layout_c, const char op_format,
	 const int order_sub_ma, const ptrdiff_t* sub_inds_b)
{
	return (const struct const_Multiarray_T*) constructor_mm_NN1_Operator_Multiarray_T(
		op,(struct Multiarray_T*)b,layout_c,op_format,order_sub_ma,sub_inds_b);
}

// General functions ************************************************************************************************ //

void set_ops_tp_n_rows_cols_T
	(int n_rows_sub[DMAX], int n_cols_sub[DMAX], const struct const_Multiarray_Matrix_T* ops_tp)
{
	const ptrdiff_t size = compute_size(ops_tp->order,ops_tp->extents);
	assert(2 <= size);
	assert(size <= DMAX);

	for (int i = 0; i < DMAX; ++i) {
		if (i < size && ops_tp->data[i]) {
			n_rows_sub[i] = (int)ops_tp->data[i]->ext_0;
			n_cols_sub[i] = (int)ops_tp->data[i]->ext_1;
		} else {
			n_rows_sub[i] = 1;
			n_cols_sub[i] = 1;
		}
	}
}

// Math functions *************************************************************************************************** //

void mm_NNC_Operator_Multiarray_T
	(const Real alpha, const Real beta, const struct Operator* op, const struct const_Multiarray_T* b,
	 struct Multiarray_T* c, const char op_format, const int order_sub_ma, const ptrdiff_t* sub_inds_b,
	 const ptrdiff_t* sub_inds_c)
{
	assert(b->layout == 'C');
	assert(c->layout == 'C');

	const ptrdiff_t* extents_b = b->extents;
	ptrdiff_t* extents_c       = c->extents;
	assert(op->op_std->ext_1 == extents_b[0]);
	assert(op->op_std->ext_0 == extents_c[0]);

	for (int i = 1; i < order_sub_ma; ++i)
		assert(extents_b[i] == extents_c[i]);

	const struct const_Multiarray_T* b_op = NULL;
	struct Multiarray_T* c_op = NULL;

	const int order_b = b->order;
	if (order_sub_ma == order_b) {
		b_op = b;
	} else {
		assert(1 <= order_sub_ma);
		assert(order_sub_ma <= order_b);

		const ptrdiff_t ind_sub_b = compute_index_sub_container(order_b,order_sub_ma,extents_b,sub_inds_b);
		b_op = constructor_move_const_Multiarray_T_T(
			b->layout,order_sub_ma,(ptrdiff_t*)extents_b,false,&b->data[ind_sub_b]); // destructed
	}

	const int order_c = c->order;
	if (order_sub_ma == order_c) {
		c_op = c;
	} else {
		assert(1 <= order_sub_ma);
		assert(order_sub_ma <= order_c);

		const ptrdiff_t ind_sub_c = compute_index_sub_container(order_c,order_sub_ma,extents_c,sub_inds_c);

		c_op = constructor_move_Multiarray_T_T(
			c->layout,order_sub_ma,extents_c,false,&c->data[ind_sub_c]); // destructed
	}

	switch (op_format) {
	case 'd': // fallthrough
	case 's': mm_NNC_Multiarray_T(alpha,beta,op->op_std,b_op,c_op);    break;
//	case 't': mm_tp_NNC_Multiarray_T(alpha,beta,op->ops_tp,b_op,c_op); break;
	case 'c':
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",op_format);
		break;
	}

	if (order_sub_ma != order_b)
		destructor_const_Multiarray_T(b_op);
	if (order_sub_ma != order_c)
		destructor_Multiarray_T(c_op);
}

void mm_NN1C_Operator_Multiarray_T
	(const struct Operator* op, const struct const_Multiarray_T* b, struct Multiarray_T* c, const char op_format,
	 const int order_sub_ma, const ptrdiff_t* sub_inds_b, const ptrdiff_t* sub_inds_c)
{
	mm_NNC_Operator_Multiarray_T(1.0,0.0,op,b,c,op_format,order_sub_ma,sub_inds_b,sub_inds_c);
}

void mm_NN1_Operator_Multiarray_T
	(const struct Operator* op, const struct const_Multiarray_T* b, struct Multiarray_T* c, const char op_format,
	 const int order_sub_ma, const ptrdiff_t* sub_inds_b, const ptrdiff_t* sub_inds_c)
{
	const bool transpose_b = ( b->layout == 'C' ? false : true ),
	           transpose_c = ( c->layout == 'C' ? false : true );

	if (transpose_b)
		transpose_Multiarray_T((struct Multiarray_T*)b,true);
	if (transpose_c)
		swap_layout(&c->layout); // Data is about to be overwritten

	mm_NN1C_Operator_Multiarray_T(op,b,c,op_format,order_sub_ma,sub_inds_b,sub_inds_c);

	if (transpose_b)
		transpose_Multiarray_T((struct Multiarray_T*)b,true);
	if (transpose_c)
		transpose_Multiarray_T(c,true);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
