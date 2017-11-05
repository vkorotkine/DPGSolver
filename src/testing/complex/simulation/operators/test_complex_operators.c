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

#include "test_complex_operators.h"

#include <assert.h>

#include "macros.h"

#include "complex_multiarray.h"
#include "matrix.h"
#include "multiarray.h"

#include "multiarray_operator.h"
#include "operator.h"

// Static function declarations ************************************************************************************* //

/// \brief `complex` version of \ref mm_NN1_Operator_Multiarray_d.
static void mm_NN1_Operator_Multiarray_c
	(const struct Operator* op,          ///< Defined for \ref mm_NN1_Operator_Multiarray_d.
	 const struct const_Multiarray_c* b, ///< Defined for \ref mm_NN1_Operator_Multiarray_d.
	 struct Multiarray_c* c,             ///< Defined for \ref mm_NN1_Operator_Multiarray_d.
	 const char op_format,               ///< Defined for \ref mm_NN1_Operator_Multiarray_d.
	 const int order_sub_ma,             ///< Defined for \ref mm_NN1_Operator_Multiarray_d.
	 const ptrdiff_t* sub_inds_b,        ///< Defined for \ref mm_NN1_Operator_Multiarray_d.
	 const ptrdiff_t* sub_inds_c         ///< Defined for \ref mm_NN1_Operator_Multiarray_d.
	);

// Interface functions ********************************************************************************************** //
// Constructors ***************************************************************************************************** //

/** \brief `complex` version of \ref constructor_mm_NN1_Operator_const_Multiarray_d.
 *  \return See brief. */
const struct const_Multiarray_c* constructor_mm_NN1_Operator_const_Multiarray_c
	(const struct Operator* op, const struct const_Multiarray_c* b, const char layout_c, const char op_format,
	 const int order_sub_ma, const ptrdiff_t* sub_inds_b)
{
	ptrdiff_t* extents = compute_extents_mm_MMa(op->op_std->ext_0,b->order,b->extents); // keep

	struct Multiarray_c* c = constructor_empty_Multiarray_c_dyn_extents(layout_c,order_sub_ma,extents); // returned

	mm_NN1_Operator_Multiarray_c(op,b,c,op_format,order_sub_ma,sub_inds_b,NULL);

	return (const struct const_Multiarray_c*) c;
}

void mm_NNC_Operator_Multiarray_c
	(const double alpha, const double beta, const struct Operator* op, const struct const_Multiarray_c* b,
	 struct Multiarray_c* c, const char op_format, const int order_sub_ma, const ptrdiff_t* sub_inds_b,
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

	const struct const_Multiarray_c* b_op = NULL;
	struct Multiarray_c* c_op = NULL;

	const int order_b = b->order;
	if (order_sub_ma == order_b) {
		b_op = b;
	} else {
		assert(1 <= order_sub_ma);
		assert(order_sub_ma <= order_b);

		const ptrdiff_t ind_sub_b = compute_index_sub_container(order_b,order_sub_ma,extents_b,sub_inds_b);
		b_op = constructor_move_const_Multiarray_c_c(
			b->layout,order_sub_ma,(ptrdiff_t*)extents_b,false,&b->data[ind_sub_b]); // destructed
	}

	const int order_c = c->order;
	if (order_sub_ma == order_c) {
		c_op = c;
	} else {
		assert(1 <= order_sub_ma);
		assert(order_sub_ma <= order_c);

		const ptrdiff_t ind_sub_c = compute_index_sub_container(order_c,order_sub_ma,extents_c,sub_inds_c);

		c_op = constructor_move_Multiarray_c_c(
			c->layout,order_sub_ma,extents_c,false,&c->data[ind_sub_c]); // destructed
	}

	assert(op_format == 'd');
	mm_NNC_Multiarray_c(alpha,beta,op->op_std,b_op,c_op);

	if (order_sub_ma != order_b)
		destructor_const_Multiarray_c(b_op);
	if (order_sub_ma != order_c)
		destructor_Multiarray_c(c_op);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief `complex` version of \ref mm_NN1C_Operator_Multiarray_d.
void mm_NN1C_Operator_Multiarray_c
	(const struct Operator* op,          ///< Defined for \ref mm_NN1C_Operator_Multiarray_d.
	 const struct const_Multiarray_c* b, ///< Defined for \ref mm_NN1C_Operator_Multiarray_d.
	 struct Multiarray_c* c,             ///< Defined for \ref mm_NN1C_Operator_Multiarray_d.
	 const char op_format,               ///< Defined for \ref mm_NN1C_Operator_Multiarray_d.
	 const int order_sub_ma,             ///< Defined for \ref mm_NN1C_Operator_Multiarray_d.
	 const ptrdiff_t* sub_inds_b,        ///< Defined for \ref mm_NN1C_Operator_Multiarray_d.
	 const ptrdiff_t* sub_inds_c         ///< Defined for \ref mm_NN1C_Operator_Multiarray_d.
	);

static void mm_NN1_Operator_Multiarray_c
	(const struct Operator* op, const struct const_Multiarray_c* b, struct Multiarray_c* c, const char op_format,
	 const int order_sub_ma, const ptrdiff_t* sub_inds_b, const ptrdiff_t* sub_inds_c)
{
	const bool transpose_b = ( b->layout == 'C' ? false : true ),
	           transpose_c = ( c->layout == 'C' ? false : true );

	if (transpose_b)
		transpose_Multiarray_c((struct Multiarray_c*)b,true);
	if (transpose_c)
		swap_layout(&c->layout); // Data is about to be overwritten

	mm_NN1C_Operator_Multiarray_c(op,b,c,op_format,order_sub_ma,sub_inds_b,sub_inds_c);

	if (transpose_b)
		transpose_Multiarray_c((struct Multiarray_c*)b,true);
	if (transpose_c)
		transpose_Multiarray_c(c,true);
}

// Level 1 ********************************************************************************************************** //

void mm_NN1C_Operator_Multiarray_c
	(const struct Operator* op, const struct const_Multiarray_c* b, struct Multiarray_c* c, const char op_format,
	 const int order_sub_ma, const ptrdiff_t* sub_inds_b, const ptrdiff_t* sub_inds_c)
{
	mm_NNC_Operator_Multiarray_c(1.0,0.0,op,b,c,op_format,order_sub_ma,sub_inds_b,sub_inds_c);
}
