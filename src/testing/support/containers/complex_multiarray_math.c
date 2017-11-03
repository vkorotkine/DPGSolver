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

#include "complex_multiarray_math.h"

#include <assert.h>

#include "complex_matrix.h"
#include "test_support_complex_multiarray.h"
#include "matrix.h"
#include "multiarray.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void transpose_Multiarray_c (struct Multiarray_c* a, const bool mem_only)
{
	const int order = a->order;

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = (order == 1 ? 1 : a->extents[1]);

assert(a->order <= 2);
// Need to do this block-wise for order > 2. Make sure to reset layout,ext_0,ext_1 before each transpose_Matrix call.

	struct Matrix_c* a_M = constructor_move_Matrix_c_c(a->layout,ext_0,ext_1,false,a->data); // destructed

	transpose_Matrix_c(a_M,mem_only);
	a->layout = a_M->layout;
	a->extents[0] = a_M->ext_0;
	a->extents[1] = a_M->ext_1;

	destructor_Matrix_c(a_M);
}

void mm_NNC_Multiarray_c
	(const double alpha, const double beta, const struct const_Matrix_d*const a,
	 const struct const_Multiarray_c*const b, struct Multiarray_c*const c)
{
	const char layout = 'C';
	assert(c->layout == layout);
	assert(b->layout == layout);

	const int order = b->order;
	assert(b->order == c->order);

	const ptrdiff_t ext_0_b = b->extents[0],
	                ext_0_c = c->extents[0];
	assert(a->ext_1 == ext_0_b);
	assert(a->ext_0 == ext_0_c);
	for (int i = 1; i < order; ++i)
		assert(b->extents[i] == c->extents[i]);

	const ptrdiff_t ext_1 = compute_size(order,b->extents)/ext_0_b;

	const struct const_Matrix_c* b_M =
		constructor_move_const_Matrix_c_c(b->layout,ext_0_b,ext_1,false,b->data); // destructed
	struct Matrix_c* c_M = constructor_move_Matrix_c_c(c->layout,ext_0_c,ext_1,false,c->data); // destructed

	mm_c('N','N',alpha,beta,a,b_M,c_M);

	destructor_const_Matrix_c(b_M);
	destructor_Matrix_c(c_M);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
