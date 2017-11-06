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
#include "complex_multiarray.h"
#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

/// \brief `complex` version of \ref permute_Multiarray_d.
void permute_Multiarray_c
	(struct Multiarray_c* a, ///< See brief.
	 const ptrdiff_t* p,     ///< See brief.
	 const char perm_layout  ///< See brief.
	);

/// \brief `complex` version of \ref reinterpret_Multiarray_as_Matrix_d.
void reinterpret_Multiarray_as_Matrix_c
	(struct Multiarray_c* a, ///< See brief.
	 struct Matrix_c* a_M,   ///< See brief.
	 const ptrdiff_t ext_0,  ///< See brief.
	 const ptrdiff_t ext_1   ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void transpose_Multiarray_c (struct Multiarray_c* a, const bool mem_only)
{
	const int order = a->order;

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = (order == 1 ? 1 : a->extents[1]);

	assert(a->order <= 2);

	struct Matrix_c* a_M = constructor_move_Matrix_c_c(a->layout,ext_0,ext_1,false,a->data); // destructed

	transpose_Matrix_c(a_M,mem_only);
	a->layout = a_M->layout;
	a->extents[0] = a_M->ext_0;
	a->extents[1] = a_M->ext_1;

	destructor_Matrix_c(a_M);
}

void permute_Multiarray_c_V (struct Multiarray_c* a, const struct const_Vector_i* p_V, const char perm_layout)
{
	const ptrdiff_t ext_0 = a->extents[0];
	assert(p_V->ext_0 == ext_0);

	ptrdiff_t p[ext_0];
	for (int i = 0; i < ext_0; ++i)
		p[i] = p_V->data[i];

	permute_Multiarray_c(a,p,perm_layout);
}

void scale_Multiarray_c (struct Multiarray_c* a, const double complex val)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] *= val;
}

void scale_Multiarray_c_by_Vector_d
	(const char side, const double alpha, struct Multiarray_c*const a, const struct const_Vector_d*const b,
	 const bool invert_diag)
{
	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = compute_size(a->order,a->extents)/ext_0;
	struct Matrix_c a_M;
	reinterpret_Multiarray_as_Matrix_c(a,&a_M,ext_0,ext_1);
	scale_Matrix_c_by_Vector_d(side,alpha,&a_M,b,invert_diag);
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

void permute_Multiarray_c (struct Multiarray_c* a, const ptrdiff_t* p, const char perm_layout)
{
	if (p == NULL)
		return;

	assert(a->order > 0);

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = compute_size(a->order,a->extents)/ext_0;

	struct Matrix_c a_M;
	reinterpret_Multiarray_as_Matrix_c(a,&a_M,ext_0,ext_1);

	if (perm_layout != a->layout)
		transpose_Matrix_c(&a_M,true);
	permute_Matrix_c(&a_M,p);
	if (perm_layout != a->layout)
		transpose_Matrix_c(&a_M,true);
}

void reinterpret_Multiarray_as_Matrix_c
	(struct Multiarray_c* a, struct Matrix_c* a_M, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	assert(compute_size(a->order,a->extents) == ext_0*ext_1);
	a_M->layout    = a->layout;
	a_M->ext_0     = ext_0;
	a_M->ext_1     = ext_1;
	a_M->owns_data = false;
	a_M->data      = a->data;
}
