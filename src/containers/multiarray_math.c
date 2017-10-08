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

#include "multiarray_math.h"

#include <assert.h>

#include "macros.h"

#include "multiarray.h"
#include "matrix.h"
#include "math_functions.h"

// Static function declarations ************************************************************************************* //

/// \brief `mutable` version of \ref reinterpret_const_Multiarray_as_Matrix_d.
void reinterpret_Multiarray_as_Matrix_d
	(struct Multiarray_d* a, ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	 struct Matrix_d* a_M,   ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	 const ptrdiff_t ext_0,  ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	 const ptrdiff_t ext_1   ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	);

/// \brief `mutable` version of \ref reinterpret_const_Matrix_as_Multiarray_d.
void reinterpret_Matrix_as_Multiarray_d
	(struct Matrix_d* a_M,   ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_d.
	 struct Multiarray_d* a, ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_d.
	 const int order,        ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_d.
	 ptrdiff_t* extents      ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_d.
	);

// Interface functions ********************************************************************************************** //

void swap_layout_Multiarray_d (struct Multiarray_d* a)
{
	a->layout = ( a->layout == 'R' ? 'C' : 'R' );
}

void transpose_Multiarray_d (struct Multiarray_d* a, const bool mem_only)
{
	const int order = a->order;

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = (order == 1 ? 1 : a->extents[1]);

assert(a->order <= 2);
// Need to do this block-wise for order > 2. Make sure to reset layout,ext_0,ext_1 before each transpose_Matrix call.

	struct Matrix_d* a_M = constructor_move_Matrix_d_d(a->layout,ext_0,ext_1,false,a->data); // destructed

	transpose_Matrix_d(a_M,mem_only);
	a->layout = a_M->layout;
	a->extents[0] = a_M->ext_0;
	a->extents[1] = a_M->ext_1;

	destructor_Matrix_d(a_M);
}

void normalize_Multiarray_d
	(struct Multiarray_d* a, const char*const norm_type, const bool store_norms, struct Multiarray_d* a_norms)
{
	assert(a->layout == 'R'); // Can be made flexible if necessary.
	assert(a->order == 2);

	if (store_norms) {
		const ptrdiff_t n_vals    = a->extents[0],
		                n_entries = a->extents[1];

		assert(a_norms->order == 1);
		resize_Multiarray_d(a_norms,1,&n_vals);
		double* norm_data = a_norms->data;
		for (ptrdiff_t i = 0; i < n_vals; ++i) {
			double* a_row = get_row_Multiarray_d(i,a);
			norm_data[i] = norm_d(n_entries,a_row,norm_type);
			for (ptrdiff_t j = 0; j < n_entries; ++j)
				a_row[j] /= norm_data[i];
		}
	} else {
		EXIT_ADD_SUPPORT;
		double norm_data = 0.0;
		UNUSED(norm_data);
	}
}

void permute_Multiarray_d (struct Multiarray_d* a, const ptrdiff_t* p)
{
	if (p == NULL)
		return;

	assert(a->order > 0);

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = compute_size(a->order,a->extents)/ext_0;

	struct Matrix_d a_M;
	reinterpret_Multiarray_as_Matrix_d(a,&a_M,ext_0,ext_1);

	permute_Matrix_d (&a_M,p);
}

void mm_NN1C_Multiarray_d
	(const struct const_Matrix_d*const a, const struct const_Multiarray_d*const b, struct Multiarray_d*const c)
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

	const struct const_Matrix_d* b_M =
		constructor_move_const_Matrix_d_d(b->layout,ext_0_b,ext_1,false,b->data); // destructed
	struct Matrix_d* c_M = constructor_move_Matrix_d_d(c->layout,ext_0_c,ext_1,false,c->data); // destructed

	mm_d('N','N',1.0,0.0,a,b_M,c_M);

	destructor_const_Matrix_d(b_M);
	destructor_Matrix_d(c_M);
}

void reinterpret_const_Multiarray_as_Matrix_d
	(const struct const_Multiarray_d* a, const struct const_Matrix_d* a_M, const ptrdiff_t ext_0,
	 const ptrdiff_t ext_1)
{
	reinterpret_Multiarray_as_Matrix_d((struct Multiarray_d*)a,(struct Matrix_d*)a_M,ext_0,ext_1);
}

void reinterpret_const_Matrix_as_Multiarray_d
	(const struct const_Matrix_d* a_M, const struct const_Multiarray_d* a, const int order, const ptrdiff_t* extents)
{
	reinterpret_Matrix_as_Multiarray_d((struct Matrix_d*)a_M,(struct Multiarray_d*)a,order,(ptrdiff_t*)extents);
}

ptrdiff_t* compute_extents_mm_MMa (const ptrdiff_t ext_0, const int order, const ptrdiff_t* extents_i)
{
	ptrdiff_t* extents = malloc(order * sizeof *extents); // returned

	extents[0] = ext_0;
	for (int i = 1; i < order; ++i)
		extents[i] = extents_i[i];

	return extents;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

void reinterpret_Multiarray_as_Matrix_d
	(struct Multiarray_d* a, struct Matrix_d* a_M, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	assert(compute_size(a->order,a->extents) == ext_0*ext_1);
	a_M->layout    = a->layout;
	a_M->ext_0     = ext_0;
	a_M->ext_1     = ext_1;
	a_M->owns_data = false;
	a_M->data      = a->data;
}

void reinterpret_Matrix_as_Multiarray_d
	(struct Matrix_d* a_M, struct Multiarray_d* a, const int order, ptrdiff_t* extents)
{
	assert(compute_size(order,extents) == ((a_M->ext_0)*(a_M->ext_1)));

	a->layout    = a_M->layout;
	a->order     = order;
	a->extents   = extents;
	a->owns_data = a_M->owns_data;
	a->data      = a_M->data;
}
