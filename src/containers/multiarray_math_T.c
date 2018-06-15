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

#include <assert.h>

#include "macros.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "def_templates_math_functions.h"

// Static function declarations ************************************************************************************* //

/// \brief `mutable` version of \ref reinterpret_const_Multiarray_as_Matrix_T.
static void reinterpret_Multiarray_as_Matrix_T
	(struct Multiarray_T* a, ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_T.
	 struct Matrix_T* a_M,   ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_T.
	 const ptrdiff_t ext_0,  ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_T.
	 const ptrdiff_t ext_1   ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_T.
	);

/// \brief `mutable` version of \ref reinterpret_const_Matrix_as_Multiarray_T.
static void reinterpret_Matrix_as_Multiarray_T
	(struct Matrix_T* a_M,   ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_T.
	 struct Multiarray_T* a, ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_T.
	 const int order,        ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_T.
	 ptrdiff_t* extents      ///< Defined for \ref reinterpret_const_Matrix_as_Multiarray_T.
	);

// Interface functions ********************************************************************************************** //

void transpose_Multiarray_T (struct Multiarray_T* a, const bool mem_only)
{
	const int order = a->order;

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = (order == 1 ? 1 : a->extents[1]);

assert(a->order <= 2);
// Need to do this block-wise for order > 2. Make sure to reset layout,ext_0,ext_1 before each transpose_Matrix call.

	struct Matrix_T* a_M = constructor_move_Matrix_T_T(a->layout,ext_0,ext_1,false,a->data); // destructed

	transpose_Matrix_T(a_M,mem_only);
	a->layout = a_M->layout;
	a->extents[0] = a_M->ext_0;
	a->extents[1] = a_M->ext_1;

	destructor_Matrix_T(a_M);
}

void scale_Multiarray_T (struct Multiarray_T* a, const Type val)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] *= val;
}

void normalize_Multiarray_T
	(struct Multiarray_T* a, const char*const norm_type, const bool store_norms, struct Multiarray_T* a_norms)
{
	assert(a->layout == 'R'); // Can be made flexible if necessary.
	assert(a->order == 2);

	const ptrdiff_t n_vals    = a->extents[0],
	                n_entries = a->extents[1];

	if (store_norms) {
		assert(a_norms->order == 1);
		resize_Multiarray_T(a_norms,1,&n_vals);
		Type* norm_data = a_norms->data;
		for (ptrdiff_t i = 0; i < n_vals; ++i) {
			Type* a_row = get_row_Multiarray_T(i,a);
			norm_data[i] = norm_T(n_entries,a_row,norm_type);
			for (ptrdiff_t j = 0; j < n_entries; ++j)
				a_row[j] /= norm_data[i];
		}
	} else {
		Type norm_data = 0.0;
		for (ptrdiff_t i = 0; i < n_vals; ++i) {
			Type* a_row = get_row_Multiarray_T(i,a);
			norm_data = norm_T(n_entries,a_row,norm_type);
			for (ptrdiff_t j = 0; j < n_entries; ++j)
				a_row[j] /= norm_data;
		}
	}
}

void permute_Multiarray_T (struct Multiarray_T* a, const ptrdiff_t* p, const char perm_layout)
{
	if (p == NULL)
		return;

	assert(a->order > 0);

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = compute_size(a->order,a->extents)/ext_0;

	struct Matrix_T a_M;
	reinterpret_Multiarray_as_Matrix_T(a,&a_M,ext_0,ext_1);

	if (perm_layout != a->layout)
		transpose_Matrix_T(&a_M,true);
	permute_Matrix_T(&a_M,p);
	if (perm_layout != a->layout)
		transpose_Matrix_T(&a_M,true);
}

void permute_Multiarray_T_V (struct Multiarray_T* a, const struct const_Vector_i* p_V, const char perm_layout)
{
	const ptrdiff_t ext_0 = a->extents[0];
	assert(p_V->ext_0 == ext_0);

	ptrdiff_t p[ext_0];
	for (int i = 0; i < ext_0; ++i)
		p[i] = p_V->data[i];

	permute_Multiarray_T(a,p,perm_layout);
}

void scale_Multiarray_T_by_Vector_R
	(const char side, const Real alpha, struct Multiarray_T*const a, const struct const_Vector_R*const b,
	 const bool invert_diag)
{
	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = compute_size(a->order,a->extents)/ext_0;
	struct Matrix_T a_M;
	reinterpret_Multiarray_as_Matrix_T(a,&a_M,ext_0,ext_1);
	scale_Matrix_T_by_Vector_R(side,alpha,&a_M,b,invert_diag);
}

void scale_Multiarray_by_Vector_T
	(const char side, const Real alpha, struct Multiarray_T*const a, const struct const_Vector_T*const b,
	 const bool invert_diag)
{
	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = compute_size(a->order,a->extents)/ext_0;
	struct Matrix_T a_M;
	reinterpret_Multiarray_as_Matrix_T(a,&a_M,ext_0,ext_1);
	scale_Matrix_by_Vector_T(side,alpha,&a_M,b,invert_diag);
}

void add_in_place_Multiarray_T (const Real alpha, struct Multiarray_T*const a, const struct const_Multiarray_T* b)
{
	assert(check_equal_order_extents(a->order,b->order,a->extents,b->extents));
	assert(a->layout == b->layout);

	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] += alpha*(b->data[i]);
}

void multiply_in_place_Multiarray_T
	(const Type alpha, struct Multiarray_T*const a, const struct const_Multiarray_T*const b)
{
	assert(check_equal_order_extents(a->order,b->order,a->extents,b->extents));
	assert(a->layout == b->layout);

	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] *= alpha*(b->data[i]);
}

void multiply_in_place_Multiarray_TR
	(const Type alpha, struct Multiarray_T*const a, const struct const_Multiarray_R*const b)
{
	assert(check_equal_order_extents(a->order,b->order,a->extents,b->extents));
	assert(a->layout == b->layout);

	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] *= alpha*(b->data[i]);
}

void subtract_in_place_Multiarray_T (struct Multiarray_T* a, const struct const_Multiarray_T* b)
{
	assert(compute_size(a->order,a->extents) == compute_size(b->order,b->extents));

	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] -= b->data[i];
}

void mm_NNC_Multiarray_T
	(const Real alpha, const Real beta, const struct const_Matrix_R*const a,
	 const struct const_Multiarray_T*const b, struct Multiarray_T*const c)
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

	const struct const_Matrix_T* b_M =
		constructor_move_const_Matrix_T_T(b->layout,ext_0_b,ext_1,false,b->data); // destructed
	struct Matrix_T* c_M = constructor_move_Matrix_T_T(c->layout,ext_0_c,ext_1,false,c->data); // destructed

	mm_RTT('N','N',alpha,beta,a,b_M,c_M);

	destructor_const_Matrix_T(b_M);
	destructor_Matrix_T(c_M);
}

void mm_NNC_Multiarray_TTT
	(const Real alpha, const Real beta, const struct const_Matrix_T*const a,
	 const struct const_Multiarray_T*const b, struct Multiarray_T*const c)
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

	const struct const_Matrix_T* b_M =
		constructor_move_const_Matrix_T_T(b->layout,ext_0_b,ext_1,false,b->data); // destructed
	struct Matrix_T* c_M = constructor_move_Matrix_T_T(c->layout,ext_0_c,ext_1,false,c->data); // destructed

	mm_T('N','N',alpha,beta,a,b_M,c_M);

	destructor_const_Matrix_T(b_M);
	destructor_Matrix_T(c_M);
}

void mm_NN1C_Multiarray_T
	(const struct const_Matrix_R*const a, const struct const_Multiarray_T*const b, struct Multiarray_T*const c)
{
	mm_NNC_Multiarray_T(1.0,0.0,a,b,c);
}

void mm_NN1C_overwrite_Multiarray_T (const struct const_Matrix_R*const a, struct Multiarray_T** b)
{
	struct Multiarray_T* c = constructor_mm_NN1C_Multiarray_T(a,(struct const_Multiarray_T*)*b); // keep
	destructor_Multiarray_T(*b);
	*b = c;
}

void reinterpret_const_Multiarray_as_Matrix_T
	(const struct const_Multiarray_T* a, const struct const_Matrix_T* a_M, const ptrdiff_t ext_0,
	 const ptrdiff_t ext_1)
{
	reinterpret_Multiarray_as_Matrix_T((struct Multiarray_T*)a,(struct Matrix_T*)a_M,ext_0,ext_1);
}

void reinterpret_const_Matrix_as_Multiarray_T
	(const struct const_Matrix_T* a_M, const struct const_Multiarray_T* a, const int order, const ptrdiff_t* extents)
{
	reinterpret_Matrix_as_Multiarray_T((struct Matrix_T*)a_M,(struct Multiarray_T*)a,order,(ptrdiff_t*)extents);
}

void update_layout_Multiarray_Matrix_T (struct Multiarray_Matrix_T* a, const char layout_o)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (int i = 0; i < size; ++i) {
		if (a->data[i]->layout != layout_o)
			transpose_Matrix_T(a->data[i],true);
	}
}

Type norm_Multiarray_T (const struct Multiarray_T*const src, const char*const norm_type)
{
	const ptrdiff_t n_entries = compute_size(src->order,src->extents);
	return norm_T(n_entries,src->data,norm_type);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void reinterpret_Multiarray_as_Matrix_T
	(struct Multiarray_T* a, struct Matrix_T* a_M, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	assert(compute_size(a->order,a->extents) == ext_0*ext_1);
	a_M->layout    = a->layout;
	a_M->ext_0     = ext_0;
	a_M->ext_1     = ext_1;
	a_M->owns_data = false;
	a_M->data      = a->data;
}

static void reinterpret_Matrix_as_Multiarray_T
	(struct Matrix_T* a_M, struct Multiarray_T* a, const int order, ptrdiff_t* extents)
{
	assert(compute_size(order,extents) == ((a_M->ext_0)*(a_M->ext_1)));

	a->layout    = a_M->layout;
	a->order     = order;
	a->extents   = extents;
	a->owns_data = a_M->owns_data;
	a->data      = a_M->data;
}

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
#include "undef_templates_math_functions.h"
