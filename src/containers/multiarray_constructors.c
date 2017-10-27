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

#include "multiarray_constructors.h"

#include <stdlib.h>
#include <assert.h>

#include "macros.h"
#include "definitions_core.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "const_cast.h"
#include "operator.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Helper functions ************************************************************************************************* //

ptrdiff_t* allocate_and_set_extents (const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t* extents = malloc(order * sizeof *extents); // returned

	for (ptrdiff_t i = 0; i < order; i++)
		extents[i] = extents_i[i];

	return extents;
}

// Default constructors ********************************************************************************************* //

struct Multiarray_d* constructor_default_Multiarray_d ()
{
	return constructor_move_Multiarray_d_dyn_extents('C',0,NULL,true,NULL);
}

const struct const_Multiarray_d* constructor_default_const_Multiarray_d ()
{
	return (const struct const_Multiarray_d*)constructor_default_Multiarray_d();
}

struct Multiarray_Matrix_d* constructor_default_Multiarray_Matrix_d ()
{
	return constructor_move_Multiarray_Matrix_d_dyn_extents(0,NULL,true,NULL);
}

const struct const_Multiarray_Matrix_d* constructor_default_const_Multiarray_Matrix_d ()
{
	return (const struct const_Multiarray_Matrix_d*) constructor_default_Multiarray_Matrix_d();
}

// Empty constructors *********************************************************************************************** //

struct Multiarray_d* constructor_empty_Multiarray_d
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_empty_Multiarray_d_dyn_extents(layout,order,extents);
}

struct Multiarray_d* constructor_empty_Multiarray_d_dyn_extents
	(const char layout, const int order, const ptrdiff_t*const extents)
{
	double* data = malloc(compute_size(order,extents) * sizeof *data); // keep

	return constructor_move_Multiarray_d_dyn_extents(layout,order,(ptrdiff_t*)extents,true,data);
}

struct Multiarray_Vector_i* constructor_empty_Multiarray_Vector_i
	(const bool alloc_V, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep

	struct Vector_i** data = NULL;
	if (alloc_V)
		data = constructor_default_Vector_i_2(compute_size(order,extents)); // keep
	else
		data = calloc(1,compute_size(order,extents) * sizeof *data); // keep

	return constructor_move_Multiarray_Vector_i_dyn_extents(order,extents,true,data);
}

const struct const_Multiarray_Vector_i* constructor_empty_const_Multiarray_Vector_i
	(const bool alloc_V, const int order, const ptrdiff_t*const extents_i)
{
	return (const struct const_Multiarray_Vector_i*) constructor_empty_Multiarray_Vector_i(alloc_V,order,extents_i);
}

const struct const_Multiarray_Vector_i* constructor_empty_const_Multiarray_Vector_i_V
	(const bool alloc_V, const struct const_Vector_i*const extents_i_V)
{
	const ptrdiff_t order = extents_i_V->ext_0;

	ptrdiff_t extents_i[order];
	for (ptrdiff_t i = 0; i < order; ++i)
		extents_i[i] = extents_i_V->data[i];

	return constructor_empty_const_Multiarray_Vector_i(alloc_V,order,extents_i);
}

struct Multiarray_Vector_d* constructor_empty_Multiarray_Vector_d
	(const bool alloc_V, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep

	struct Vector_d** data = NULL;
	if (alloc_V)
		EXIT_ADD_SUPPORT;
	else
		data = calloc(1,compute_size(order,extents) * sizeof *data); // keep

	return constructor_move_Multiarray_Vector_d_dyn_extents(order,extents,true,data);
}

const struct const_Multiarray_Vector_d* constructor_empty_const_Multiarray_Vector_d
	(const bool alloc_V, const int order, const ptrdiff_t*const extents_i)
{
	return (const struct const_Multiarray_Vector_d*) constructor_empty_Multiarray_Vector_d(alloc_V,order,extents_i);
}

const struct const_Multiarray_Vector_d* constructor_empty_const_Multiarray_Vector_d_V
	(const bool alloc_V, const struct const_Vector_i*const extents_i_V)
{
	const ptrdiff_t order = extents_i_V->ext_0;

	ptrdiff_t extents_i[order];
	for (ptrdiff_t i = 0; i < order; ++i)
		extents_i[i] = extents_i_V->data[i];

	return constructor_empty_const_Multiarray_Vector_d(alloc_V,order,extents_i);
}

struct Multiarray_Matrix_d* constructor_empty_Multiarray_Matrix_d
	(const bool alloc_M, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep

	struct Matrix_d** data = NULL;
	if (alloc_M)
		EXIT_ADD_SUPPORT;
	else
		data = calloc(compute_size(order,extents) , sizeof *data); // keep

	return constructor_move_Multiarray_Matrix_d_dyn_extents(order,extents,true,data);
}

const struct const_Multiarray_Matrix_d* constructor_empty_const_Multiarray_Matrix_d
	(const bool alloc_M, const int order, const ptrdiff_t*const extents_i)
{
	return (const struct const_Multiarray_Matrix_d*) constructor_empty_Multiarray_Matrix_d(alloc_M,order,extents_i);
}


struct Multiarray_Matrix_d* constructor_empty_Multiarray_Matrix_d_V
	(const bool alloc_M, const struct const_Vector_i*const extents_i_V)
{
	const ptrdiff_t order = extents_i_V->ext_0;

	ptrdiff_t extents_i[order];
	for (ptrdiff_t i = 0; i < order; ++i)
		extents_i[i] = extents_i_V->data[i];

	return constructor_empty_Multiarray_Matrix_d(alloc_M,order,extents_i);
}

const struct const_Multiarray_Matrix_d* constructor_empty_const_Multiarray_Matrix_d_V
	(const bool alloc_M, const struct const_Vector_i*const extents_i_V)
{
	return (const struct const_Multiarray_Matrix_d*) constructor_empty_Multiarray_Matrix_d_V(alloc_M,extents_i_V);
}

// Zero constructors ************************************************************************************************ //

struct Multiarray_i* constructor_zero_Multiarray_i
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_zero_Multiarray_i_dyn_extents(layout,order,extents);
}

struct Multiarray_i* constructor_zero_Multiarray_i_dyn_extents
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	int* data = calloc(compute_size(order,extents_i) , sizeof *data); // keep
	return constructor_move_Multiarray_i_dyn_extents(layout,order,(ptrdiff_t*)extents_i,true,data);
}

struct Multiarray_d* constructor_zero_Multiarray_d
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_zero_Multiarray_d_dyn_extents(layout,order,extents);
}

struct Multiarray_d* constructor_zero_Multiarray_d_dyn_extents
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	double* data = calloc(compute_size(order,extents_i) , sizeof *data); // keep
	return constructor_move_Multiarray_d_dyn_extents(layout,order,(ptrdiff_t*)extents_i,true,data);
}

// Copy constructors ************************************************************************************************ //

struct Multiarray_Vector_i* constructor_copy_Multiarray_Vector_i_i
	(const int* data_V, const int*const ext_V, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // free

	struct Multiarray_Vector_i* dest = constructor_empty_Multiarray_Vector_i(true,order,extents); // returned
	free(extents);

	set_Multiarray_Vector_i_i(dest,data_V,ext_V);

	return dest;
}

struct Multiarray_d* constructor_copy_Multiarray_d (struct Multiarray_d* src)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);
	double* data = malloc(size * sizeof *data); // moved
	for (int i = 0; i < size; ++i)
		data[i] = src->data[i];

	return constructor_move_Multiarray_d_d(src->layout,src->order,src->extents,true,data);
}

const struct const_Multiarray_d* constructor_copy_const_Multiarray_d (const struct const_Multiarray_d*const src)
{
	return (const struct const_Multiarray_d*) constructor_copy_Multiarray_d((struct Multiarray_d*)src);
}

void const_constructor_copy_Multiarray_d
	(const struct const_Multiarray_d*const* dest, const struct const_Multiarray_d*const src)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);
	double* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] = src->data[i];

	struct Multiarray_d* dest_m = constructor_move_Multiarray_d_d(src->layout,src->order,src->extents,true,data);
	const_constructor_move_Multiarray_d(dest,dest_m);
}

// Move constructors ************************************************************************************************ //

struct Multiarray_d* constructor_move_Multiarray_d_d
	(const char layout, const int order, const ptrdiff_t*const extents_i, const bool owns_data, double*const data)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_move_Multiarray_d_dyn_extents(layout,order,extents,owns_data,data);
}

const struct const_Multiarray_d* constructor_move_const_Multiarray_d_d
	(const char layout, const int order, const ptrdiff_t*const extents_i, const bool owns_data,
	 const double*const data)
{
	return (const struct const_Multiarray_d*)
		constructor_move_Multiarray_d_d(layout,order,extents_i,owns_data,(double*)data);
}

struct Multiarray_i* constructor_move_Multiarray_i_dyn_extents
	(const char layout, const int order, ptrdiff_t*const extents, const bool owns_data, int*const data)
{
	struct Multiarray_i* dest = calloc(1,sizeof *dest); // returned

	dest->layout    = layout;
	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct Multiarray_d* constructor_move_Multiarray_d_dyn_extents
	(const char layout, const int order, ptrdiff_t*const extents, const bool owns_data, double*const data)
{
	struct Multiarray_d* dest = calloc(1,sizeof *dest); // returned

	dest->layout    = layout;
	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

const struct const_Multiarray_d* constructor_move_const_Multiarray_d_dyn_extents
	(const char layout, const int order, ptrdiff_t*const extents, const bool owns_data, const double*const data)
{
	return (const struct const_Multiarray_d*)
		constructor_move_Multiarray_d_dyn_extents(layout,order,extents,owns_data,(double*)data);
}

struct Multiarray_Vector_i* constructor_move_Multiarray_Vector_i_dyn_extents
	(const int order, ptrdiff_t*const extents, const bool owns_data, struct Vector_i**const data)
{
	struct Multiarray_Vector_i* dest = calloc(1,sizeof *dest); // returned

	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct Multiarray_Vector_d* constructor_move_Multiarray_Vector_d_dyn_extents
	(const int order, ptrdiff_t*const extents, const bool owns_data, struct Vector_d**const data)
{
	struct Multiarray_Vector_d* dest = calloc(1,sizeof *dest); // returned

	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct Multiarray_Matrix_d* constructor_move_Multiarray_Matrix_d_dyn_extents
	(const int order, ptrdiff_t*const extents, const bool owns_data, struct Matrix_d**const data)
{
	struct Multiarray_Matrix_d* dest = calloc(1,sizeof *dest); // returned

	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct Multiarray_d* constructor_move_Multiarray_d_Matrix_d (struct Matrix_d* src)
{
	src->owns_data = false;
	return constructor_move_Multiarray_d_d(src->layout,2,(ptrdiff_t[]){src->ext_0,src->ext_1},true,src->data);
}

const struct const_Multiarray_d* constructor_move_const_Multiarray_d_Matrix_d (const struct const_Matrix_d* src)
{
	return (const struct const_Multiarray_d*) constructor_move_Multiarray_d_Matrix_d((struct Matrix_d*)src);
}

void const_constructor_move_Multiarray_d (const struct const_Multiarray_d*const* dest, struct Multiarray_d* src)
{
	*(struct const_Multiarray_d**) dest = (struct const_Multiarray_d*) src;
}

void const_constructor_move_const_Multiarray_d
	(const struct const_Multiarray_d*const* dest, const struct const_Multiarray_d*const src)
{
	const_constructor_move_Multiarray_d(dest,(struct Multiarray_d*)src);
}

void const_constructor_move_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const* dest, struct Multiarray_Vector_i* src)
{
	*(struct const_Multiarray_Vector_i**) dest = (struct const_Multiarray_Vector_i*) src;
}

void const_constructor_move_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d*const* dest, struct Multiarray_Matrix_d* src)
{
	*(struct const_Multiarray_Matrix_d**) dest = (struct const_Multiarray_Matrix_d*) src;
}

// Special constructors ********************************************************************************************* //

const struct const_Multiarray_d* constructor_MaM1_V_const_Multiarray_d
	(const char layout, const char trans_a, const double alpha, const double beta,
	 const struct const_Multiarray_Matrix_d*const A, const struct const_Vector_d* b)
{
	assert(A->order == 1);

	const ptrdiff_t order = 2,
	                ext_0 = A->data[0]->ext_0,
	                ext_1 = compute_size(A->order,A->extents);

	struct Multiarray_d* dest = constructor_empty_Multiarray_d(layout,order,(ptrdiff_t[]){ext_0,ext_1}); // returned

	struct Vector_d dest_V;
	for (ptrdiff_t i = 0; i < ext_1; ++i) {
		const struct const_Matrix_d*const A_M = A->data[i];
		set_Vector_from_Multiarray_d(&dest_V,dest,&i);

		mv_d(trans_a,alpha,beta,A_M,b,&dest_V);
	}
	return (const struct const_Multiarray_d*) dest;
}

void set_Multiarray_Matrix_from_Multiarray_Matrix_d
	(struct Multiarray_Matrix_d* dest, struct Multiarray_Matrix_d* src, const int order_o,
	 const ptrdiff_t*const sub_indices)
{
	dest->owns_data = false;
	dest->order     = order_o;
	dest->extents   = src->extents;
	dest->data      = &src->data[compute_index_sub_container(src->order,dest->order,src->extents,sub_indices)];
}

void set_const_Multiarray_Matrix_from_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d* dest, const struct const_Multiarray_Matrix_d* src, const int order_o,
	 const ptrdiff_t*const sub_indices)
{
	set_Multiarray_Matrix_from_Multiarray_Matrix_d(
		(struct Multiarray_Matrix_d*)dest,(struct Multiarray_Matrix_d*)src,order_o,sub_indices);
}

struct Multiarray_d* constructor_mm_NN1C_Multiarray_d
	(const struct const_Matrix_d*const a, const struct const_Multiarray_d*const b)
{
	const char layout  = 'C';
	const int order    = b->order;
	ptrdiff_t* extents = compute_extents_mm_MMa(a->ext_0,b->order,b->extents); // keep

	struct Multiarray_d* c = constructor_empty_Multiarray_d_dyn_extents(layout,order,extents); // returned

	mm_NN1C_Multiarray_d(a,b,c);

	return c;
}

const struct const_Multiarray_d* constructor_mm_NN1C_const_Multiarray_d
	(const struct const_Matrix_d*const a, const struct const_Multiarray_d*const b)
{
	return (const struct const_Multiarray_d*) constructor_mm_NN1C_Multiarray_d(a,b);
}

/*const struct const_Multiarray_d* constructor_mm_NN1_const_Multiarray_d
	(const struct const_Matrix_d*const a, const struct const_Multiarray_d*const b, const char layout_c)
{
	ptrdiff_t* extents = compute_extents_mm_MMa(a->ext_0,b->order,b->extents); // keep

	struct Multiarray_d* c = constructor_empty_Multiarray_d_dyn_extents(layout_c,b->order,extents); // returned

	mm_NN1_Multiarray_d(a,b,c);

	return (const struct const_Multiarray_d*) c;
}*/

const struct const_Multiarray_d* constructor_mm_tp_NN1C_const_Multiarray_d
	(const struct const_Multiarray_Matrix_d* a_tp, const struct const_Multiarray_d* b)
{
	/**
	 *  This function works by applying the tensor-product sub-operators to the appropriately rearranged memory of the
	 *  input multiarray and its intermediate results.
	 *
	 *  In the interest of efficiency, memory is **always** rearranged using matrix transpose calls, which **require**
	 *  that the operators be applied to a column-major inputs. The basic idea is to reinterpret the matrix input as
	 *  having extents such that the desired leading dimension will be obtained after transposing.
	 *
	 *  \todo Add example documentation for this process.
	 */
	assert(a_tp->order == 1);

	const ptrdiff_t d_op = a_tp->extents[0];
	assert(1 < d_op);
	assert(d_op <= DMAX);

	assert(b->layout == 'C');

	int n_rows_sub[DMAX] = { 0, 0, 0, },
	    n_cols_sub[DMAX] = { 0, 0, 0, };
	set_ops_tp_n_rows_cols(n_rows_sub,n_cols_sub,a_tp);

	const int order_i          = b->order;
	const ptrdiff_t* extents_i = b->extents;
	assert(extents_i[0] == (n_cols_sub[0]*n_cols_sub[1]*n_cols_sub[2]));

	const ptrdiff_t ext_1_Ma = compute_size(order_i,extents_i)/extents_i[0];

	const struct const_Matrix_d* sub_op = NULL;

	// r-direction
	int ind_sub_op = 0;
	sub_op = a_tp->data[ind_sub_op];
	assert(sub_op != NULL);

	const struct const_Matrix_d* b_M = constructor_default_const_Matrix_d(); // destructed
	reinterpret_const_Multiarray_as_Matrix_d(b,b_M,n_cols_sub[0],n_cols_sub[1]*n_cols_sub[2]*ext_1_Ma);

	const struct const_Matrix_d* c_r = constructor_mm_NN1C_const_Matrix_d(sub_op,b_M); // destructed/moved
	destructor_const_Matrix_d(b_M);

	// s-direction
	++ind_sub_op;
	sub_op = a_tp->data[ind_sub_op];

	const struct const_Matrix_d* c_rs = NULL;
	if (sub_op) {
		transpose_Matrix_d((struct Matrix_d*)c_r,false);
		reinterpret_const_Matrix_d(c_r,n_cols_sub[1],n_rows_sub[0]*n_cols_sub[2]*ext_1_Ma);

		c_rs = constructor_mm_NN1C_const_Matrix_d(sub_op,c_r); // destructed/moved
		destructor_const_Matrix_d(c_r);

		reinterpret_const_Matrix_d(c_rs,n_rows_sub[1]*n_cols_sub[2]*ext_1_Ma,n_rows_sub[0]);
		transpose_Matrix_d((struct Matrix_d*)c_rs,false);
	} else {
		assert(d_op == 3);
		c_rs = c_r;
	}

	// t-direction
	const struct const_Matrix_d* c_rst = NULL;
	if (d_op > 2) {
		++ind_sub_op;
		sub_op = a_tp->data[ind_sub_op];

		reinterpret_const_Matrix_d(c_rs,n_rows_sub[0]*n_rows_sub[1],n_cols_sub[2]*ext_1_Ma);
		transpose_Matrix_d((struct Matrix_d*)c_rs,false);
		reinterpret_const_Matrix_d(c_rs,n_cols_sub[2],n_rows_sub[0]*n_rows_sub[1]*ext_1_Ma);

		c_rst = constructor_mm_NN1C_const_Matrix_d(sub_op,c_rs); // destructed/moved
		destructor_const_Matrix_d(c_rs);

		reinterpret_const_Matrix_d(c_rst,n_rows_sub[2]*ext_1_Ma,n_rows_sub[0]*n_rows_sub[1]);
		transpose_Matrix_d((struct Matrix_d*)c_rst,false);
	} else {
		assert(d_op == 2);
		c_rst = c_rs;
	}

	ptrdiff_t* extents_o = malloc(order_i * sizeof *extents_o); // keep
	extents_o[0] = n_rows_sub[0]*n_rows_sub[1]*n_rows_sub[2];
	for (int i = 1; i < order_i; ++i)
		extents_o[i] = extents_i[i];

	const struct const_Multiarray_d* c = constructor_default_const_Multiarray_d(); // returned
	reinterpret_const_Matrix_as_Multiarray_d(c_rst,c,order_i,extents_o);

	const_cast_bool(&c_rst->owns_data,false);
	destructor_const_Matrix_d(c_rst);

	return c;
}

// Destructors ****************************************************************************************************** //

void destructor_Multiarray_i (struct Multiarray_i* a)
{
	assert(a != NULL);

	free(a->extents);
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Multiarray_i (const struct const_Multiarray_i* a)
{
	destructor_Multiarray_i((struct Multiarray_i*)a);
}

void destructor_Multiarray_d (struct Multiarray_d* a)
{
	assert(a != NULL);

	free(a->extents);
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Multiarray_d (const struct const_Multiarray_d* a)
{
	destructor_Multiarray_d((struct Multiarray_d*)a);
}

void destructor_Multiarray_Vector_i (struct Multiarray_Vector_i* a)
{
	assert(a != NULL);
	assert(a->data != NULL);

	if (a->owns_data) {
		const ptrdiff_t size = compute_size(a->order,a->extents);
		for (ptrdiff_t i = 0; i < size; ++i) {
			// It is possible that not all vectors present in the Multiarray were set.
			if (a->data[i])
				destructor_Vector_i(a->data[i]);
		}
	}
	free(a->data);
	free(a->extents);
	free(a);
}

void destructor_const_Multiarray_Vector_i (const struct const_Multiarray_Vector_i* a)
{
	destructor_Multiarray_Vector_i((struct Multiarray_Vector_i*)a);
}

void destructor_Multiarray_Vector_d (struct Multiarray_Vector_d* a)
{
	assert(a != NULL);
	assert(a->data != NULL);

	if (a->owns_data) {
		const ptrdiff_t size = compute_size(a->order,a->extents);
		for (ptrdiff_t i = 0; i < size; ++i) {
			// It is possible that not all vectors present in the Multiarray were set.
			if (a->data[i])
				destructor_Vector_d(a->data[i]);
		}
	}
	free(a->data);
	free(a->extents);
	free(a);
}

void destructor_const_Multiarray_Vector_d (const struct const_Multiarray_Vector_d* a)
{
	destructor_Multiarray_Vector_d((struct Multiarray_Vector_d*)a);
}

void destructor_Multiarray_Matrix_d (struct Multiarray_Matrix_d* a)
{
	assert(a != NULL);
	assert(a->data != NULL);

	if (a->owns_data) {
		const ptrdiff_t size = compute_size(a->order,a->extents);
		for (ptrdiff_t i = 0; i < size; ++i) {
			// It is common that not all operators present in the Multiarray were set.
			if (a->data[i])
				destructor_Matrix_d(a->data[i]);
		}
	}
	free(a->data);
	free(a->extents);
	free(a);
}

void destructor_const_Multiarray_Matrix_d (const struct const_Multiarray_Matrix_d* a)
{
	destructor_Multiarray_Matrix_d((struct Multiarray_Matrix_d*)a);
}

void destructor_const_Multiarray2_Matrix_d (const struct const_Multiarray_Matrix_d* a[2])
{
	destructor_const_Multiarray_Matrix_d(a[0]);
	destructor_const_Multiarray_Matrix_d(a[1]);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
