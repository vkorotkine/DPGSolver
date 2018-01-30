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
/// \file

#include <assert.h>
#include <string.h>
#include "definitions_mkl.h"
#include "mkl.h"

#include "macros.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

struct Vector_T* constructor_default_Vector_T ()
{
	struct Vector_T* dest = calloc(1,sizeof *dest); // returned
	dest->ext_0     = 0;
	dest->owns_data = true;
	dest->data      = NULL;

	return dest;
}

const struct const_Vector_T* constructor_default_const_Vector_T ()
{
	return (const struct const_Vector_T*) constructor_default_Vector_T();
}

struct Vector_T** constructor_default_Vector_T_2 (const ptrdiff_t n_dest)
{
	struct Vector_T** dest = malloc((size_t)n_dest * sizeof *dest); // returned;

	for (ptrdiff_t n = 0; n < n_dest; n++)
		dest[n] = constructor_default_Vector_T();

	return dest;
}

// Empty constructors *********************************************************************************************** //

struct Vector_T* constructor_empty_Vector_T (const ptrdiff_t ext_0)
{
	Type* data = malloc((size_t)ext_0 * sizeof *data); // keep

	return constructor_move_Vector_T_T(ext_0,true,data);
}

// Zero constructors ************************************************************************************************ //

struct Vector_T* constructor_zero_Vector_T (const ptrdiff_t ext_0)
{
	Type* data = calloc((size_t)ext_0 , sizeof *data); // keep

	return constructor_move_Vector_T_T(ext_0,true,data);
}

// Copy constructors ************************************************************************************************ //

struct Vector_T* constructor_copy_Vector_T (const struct Vector_T*const src)
{
	const ptrdiff_t ext_0 = src->ext_0;
	const Type*const data_src = src->data;

	Type* data = malloc((size_t)ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	return constructor_move_Vector_T_T(ext_0,true,data);
}

const struct const_Vector_T* constructor_copy_const_Vector_T (const struct const_Vector_T*const src)
{
	return (struct const_Vector_T*) constructor_copy_Vector_T((struct Vector_T*)src);
}

struct Vector_T* constructor_copy_Vector_T_T (const ptrdiff_t ext_0, const Type*const data_src)
{
	Type* data = malloc((size_t)ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	return constructor_move_Vector_T_T(ext_0,true,data);
}

const struct const_Vector_T* constructor_copy_const_Vector_T_T (const ptrdiff_t ext_0, const Type*const data_src)
{
	return (struct const_Vector_T*) constructor_copy_Vector_T_T(ext_0,data_src);
}

struct Vector_T* constructor_copy_Vector_T_Vector_R (struct Vector_R* src)
{
	const ptrdiff_t size = src->ext_0;
	const Real*const data_src = src->data;

	Type* data = calloc((size_t)size , sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = (Type)data_src[i];

	return constructor_move_Vector_T_T(src->ext_0,true,data);
}

const struct const_Vector_T* constructor_copy_const_Vector_T_Vector_R (const struct const_Vector_R* src)
{
	return (struct const_Vector_T*) constructor_copy_Vector_T_Vector_R((struct Vector_R*)src);
}

// Move constructors ************************************************************************************************ //

struct Vector_T* constructor_move_Vector_T_T (const ptrdiff_t ext_0, const bool owns_data, Type*const data)
{
	struct Vector_T* dest = calloc(1,sizeof *dest); // returned

	dest->ext_0     = ext_0;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct const_Vector_T* constructor_move_const_Vector_T_T
	(const ptrdiff_t ext_0, const bool owns_data, const Type*const data)
{
	return (struct const_Vector_T*) constructor_move_Vector_T_T(ext_0,owns_data,(Type*)data);
}

const struct const_Vector_T* constructor_move_const_Vector_Matrix_row_T
	(const int row, const struct const_Matrix_T* src, const int owns_data)
{
	return constructor_move_const_Vector_T_T(src->ext_1,owns_data,get_row_const_Matrix_T(row,src));
}

struct Vector_T* constructor_move_Vector_T_Matrix_T (struct Matrix_T* src)
{
	src->owns_data = false;
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	return constructor_move_Vector_T_T(size,true,src->data);
}

void const_constructor_move_Vector_T (const struct const_Vector_T*const* dest, struct Vector_T* src)
{
	*(struct const_Vector_T**) dest = (struct const_Vector_T*) src;
}

void const_constructor_move_const_Vector_T
	(const struct const_Vector_T*const* dest, const struct const_Vector_T* src)
{
	const_constructor_move_Vector_T(dest,(struct Vector_T*)src);
}

// Set constructors ************************************************************************************************* //
#ifdef TYPE_RC
struct Vector_T* constructor_set_Vector_T_Multiarray_T (struct Multiarray_T* src, const ptrdiff_t* sub_indices)
{
	struct Vector_T* dest = constructor_default_Vector_T(); // returned
	set_Vector_from_Multiarray_T(dest,src,sub_indices);

	return dest;
}

const struct const_Vector_T* constructor_set_const_Vector_T_Multiarray_T
	(const struct const_Multiarray_T* src, const ptrdiff_t* sub_indices)
{
	return (const struct const_Vector_T*)
		constructor_set_Vector_T_Multiarray_T((struct Multiarray_T*)src,sub_indices);
}

// Special constructors ********************************************************************************************* //
struct Vector_T* constructor_inverse_Vector_T (const struct const_Vector_T* src)
{
	const ptrdiff_t ext_0 = src->ext_0;
	const Type*const data_src = src->data;

	Type* data = malloc((size_t)ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++) {
		assert(data_src[i] != 0.0);
		data[i] = 1.0/data_src[i];
	}

	return constructor_move_Vector_T_T(ext_0,true,data);
}

const struct const_Vector_T* constructor_inverse_const_Vector_T (const struct const_Vector_T* src)
{
	return (struct const_Vector_T*) constructor_inverse_Vector_T(src);
}

const struct const_Vector_T* constructor_dot_mult_const_Vector_T
	(const Type alpha, const struct const_Vector_T* a, const struct const_Vector_T* b, const int n_repeated)
{
	assert(a->ext_0 == b->ext_0);

	const ptrdiff_t ext_0     = a->ext_0,
	                ext_0_rep = n_repeated*ext_0;
	Type* data_c = malloc((size_t)ext_0_rep * sizeof *data_c); // moved

	for (ptrdiff_t i = 0; i < ext_0; ++i)
		data_c[i] = alpha*a->data[i]*b->data[i];

	for (int n = 1; n < n_repeated; ++n) {
		const ptrdiff_t ind_c = ext_0*n;
		for (ptrdiff_t i = 0; i < ext_0; ++i)
			data_c[ind_c+i] = data_c[i];
	}

	return (const struct const_Vector_T*) constructor_move_Vector_T_T(ext_0_rep,true,data_c);
}

struct Vector_T* constructor_sum_Vectors_Vector_T
	(const Type alpha_0, struct Vector_T*const src_0, const Type alpha_1, struct Vector_T*const src_1)
{
	const ptrdiff_t ext_0 = src_0->ext_0;
	assert(ext_0 == src_1->ext_0);

	struct Vector_T*const dest = constructor_empty_Vector_T(ext_0); // returned
	for (int i = 0; i < ext_0; ++i)
		dest->data[i] = alpha_0*src_0->data[i] + alpha_1*src_1->data[i];

	return dest;
}

const struct const_Vector_T* constructor_sum_Vectors_const_Vector_T
	(const Type alpha_0, const struct const_Vector_T*const src_0, const Type alpha_1,
	 const struct const_Vector_T*const src_1)
{
	return (struct const_Vector_T*)
		constructor_sum_Vectors_Vector_T(alpha_0,(struct Vector_T*)src_0,alpha_1,(struct Vector_T*)src_1);
}

struct Vector_T* constructor_sum_Vector_T_const_Matrix_T (const char sum_dir, const struct const_Matrix_T*const src)
{
	assert((sum_dir == 'R' || sum_dir == 'C'));
	if (sum_dir != src->layout) {
		transpose_Matrix_T((struct Matrix_T*)src,true);
		struct Vector_T* dest = constructor_sum_Vector_T_const_Matrix_T(sum_dir,src);
		transpose_Matrix_T((struct Matrix_T*)src,true);
		return dest;
	}

	const ptrdiff_t ext_0 = ( sum_dir == 'R' ? src->ext_0 : src->ext_1 );
	struct Vector_T* dest = constructor_empty_Vector_T(ext_0); // returned

	set_to_sum_Vector_T(sum_dir,src,dest);

	return dest;
}

const struct const_Vector_T* constructor_sum_const_Vector_T_const_Matrix_T
	(const char sum_dir, const struct const_Matrix_T*const src)
{
	return (const struct const_Vector_T*) constructor_sum_Vector_T_const_Matrix_T(sum_dir,src);
}

struct Vector_T* constructor_mv_Vector_T
	(const char trans_a_i, const Type alpha, const struct const_Matrix_T*const a,
	 const struct const_Vector_T*const b)
{
	const MKL_INT m = (MKL_INT) ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 );

	struct Vector_T* c = constructor_empty_Vector_T(m); // returned

	mv_T(trans_a_i,alpha,0.0,a,b,c);

	return c;
}

const struct const_Vector_T* constructor_mv_const_Vector_T
	(const char trans_a_i, const Type alpha, const struct const_Matrix_T*const a,
	 const struct const_Vector_T*const b)
{
	return (const struct const_Vector_T*) constructor_mv_Vector_T(trans_a_i,alpha,a,b);
}

struct Vector_T* constructor_sgesv_Vector_T (struct Matrix_T* A_i, struct Vector_T* B_i)
{
	assert(A_i->ext_0 == A_i->ext_1);

	// The source matrix is copied as the entries would otherwise be modified while solving the linear system.
	struct Matrix_T* A = constructor_copy_Matrix_T(A_i);         // destructed;
	struct Vector_T* X = constructor_empty_Vector_T(B_i->ext_0); // returned;

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = (lapack_int) A->ext_0,
	                 nrhs   = 1;
	Type* a               = A->data,
	      * b               = B_i->data,
	      * x               = X->data;
	const lapack_int lda    = (lapack_int) A->ext_0,
	                 ldb    = ( matrix_layout == LAPACK_COL_MAJOR ? n : nrhs ),
	                 ldx    = ldb;
	lapack_int ipiv[n],
	           iter         = 0;

	const int info = LAPACKE_Tsgesv(matrix_layout,n,nrhs,a,lda,ipiv,b,ldb,x,ldx,&iter);
	assert(info == 0);

	destructor_Matrix_T(A);
	return X;
}

const struct const_Vector_T* constructor_sgesv_const_Vector_T
	(const struct const_Matrix_T* A_i, const struct const_Vector_T* B_i)
{
	return (const struct const_Vector_T*) constructor_sgesv_Vector_T((struct Matrix_T*)A_i,(struct Vector_T*)B_i);
}

void set_Vector_from_Matrix_T (struct Vector_T* dest, struct Matrix_T* src, const ptrdiff_t*const sub_indices)
{
	dest->owns_data = false;
	if (src->layout == 'R') {
		dest->ext_0 = src->ext_1;
		dest->data  = get_row_Matrix_T(*sub_indices,src);
	} else {
		dest->ext_0 = src->ext_0;
		dest->data  = get_col_Matrix_T(*sub_indices,src);
	}
}

void set_const_Vector_from_Matrix_T
	(const struct const_Vector_T* dest, const struct const_Matrix_T* src, const ptrdiff_t*const sub_indices)
{
	set_Vector_from_Matrix_T((struct Vector_T*)dest,(struct Matrix_T*)src,sub_indices);
}

void set_Vector_from_Multiarray_T (struct Vector_T* dest, struct Multiarray_T* src, const ptrdiff_t*const sub_indices)
{
	dest->owns_data = false;
	dest->ext_0 = src->extents[0];
	dest->data  = &src->data[compute_index_sub_vector(src->order,src->extents,sub_indices)];
}

void set_const_Vector_from_Multiarray_T
	(const struct const_Vector_T* dest, const struct const_Multiarray_T* src, const ptrdiff_t*const sub_indices)
{
	set_Vector_from_Multiarray_T((struct Vector_T*)dest,(struct Multiarray_T*)src,sub_indices);
}
#endif
// Destructors ****************************************************************************************************** //

void destructor_Vector_T (struct Vector_T* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Vector_T (const struct const_Vector_T* a)
{
	destructor_Vector_T((struct Vector_T*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
