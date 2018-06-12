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

#include <stdlib.h>
#include <assert.h>
#include "definitions_mkl.h"
#include "mkl.h"

#include "macros.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

struct Matrix_T* constructor_default_Matrix_T ()
{
	return constructor_move_Matrix_T_T('R',0,0,true,NULL);
}

const struct const_Matrix_T* constructor_default_const_Matrix_T ()
{
	return (const struct const_Matrix_T*) constructor_default_Matrix_T();
}

// Empty constructors *********************************************************************************************** //

struct Matrix_T* constructor_empty_Matrix_T (const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	Type* data = malloc((size_t)(ext_0*ext_1) * sizeof *data); // keep
	struct Matrix_T* a = constructor_move_Matrix_T_T(layout,ext_0,ext_1,true,data); // returned

	return a;
}

const struct const_Matrix_T* constructor_empty_const_Matrix_T
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	return (struct const_Matrix_T*) constructor_empty_Matrix_T(layout,ext_0,ext_1);
}

// Empty constructors *********************************************************************************************** //

struct Matrix_T* constructor_zero_Matrix_T (const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	Type* data = calloc((size_t)(ext_0*ext_1) , sizeof *data); // keep
	struct Matrix_T* a = constructor_move_Matrix_T_T(layout,ext_0,ext_1,true,data); // returned

	return a;
}

// Copy constructors ************************************************************************************************ //

struct Matrix_T* constructor_copy_Matrix_T (const struct Matrix_T* src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	const Type*const data_src = src->data;

	Type* data = malloc((size_t)size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_T_T(src->layout,src->ext_0,src->ext_1,true,data);
}

struct Matrix_T* constructor_copy_scale_Matrix_T (const struct Matrix_T*const src, const Type scale)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	const Type*const data_src = src->data;

	Type* data = malloc((size_t)size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = scale*data_src[i];

	return constructor_move_Matrix_T_T(src->layout,src->ext_0,src->ext_1,true,data);
}

const struct const_Matrix_T* constructor_copy_scale_const_Matrix_T
	(const struct const_Matrix_T*const src, const Type scale)
{
	return (struct const_Matrix_T*) constructor_copy_scale_Matrix_T((struct Matrix_T*)src,scale);
}

struct Matrix_T* constructor_copy_Matrix_T_T
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const Type*const data_src)
{
	const ptrdiff_t size = ext_0*ext_1;

	Type* data = malloc((size_t)size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_T_T(layout,ext_0,ext_1,true,data);
}

const struct const_Matrix_T* constructor_copy_const_Matrix_T_T
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const Type*const data_src)
{
	return (const struct const_Matrix_T*) constructor_copy_Matrix_T_T(layout,ext_0,ext_1,data_src);
}

struct Matrix_T* constructor_copy_Matrix_T_Matrix_R (struct Matrix_R* src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	const Real*const data_src = src->data;

	Type* data = calloc((size_t)size , sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = (Type)data_src[i];

	return constructor_move_Matrix_T_T(src->layout,src->ext_0,src->ext_1,true,data);
}

const struct const_Matrix_T* constructor_copy_const_Matrix_T_Matrix_R (const struct const_Matrix_R* src)
{
	return (struct const_Matrix_T*) constructor_copy_Matrix_T_Matrix_R((struct Matrix_R*)src);
}

const struct const_Matrix_T* constructor_copy_const_Matrix_T (const struct const_Matrix_T*const src)
{
	return (const struct const_Matrix_T*) constructor_copy_Matrix_T((struct Matrix_T*)src);
}

const struct const_Matrix_T* constructor_copy_extract_const_Matrix_T
	(const struct const_Matrix_T*const src, const struct const_Vector_i*const indices)
{
	const char layout = src->layout;

	const ptrdiff_t i_max = indices->ext_0,
	                j_max = ( layout == 'R' ? src->ext_1 : src->ext_0 ),
	                size  = i_max * j_max;

	Type* data = malloc((size_t)size * sizeof *data); // keep
	ptrdiff_t ind = 0;
	for (ptrdiff_t i = 0; i < i_max; ++i) {
		for (ptrdiff_t j = 0; j < j_max; ++j) {
			data[ind] = src->data[(indices->data[i])*j_max+j];
			++ind;
		}
	}

	ptrdiff_t ext_0 = 0,
	          ext_1 = 0;
	if (layout == 'R') {
		ext_0 = i_max;
		ext_1 = j_max;
	} else {
		ext_0 = j_max;
		ext_1 = i_max;
	}

	return (struct const_Matrix_T*) constructor_move_Matrix_T_T(layout,ext_0,ext_1,true,data);
}

void const_constructor_copy_Matrix_T (const struct const_Matrix_T*const* dest, const struct const_Matrix_T*const src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	Type* data = malloc((size_t)size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] = src->data[i];

	struct Matrix_T* dest_m = constructor_move_Matrix_T_T(src->layout,src->ext_0,src->ext_1,true,data);
	const_constructor_move_Matrix_T(dest,dest_m);
}

// Move constructors ************************************************************************************************ //

struct Matrix_T* constructor_move_Matrix_T_T
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, Type*const data)
{
    struct Matrix_T* dest = calloc(1,sizeof *dest); // returned

    dest->layout    = layout;
    dest->ext_0     = ext_0;
    dest->ext_1     = ext_1;
    dest->owns_data = owns_data;
    dest->data      = data;

    return dest;
}

const struct const_Matrix_T* constructor_move_const_Matrix_T_T
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, const Type*const data)
{
	return (const struct const_Matrix_T*) constructor_move_Matrix_T_T(layout,ext_0,ext_1,owns_data,(Type*)data);
}

void const_constructor_move_Matrix_T (const struct const_Matrix_T*const* dest, struct Matrix_T* src)
{
	*(struct const_Matrix_T**) dest = (struct const_Matrix_T*) src;
}

void const_constructor_move_const_Matrix_T (const struct const_Matrix_T*const* dest, const struct const_Matrix_T* src)
{
	*(const struct const_Matrix_T**) dest = src;
}

// Special constructors (only available for real/complex types) ***************************************************** //
#ifdef TYPE_RC

struct Matrix_T* constructor_copy_permute_Matrix_T
	(const struct Matrix_T*const src, const struct const_Vector_i* p_V, const char perm_layout)
{
	assert(perm_layout == 'R'); // Check that all is as expected if removed.
	assert(src->layout == 'R');

	struct Matrix_T*const dest = constructor_copy_Matrix_T(src); // returned
	permute_Matrix_T_V(dest,p_V);

	return dest;
}

const struct const_Matrix_T* constructor_copy_permute_const_Matrix_T
	(const struct const_Matrix_T*const src, const struct const_Vector_i*const p_V, const char perm_layout)
{
	return (struct const_Matrix_T*) constructor_copy_permute_Matrix_T((struct Matrix_T*)src,p_V,perm_layout);
}

struct Matrix_T* constructor_sub_block_Matrix_T
	(const ptrdiff_t row0, const ptrdiff_t col0, const ptrdiff_t n_row, const ptrdiff_t n_col,
	 const struct Matrix_T* src)
{
	const char layout = src->layout;

	ptrdiff_t ext_0 = n_row,
	          ext_1 = n_col;
	Type* data = malloc((size_t)(n_row*n_col) * sizeof *data); // keep
	if (layout == 'R') {
		Type* data_ptr = data;
		for (int i = 0; i < ext_0; ++i) {
			const Type* data_src = get_row_Matrix_T(row0+i,src)+col0;
			for (int j = 0; j < ext_1; ++j)
				*data_ptr++ = *data_src++;
		}
	} else {
		EXIT_ADD_SUPPORT;
	}
	return constructor_move_Matrix_T_T(layout,ext_0,ext_1,true,data);
}


const struct const_Matrix_T* constructor_subset_const_Matrix_T
	(const struct const_Matrix_T* src, const struct const_Vector_i* ind_subset)
{
	const char layout = src->layout;
	ptrdiff_t ext_0 = -1,
	          ext_1 = -1;
	Type* data = NULL;
	if (layout == 'R') {
		ext_0 = ind_subset->ext_0,
		ext_1 = src->ext_1;

		data = malloc((size_t)(ext_0*ext_1) * sizeof *data); // keep
		Type* data_ptr = data;
		for (int i = 0; i < ext_0; ++i) {
			const Type* data_src = get_row_const_Matrix_T(ind_subset->data[i],src);
			for (int j = 0; j < ext_1; ++j)
				*data_ptr++ = *data_src++;
		}
	} else {
		EXIT_ADD_SUPPORT;
	}

	return constructor_move_const_Matrix_T_T(layout,ext_0,ext_1,true,data);
}

struct Matrix_T* constructor_copy_transpose_Matrix_T (struct Matrix_T* a, const bool mem_only)
{
	struct Matrix_T* a_t = constructor_copy_Matrix_T(a); // returned
	transpose_Matrix_T(a_t,mem_only);

	return a_t;
}

const struct const_Matrix_T* constructor_block_diagonal_const_Matrix_T
	(const struct const_Matrix_T* src_b, const ptrdiff_t n_blocks)
{
	const ptrdiff_t b_ext_0 = src_b->ext_0,
	                b_ext_1 = src_b->ext_1;

	const char layout = src_b->layout;
	const ptrdiff_t ext_0 = n_blocks*b_ext_0,
	                ext_1 = n_blocks*b_ext_1;

	struct Matrix_T* dest = constructor_empty_Matrix_T(layout,ext_0,ext_1); // returned
	set_to_value_Matrix_T(dest,0.0);

	assert(layout == 'R'); // Can be made flexible if necessary.
	for (ptrdiff_t n = 0; n < n_blocks; ++n) {
		const ptrdiff_t ind_row = n*b_ext_0,
		                ind_col = n*b_ext_1;
		for (ptrdiff_t i = 0; i < b_ext_0; ++i) {
			const ptrdiff_t ind_dest = compute_index_Matrix(ind_row+i,ind_col+0,ext_0,ext_1,layout),
			                ind_src  = compute_index_Matrix(i,0,b_ext_0,b_ext_1,layout);
			for (ptrdiff_t j = 0; j < b_ext_1; ++j)
				dest->data[ind_dest+j] = src_b->data[ind_src+j];
		}
	}

	return (const struct const_Matrix_T*) dest;
}

struct Matrix_T* constructor_diagonal_Matrix_T_T (const char layout, const ptrdiff_t ext_0, const Type val)
{
	const ptrdiff_t size = ext_0*ext_0;

	Type* data = calloc((size_t)size , sizeof *data); // moved
	for (ptrdiff_t i = 0, ind = 0; i < size; i += ext_0+1, ++ind)
		data[i] = val;

	return constructor_move_Matrix_T_T(layout,ext_0,ext_0,true,data);
}

struct Matrix_T* constructor_identity_Matrix_T (const char layout, const ptrdiff_t ext_0)
{
	return constructor_diagonal_Matrix_T_T(layout,ext_0,1.0);
}

const struct const_Matrix_T* constructor_identity_const_Matrix_T (const char layout, const ptrdiff_t ext_0)
{
	return (const struct const_Matrix_T*) constructor_identity_Matrix_T(layout,ext_0);
}

struct Matrix_T* constructor_inverse_Matrix_T (struct Matrix_T* src)
{
	// The source matrix is copied as the entries would otherwise be modified while solving for the inverse.
	struct Matrix_T* A = constructor_copy_Matrix_T(src);                                // destructed;
	struct Matrix_T* B = constructor_identity_Matrix_T(src->layout,src->ext_0);         // destructed;
	struct Matrix_T* X = constructor_empty_Matrix_T(src->layout,src->ext_0,src->ext_1); // returned;

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = (lapack_int)A->ext_0,
	                 nrhs   = (lapack_int)A->ext_0;
	Type* a                 = A->data,
	    * b                 = B->data,
	    * x                 = X->data;
	const lapack_int lda    = (lapack_int)A->ext_0,
	                 ldb    = (lapack_int)A->ext_0,
	                 ldx    = (lapack_int)A->ext_0;
	lapack_int ipiv[n],
	           iter         = 0;

	const int info = LAPACKE_Tsgesv(matrix_layout,n,nrhs,a,lda,ipiv,b,ldb,x,ldx,&iter);
	assert(info == 0);

	destructor_Matrix_T(A);
	destructor_Matrix_T(B);
	return X;
}

const struct const_Matrix_T* constructor_inverse_const_Matrix_T (const struct const_Matrix_T* src)
{
	return (const struct const_Matrix_T*) constructor_inverse_Matrix_T((struct Matrix_T*)src);
}

struct Matrix_T* constructor_sgesv_Matrix_T (struct Matrix_T* A_i, struct Matrix_T* B_i)
{
	assert(A_i->layout == B_i->layout); // Can be made flexible in future if necessary.
	assert(A_i->ext_0 == A_i->ext_1);

	// The source matrix is copied as the entries would otherwise be modified while solving the linear system.
	struct Matrix_T* A = constructor_copy_Matrix_T(A_i); // destructed
	struct Matrix_T* X = constructor_empty_Matrix_T(A_i->layout,A_i->ext_0,B_i->ext_1); // returned

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = (lapack_int)A->ext_0,
	                 nrhs   = (lapack_int)B_i->ext_1;
	Type* a               = A->data,
	      * b               = B_i->data,
	      * x               = X->data;
	const lapack_int lda    = (lapack_int)A->ext_0,
	                 ldb    = ( matrix_layout == LAPACK_COL_MAJOR ? n : nrhs ),
	                 ldx    = ldb;
	lapack_int ipiv[n],
	           iter         = 0;

	const int info = LAPACKE_Tsgesv(matrix_layout,n,nrhs,a,lda,ipiv,b,ldb,x,ldx,&iter);
	assert(info == 0);

	destructor_Matrix_T(A);
	return X;
}

const struct const_Matrix_T* constructor_sgesv_const_Matrix_T
	(const struct const_Matrix_T* A_i, const struct const_Matrix_T* B_i)
{
	return (const struct const_Matrix_T*) constructor_sgesv_Matrix_T((struct Matrix_T*)A_i,(struct Matrix_T*)B_i);
}

struct Matrix_T* constructor_sysv_Matrix_T (struct Matrix_T* A_i, struct Matrix_T* B_i)
{
	assert(A_i->ext_0 == A_i->ext_1);

	// The source matrix is copied as the entries would otherwise be modified while solving the linear system.
	struct Matrix_T* A = constructor_copy_Matrix_T(A_i); // destructed
	struct Matrix_T* X = constructor_copy_Matrix_T(B_i); // returned

	if (A->layout != X->layout)
		swap_layout(&A->layout); // As A is symmetric, we simply swap the layout here.

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = (lapack_int)A->ext_0,
	                 nrhs   = (lapack_int)B_i->ext_1;
	Type* a                 = A->data,
	    * x                 = X->data;
	const lapack_int lda    = (lapack_int)A->ext_0,
	                 ldx    = ( matrix_layout == LAPACK_COL_MAJOR ? n : nrhs );
	lapack_int ipiv[n];

	const lapack_int info = LAPACKE_Tsysv(matrix_layout,'U',n,nrhs,a,lda,ipiv,x,ldx);
	assert(info == 0);

	destructor_Matrix_T(A);
	return X;
}

const struct const_Matrix_T* constructor_sysv_const_Matrix_T
	(const struct const_Matrix_T* A_i, const struct const_Matrix_T* B_i)
{
	return (const struct const_Matrix_T*) constructor_sysv_Matrix_T((struct Matrix_T*)A_i,(struct Matrix_T*)B_i);
}

struct Matrix_T* constructor_mm_Matrix_T
	(const char trans_a_i, const char trans_b_i, const Real alpha,
	 const struct const_Matrix_T*const a, const struct const_Matrix_T*const b, const char layout)
{
	const MKL_INT m = (MKL_INT) ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ),
	              n = (MKL_INT) ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 );

	struct Matrix_T* c = constructor_empty_Matrix_T(layout,m,n); // returned

	mm_T(trans_a_i,trans_b_i,alpha,0.0,a,b,c);

	return c;
}

const struct const_Matrix_T* constructor_mm_const_Matrix_T
	(const char trans_a_i, const char trans_b_i, const Real alpha,
	 const struct const_Matrix_T*const a, const struct const_Matrix_T*const b, const char layout)
{
	return (const struct const_Matrix_T*) constructor_mm_Matrix_T(trans_a_i,trans_b_i,alpha,a,b,layout);
}

struct Matrix_T* constructor_mm_RT_Matrix_T
	(const char trans_a_i, const char trans_b_i, const Real alpha,
	 const struct const_Matrix_R*const a, const struct const_Matrix_T*const b, const char layout)
{
	const MKL_INT m = (MKL_INT) ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ),
	              n = (MKL_INT) ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 );

	struct Matrix_T* c = constructor_empty_Matrix_T(layout,m,n); // returned

	mm_RTT(trans_a_i,trans_b_i,alpha,0.0,a,b,c);

	return c;
}

const struct const_Matrix_T* constructor_mm_RT_const_Matrix_T
	(const char trans_a_i, const char trans_b_i, const Real alpha,
	 const struct const_Matrix_R*const a, const struct const_Matrix_T*const b, const char layout)
{
	return (const struct const_Matrix_T*) constructor_mm_RT_Matrix_T(trans_a_i,trans_b_i,alpha,a,b,layout);
}

struct Matrix_T* constructor_mm_TR_Matrix_T
	(const char trans_a_i, const char trans_b_i, const Real alpha,
	 const struct const_Matrix_T*const a, const struct const_Matrix_R*const b, const char layout)
{
	const MKL_INT m = (MKL_INT) ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ),
	              n = (MKL_INT) ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 );

	struct Matrix_T* c = constructor_empty_Matrix_T(layout,m,n); // returned

	mm_TRT(trans_a_i,trans_b_i,alpha,0.0,a,b,c);

	return c;
}

const struct const_Matrix_T* constructor_mm_TR_const_Matrix_T
	(const char trans_a_i, const char trans_b_i, const Real alpha,
	 const struct const_Matrix_T*const a, const struct const_Matrix_R*const b, const char layout)
{
	return (const struct const_Matrix_T*) constructor_mm_TR_Matrix_T(trans_a_i,trans_b_i,alpha,a,b,layout);
}

struct Matrix_T* constructor_mm_NN1R_Matrix_T
	(const struct const_Matrix_T*const a, const struct const_Matrix_T*const b)
{
	return constructor_mm_Matrix_T('N','N',1.0,a,b,'R');
}

const struct const_Matrix_T* constructor_mm_NN1R_const_Matrix_T
	(const struct const_Matrix_T*const a, const struct const_Matrix_T*const b)
{
	return (const struct const_Matrix_T*) constructor_mm_NN1R_Matrix_T(a,b);
}

struct Matrix_T* constructor_mm_NN1C_Matrix_T
	(const struct const_Matrix_T*const a, const struct const_Matrix_T*const b)
{
	return constructor_mm_Matrix_T('N','N',1.0,a,b,'C');
}

const struct const_Matrix_T* constructor_mm_NN1C_const_Matrix_T
	(const struct const_Matrix_T*const a, const struct const_Matrix_T*const b)
{
	return (const struct const_Matrix_T*) constructor_mm_NN1C_Matrix_T(a,b);
}

struct Matrix_T* constructor_mm_diag_Matrix_T_R
	(const Real alpha, const struct const_Matrix_T*const a, const struct const_Vector_R*const b, const char side,
	 const bool invert_diag)
{
	struct Matrix_T* c = constructor_copy_Matrix_T((struct Matrix_T*)a);
	scale_Matrix_T_by_Vector_R(side,alpha,c,b,invert_diag);
	return c;
}

const struct const_Matrix_T* constructor_mm_diag_const_Matrix_T_R
	(const Real alpha, const struct const_Matrix_T*const a, const struct const_Vector_R*const b, const char side,
	 const bool invert_diag)
{
	return (const struct const_Matrix_T*) constructor_mm_diag_Matrix_T_R(alpha,a,b,side,invert_diag);
}

struct Matrix_T* constructor_mm_diag_Matrix_R_T
	(const Real alpha, const struct const_Matrix_R*const a, const struct const_Vector_T*const b, const char side,
	 const bool invert_diag)
{
	struct Matrix_T* c = constructor_copy_Matrix_T_Matrix_R((struct Matrix_R*)a); // returned
	scale_Matrix_by_Vector_T(side,alpha,c,b,invert_diag);
	return c;
}

const struct const_Matrix_T* constructor_mm_diag_const_Matrix_R_T
	(const Real alpha, const struct const_Matrix_R*const a, const struct const_Vector_T*const b, const char side,
	 const bool invert_diag)
{
	return (const struct const_Matrix_T*) constructor_mm_diag_Matrix_R_T(alpha,a,b,side,invert_diag);
}

struct Matrix_T* constructor_mm_diag_Matrix_T
	(const Real alpha, const struct const_Matrix_T*const a, const struct const_Vector_T*const b, const char side,
	 const bool invert_diag)
{
	struct Matrix_T* c = constructor_copy_Matrix_T((struct Matrix_T*)a);
	scale_Matrix_by_Vector_T(side,alpha,c,b,invert_diag);
	return c;
}

const struct const_Matrix_T* constructor_mm_diag_const_Matrix_T
	(const Real alpha, const struct const_Matrix_T*const a, const struct const_Vector_T*const b, const char side,
	 const bool invert_diag)
{
	return (const struct const_Matrix_T*) constructor_mm_diag_Matrix_T(alpha,a,b,side,invert_diag);
}
#endif

void set_Matrix_from_Multiarray_T (struct Matrix_T* dest, struct Multiarray_T* src, const ptrdiff_t*const sub_indices)
{
	dest->layout    = src->layout;
	dest->owns_data = false;
	dest->ext_0     = src->extents[0];
	dest->ext_1     = src->extents[1];
	dest->data      = &src->data[compute_index_sub_container(src->order,2,src->extents,sub_indices)];
}

void set_const_Matrix_from_Multiarray_T
	(const struct const_Matrix_T* dest, const struct const_Multiarray_T* src, const ptrdiff_t*const sub_indices)
{
	set_Matrix_from_Multiarray_T((struct Matrix_T*)dest,(struct Multiarray_T*)src,sub_indices);
}

void set_Matrix_from_Multiarray_Matrix_T
	(struct Matrix_T* dest, struct Multiarray_Matrix_T* src, const ptrdiff_t*const sub_indices)
{
	struct Matrix_T* src_M = src->data[compute_index_sub_container(src->order,0,src->extents,sub_indices)];

	dest->layout    = src_M->layout;
	dest->ext_0     = src_M->ext_0;
	dest->ext_1     = src_M->ext_1;
	dest->owns_data = false;
	dest->data      = src_M->data;
}

void set_const_Matrix_from_Multiarray_Matrix_T
	(const struct const_Matrix_T* dest, const struct const_Multiarray_Matrix_T* src,
	 const ptrdiff_t*const sub_indices)
{
	set_Matrix_from_Multiarray_Matrix_T((struct Matrix_T*)dest,(struct Multiarray_Matrix_T*)src,sub_indices);
}

// Destructors ****************************************************************************************************** //

void destructor_Matrix_T (struct Matrix_T* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Matrix_T (const struct const_Matrix_T* a)
{
	destructor_Matrix_T((struct Matrix_T*)a);
}

void destructor_conditional_Matrix_T (struct Matrix_T* a)
{
	if (a)
		destructor_Matrix_T(a);
}

void destructor_conditional_const_Matrix_T (const struct const_Matrix_T* a)
{
	destructor_conditional_Matrix_T ((struct Matrix_T*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
