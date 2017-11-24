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

#include "matrix_constructors.h"

#include <stdlib.h>
#include <assert.h>
#include "mkl.h"

#include "macros.h"
#include "definitions_mkl.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

struct Matrix_d* constructor_default_Matrix_d ()
{
	return constructor_move_Matrix_d_d('R',0,0,true,NULL);
}

const struct const_Matrix_d* constructor_default_const_Matrix_d ()
{
	return (const struct const_Matrix_d*) constructor_default_Matrix_d();
}

// Empty constructors *********************************************************************************************** //

struct Matrix_d* constructor_empty_Matrix_d (const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	double* data = malloc(ext_0*ext_1 * sizeof *data); // keep
	struct Matrix_d* a = constructor_move_Matrix_d_d(layout,ext_0,ext_1,true,data); // returned

	return a;
}

struct Matrix_i* constructor_empty_Matrix_i (const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	int* data = malloc(ext_0*ext_1 * sizeof *data); // keep
	struct Matrix_i* a = constructor_move_Matrix_i_i(layout,ext_0,ext_1,true,data); // returned

	return a;
}

// Copy constructors ************************************************************************************************ //

struct Matrix_i* constructor_copy_Matrix_i_i
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const int*const data_src)
{
	const ptrdiff_t size = ext_0*ext_1;

	int* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_i_i(layout,ext_0,ext_1,true,data);
}

struct Matrix_d* constructor_copy_Matrix_d_d
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const double*const data_src)
{
	const ptrdiff_t size = ext_0*ext_1;

	double* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_d_d(layout,ext_0,ext_1,true,data);
}

const struct const_Matrix_d* constructor_copy_const_Matrix_d_d
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const double*const data_src)
{
	return (const struct const_Matrix_d*) constructor_copy_Matrix_d_d(layout,ext_0,ext_1,data_src);
}

struct Matrix_d* constructor_copy_Matrix_d (const struct Matrix_d* src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	const double*const data_src = src->data;

	double* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_d_d(src->layout,src->ext_0,src->ext_1,true,data);
}

const struct const_Matrix_d* constructor_copy_const_Matrix_d (const struct const_Matrix_d*const src)
{
	return (const struct const_Matrix_d*) constructor_copy_Matrix_d((struct Matrix_d*)src);
}

const struct const_Matrix_d* constructor_copy_extract_const_Matrix_d
	(const struct const_Matrix_d*const src, const struct const_Vector_i*const indices)
{
	const char layout = src->layout;

	const ptrdiff_t i_max = indices->ext_0,
	                j_max = ( layout == 'R' ? src->ext_1 : src->ext_0 ),
	                size  = i_max * j_max;

	double* data = malloc(size * sizeof *data); // keep
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

	struct Matrix_d* dest = constructor_move_Matrix_d_d(layout,ext_0,ext_1,true,data); // returned
	const struct const_Matrix_d*const dest_c = NULL;
	const_constructor_move_Matrix_d(&dest_c,dest);
	return dest_c;
}

void const_constructor_copy_Matrix_d (const struct const_Matrix_d*const* dest, const struct const_Matrix_d*const src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	double* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] = src->data[i];

	struct Matrix_d* dest_m = constructor_move_Matrix_d_d(src->layout,src->ext_0,src->ext_1,true,data);
	const_constructor_move_Matrix_d(dest,dest_m);
}

// Move constructors ************************************************************************************************ //

struct Matrix_d* constructor_move_Matrix_d_d
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, double*const data)
{
    struct Matrix_d* dest = calloc(1,sizeof *dest); // returned

    dest->layout    = layout;
    dest->ext_0     = ext_0;
    dest->ext_1     = ext_1;
    dest->owns_data = owns_data;
    dest->data      = data;

    return dest;
}

const struct const_Matrix_d* constructor_move_const_Matrix_d_d
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, const double*const data)
{
	return (const struct const_Matrix_d*) constructor_move_Matrix_d_d(layout,ext_0,ext_1,owns_data,(double*)data);
}

struct Matrix_i* constructor_move_Matrix_i_i
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, int*const data)
{
    struct Matrix_i* dest = calloc(1,sizeof *dest); // returned

    dest->layout    = layout;
    dest->ext_0     = ext_0;
    dest->ext_1     = ext_1;
    dest->owns_data = owns_data;
    dest->data      = data;

    return dest;
}

void const_constructor_move_Matrix_d (const struct const_Matrix_d*const* dest, struct Matrix_d* src)
{
	*(struct const_Matrix_d**) dest = (struct const_Matrix_d*) src;
}

void const_constructor_move_const_Matrix_d (const struct const_Matrix_d*const* dest, const struct const_Matrix_d* src)
{
	*(const struct const_Matrix_d**) dest = src;
}

void const_constructor_move_Matrix_i (const struct const_Matrix_i*const* dest, struct Matrix_i* src)
{
	*(struct const_Matrix_i**) dest = (struct const_Matrix_i*) src;
}

// Special constructors ********************************************************************************************* //

struct Matrix_d* constructor_sub_block_Matrix_d
	(const ptrdiff_t row0, const ptrdiff_t col0, const ptrdiff_t n_row, const ptrdiff_t n_col,
	 const struct Matrix_d* src)
{
	const char layout = src->layout;

	ptrdiff_t ext_0 = n_row,
	          ext_1 = n_col;
	double* data = malloc(n_row*n_col * sizeof *data); // keep
	if (layout == 'R') {
		double* data_ptr = data;
		for (int i = 0; i < ext_0; ++i) {
			const double* data_src = get_row_Matrix_d(row0+i,src)+col0;
			for (int j = 0; j < ext_1; ++j)
				*data_ptr++ = *data_src++;
		}
	} else {
		EXIT_ADD_SUPPORT;
	}
	return constructor_move_Matrix_d_d(layout,ext_0,ext_1,true,data);
}


const struct const_Matrix_d* constructor_subset_const_Matrix_d
	(const struct const_Matrix_d* src, const struct const_Vector_i* ind_subset)
{
	const char layout = src->layout;
	ptrdiff_t ext_0 = -1,
	          ext_1 = -1;
	double* data = NULL;
	if (layout == 'R') {
		ext_0 = ind_subset->ext_0,
		ext_1 = src->ext_1;

		data = malloc(ext_0*ext_1 * sizeof *data); // keep
		double* data_ptr = data;
		for (int i = 0; i < ext_0; ++i) {
			const double* data_src = get_row_const_Matrix_d(ind_subset->data[i],src);
			for (int j = 0; j < ext_1; ++j)
				*data_ptr++ = *data_src++;
		}
	} else {
		EXIT_ADD_SUPPORT;
	}

	return constructor_move_const_Matrix_d_d(layout,ext_0,ext_1,true,data);
}

struct Matrix_d* constructor_copy_transpose_Matrix_d (struct Matrix_d* a, const bool mem_only)
{
	struct Matrix_d* a_t = constructor_copy_Matrix_d(a); // returned
	transpose_Matrix_d(a_t,mem_only);

	return a_t;
}

const struct const_Matrix_d* constructor_block_diagonal_const_Matrix_d
	(const struct const_Matrix_d* src_b, const ptrdiff_t n_blocks)
{
	const ptrdiff_t b_ext_0 = src_b->ext_0,
	                b_ext_1 = src_b->ext_1;

	const char layout = src_b->layout;
	const ptrdiff_t ext_0 = n_blocks*b_ext_0,
	                ext_1 = n_blocks*b_ext_1;

	struct Matrix_d* dest = constructor_empty_Matrix_d(layout,ext_0,ext_1); // returned
	set_to_value_Matrix_d(dest,0.0);

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

	return (const struct const_Matrix_d*) dest;
}

/*struct Matrix_d* constructor_diagonal_Matrix_from_Vector_d (const char layout, const struct Vector_d*const src)
{
	const ptrdiff_t ext_0 = src->ext_0,
	                size  = ext_0*ext_0;

	double* data = calloc(size , sizeof *data); // moved
	for (ptrdiff_t i = 0, ind = 0; i < size; i += ext_0+1, ++ind)
		data[i] = src->data[ind];

	return constructor_move_Matrix_d_d(layout,ext_0,ext_0,true,data);
}*/

struct Matrix_d* constructor_diagonal_Matrix_d_d (const char layout, const ptrdiff_t ext_0, const double val)
{
	const ptrdiff_t size = ext_0*ext_0;

	double* data = calloc(size , sizeof *data); // moved
	for (ptrdiff_t i = 0, ind = 0; i < size; i += ext_0+1, ++ind)
		data[i] = val;

	return constructor_move_Matrix_d_d(layout,ext_0,ext_0,true,data);
}

struct Matrix_d* constructor_identity_Matrix_d (const char layout, const ptrdiff_t ext_0)
{
	return constructor_diagonal_Matrix_d_d(layout,ext_0,1.0);
}

const struct const_Matrix_d* constructor_identity_const_Matrix_d (const char layout, const ptrdiff_t ext_0)
{
	return (const struct const_Matrix_d*) constructor_identity_Matrix_d(layout,ext_0);
}

struct Matrix_d* constructor_inverse_Matrix_d (struct Matrix_d* src)
{
	// The source matrix is copied as the entries would otherwise be modified while solving for for the inverse.
	struct Matrix_d* A = constructor_copy_Matrix_d(src);                                // destructed;
	struct Matrix_d* B = constructor_identity_Matrix_d(src->layout,src->ext_0);         // destructed;
	struct Matrix_d* X = constructor_empty_Matrix_d(src->layout,src->ext_0,src->ext_1); // returned;

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = A->ext_0,
	                 nrhs   = A->ext_0;
	double* a               = A->data,
	      * b               = B->data,
	      * x               = X->data;
	const lapack_int lda    = A->ext_0,
	                 ldb    = A->ext_0,
	                 ldx    = A->ext_0;
	lapack_int ipiv[n],
	           iter         = 0;

	const int info = LAPACKE_dsgesv(matrix_layout,n,nrhs,a,lda,ipiv,b,ldb,x,ldx,&iter);
	assert(info == 0);

	destructor_Matrix_d(A);
	destructor_Matrix_d(B);
	return X;
}

const struct const_Matrix_d* constructor_inverse_const_Matrix_d (const struct const_Matrix_d* src)
{
	return (const struct const_Matrix_d*) constructor_inverse_Matrix_d((struct Matrix_d*)src);
}

struct Matrix_d* constructor_sgesv_Matrix_d (struct Matrix_d* A_i, struct Matrix_d* B_i)
{
	assert(A_i->layout == B_i->layout); // Can be made flexible in future if necessary.
	assert(A_i->ext_0 == A_i->ext_1);

	// The source matrix is copied as the entries would otherwise be modified while solving the linear system.
	struct Matrix_d* A = constructor_copy_Matrix_d(A_i); // destructed;
	struct Matrix_d* X = constructor_empty_Matrix_d(A_i->layout,A_i->ext_0,B_i->ext_1); // returned;

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = A->ext_0,
	                 nrhs   = B_i->ext_1;
	double* a               = A->data,
	      * b               = B_i->data,
	      * x               = X->data;
	const lapack_int lda    = A->ext_0,
	                 ldb    = ( matrix_layout == LAPACK_COL_MAJOR ? n : nrhs ),
	                 ldx    = ldb;
	lapack_int ipiv[n],
	           iter         = 0;

	const int info = LAPACKE_dsgesv(matrix_layout,n,nrhs,a,lda,ipiv,b,ldb,x,ldx,&iter);
	assert(info == 0);

	destructor_Matrix_d(A);
	return X;
}

const struct const_Matrix_d* constructor_sgesv_const_Matrix_d
	(const struct const_Matrix_d* A_i, const struct const_Matrix_d* B_i)
{
	return (const struct const_Matrix_d*) constructor_sgesv_Matrix_d((struct Matrix_d*)A_i,(struct Matrix_d*)B_i);
}

struct Matrix_d* constructor_sysv_Matrix_d (struct Matrix_d* A_i, struct Matrix_d* B_i)
{
	assert(A_i->layout == B_i->layout); // Can be made flexible in future if necessary.
	assert(A_i->ext_0 == A_i->ext_1);

	// The source matrix is copied as the entries would otherwise be modified while solving the linear system.
/// \todo check if only half of the A matrix needs to be copied.
	struct Matrix_d* A = constructor_copy_Matrix_d(A_i); // destructed
	struct Matrix_d* X = constructor_copy_Matrix_d(B_i); // returned

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = A->ext_0,
	                 nrhs   = B_i->ext_1;
	double* a               = A->data,
	      * x               = X->data;
	const lapack_int lda    = A->ext_0,
	                 ldx    = ( matrix_layout == LAPACK_COL_MAJOR ? n : nrhs );
	lapack_int ipiv[n];

	const lapack_int info = LAPACKE_dsysv(matrix_layout,'U',n,nrhs,a,lda,ipiv,x,ldx);
	assert(info == 0);

	destructor_Matrix_d(A);
	return X;
}

const struct const_Matrix_d* constructor_sysv_const_Matrix_d
	(const struct const_Matrix_d* A_i, const struct const_Matrix_d* B_i)
{
	return (const struct const_Matrix_d*) constructor_sysv_Matrix_d((struct Matrix_d*)A_i,(struct Matrix_d*)B_i);
}

struct Matrix_d* constructor_mm_Matrix_d
	(const char trans_a_i, const char trans_b_i, const double alpha,
	 const struct const_Matrix_d*const a, const struct const_Matrix_d*const b, const char layout)
{
	const MKL_INT m = ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ),
	              n = ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 );

	struct Matrix_d* c = constructor_empty_Matrix_d(layout,m,n); // returned

	mm_d(trans_a_i,trans_b_i,alpha,0.0,a,b,c);

	return c;
}

const struct const_Matrix_d* constructor_mm_const_Matrix_d
	(const char trans_a_i, const char trans_b_i, const double alpha,
	 const struct const_Matrix_d*const a, const struct const_Matrix_d*const b, const char layout)
{
	return (const struct const_Matrix_d*) constructor_mm_Matrix_d(trans_a_i,trans_b_i,alpha,a,b,layout);
}

struct Matrix_d* constructor_mm_NN1R_Matrix_d
	(const struct const_Matrix_d*const a, const struct const_Matrix_d*const b)
{
	return constructor_mm_Matrix_d('N','N',1.0,a,b,'R');
}

const struct const_Matrix_d* constructor_mm_NN1R_const_Matrix_d
	(const struct const_Matrix_d*const a, const struct const_Matrix_d*const b)
{
	return (const struct const_Matrix_d*) constructor_mm_NN1R_Matrix_d(a,b);
}

struct Matrix_d* constructor_mm_NN1C_Matrix_d
	(const struct const_Matrix_d*const a, const struct const_Matrix_d*const b)
{
	return constructor_mm_Matrix_d('N','N',1.0,a,b,'C');
}

const struct const_Matrix_d* constructor_mm_NN1C_const_Matrix_d
	(const struct const_Matrix_d*const a, const struct const_Matrix_d*const b)
{
	return (const struct const_Matrix_d*) constructor_mm_NN1C_Matrix_d(a,b);
}

struct Matrix_d* constructor_mm_diag_Matrix_d
	(const double alpha, const struct const_Matrix_d*const a, const struct const_Vector_d*const b, const char side,
	 const bool invert_diag)
{
	struct Matrix_d* c = constructor_copy_Matrix_d((struct Matrix_d*)a);
	scale_Matrix_by_Vector_d(side,alpha,c,b,invert_diag);
	return c;
}

const struct const_Matrix_d* constructor_mm_diag_const_Matrix_d
	(const double alpha, const struct const_Matrix_d*const a, const struct const_Vector_d*const b, const char side,
	 const bool invert_diag)
{
	return (const struct const_Matrix_d*) constructor_mm_diag_Matrix_d(alpha,a,b,side,invert_diag);
}

void set_Matrix_from_Multiarray_d (struct Matrix_d* dest, struct Multiarray_d* src, const ptrdiff_t*const sub_indices)
{
	dest->layout    = src->layout;
	dest->owns_data = false;
	dest->ext_0     = src->extents[0];
	dest->ext_1     = src->extents[1];
	dest->data      = &src->data[compute_index_sub_container(src->order,2,src->extents,sub_indices)];
}
void set_const_Matrix_from_Multiarray_d
	(const struct const_Matrix_d* dest, const struct const_Multiarray_d* src, const ptrdiff_t*const sub_indices)
{
	set_Matrix_from_Multiarray_d((struct Matrix_d*)dest,(struct Multiarray_d*)src,sub_indices);
}

void set_Matrix_from_Multiarray_Matrix_d
	(struct Matrix_d* dest, struct Multiarray_Matrix_d* src, const ptrdiff_t*const sub_indices)
{
	struct Matrix_d* src_M = src->data[compute_index_sub_container(src->order,0,src->extents,sub_indices)];

	dest->layout    = src_M->layout;
	dest->ext_0     = src_M->ext_0;
	dest->ext_1     = src_M->ext_1;
	dest->owns_data = false;
	dest->data      = src_M->data;
}

void set_const_Matrix_from_Multiarray_Matrix_d
	(const struct const_Matrix_d* dest, const struct const_Multiarray_Matrix_d* src,
	 const ptrdiff_t*const sub_indices)
{
	set_Matrix_from_Multiarray_Matrix_d((struct Matrix_d*)dest,(struct Multiarray_Matrix_d*)src,sub_indices);
}

// Destructors ****************************************************************************************************** //

void destructor_Matrix_d (struct Matrix_d* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Matrix_d (const struct const_Matrix_d* a)
{
	destructor_Matrix_d((struct Matrix_d*)a);
}

void destructor_Matrix_i (struct Matrix_i* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Matrix_i (const struct const_Matrix_i* a)
{
	destructor_Matrix_i((struct Matrix_i*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
