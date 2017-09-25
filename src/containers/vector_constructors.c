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

#include "vector_constructors.h"

#include <assert.h>
#include <string.h>
#include "mkl.h"

#include "macros.h"
#include "definitions_mkl.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

struct Vector_d* constructor_default_Vector_d ()
{
	struct Vector_d* dest = calloc(1,sizeof *dest); // returned
	dest->ext_0     = 0;
	dest->owns_data = true;
	dest->data      = NULL;

	return dest;
}

const struct const_Vector_d* constructor_default_const_Vector_d ()
{
	return (const struct const_Vector_d*) constructor_default_Vector_d();
}

struct Vector_i* constructor_default_Vector_i ()
{
	struct Vector_i* dest = calloc(1,sizeof *dest); // returned
	dest->ext_0     = 0;
	dest->owns_data = true;
	dest->data      = NULL;

	return dest;
}

struct Vector_i** constructor_default_Vector_i_2 (const ptrdiff_t n_dest)
{
	struct Vector_i** dest = malloc(n_dest * sizeof *dest); // returned;

	for (ptrdiff_t n = 0; n < n_dest; n++)
		dest[n] = constructor_default_Vector_i();

	return dest;
}

// Empty constructors *********************************************************************************************** //

struct Vector_d* constructor_empty_Vector_d (const ptrdiff_t ext_0)
{
	double* data = malloc(ext_0 * sizeof *data); // keep

	return constructor_move_Vector_d_d(ext_0,true,data);
}

struct Vector_i* constructor_empty_Vector_i (const ptrdiff_t ext_0)
{
	int* data = malloc(ext_0 * sizeof *data); // keep

	return constructor_move_Vector_i_i(ext_0,true,data);
}

// Copy constructors ************************************************************************************************ //

struct Vector_i* constructor_copy_Vector_i (const struct Vector_i*const src)
{
	const ptrdiff_t ext_0 = src->ext_0;
	const int*const data_src = src->data;

	int* data = malloc(ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	return constructor_move_Vector_i_i(ext_0,true,data);
}

struct Vector_d* constructor_copy_Vector_d_d (const ptrdiff_t ext_0, const double*const data_src)
{
	double* data = malloc(ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	return constructor_move_Vector_d_d(ext_0,true,data);
}

struct Vector_i* constructor_copy_Vector_i_i (const ptrdiff_t ext_0, const int*const data_src)
{
	int* data = malloc(ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	return constructor_move_Vector_i_i(ext_0,true,data);
}

const struct const_Vector_i* constructor_copy_const_Vector_i_i (const ptrdiff_t ext_0, const int*const data_src)
{
	return (const struct const_Vector_i*) constructor_copy_Vector_i_i(ext_0,data_src);
}

// Move constructors ************************************************************************************************ //

struct Vector_i* constructor_move_Vector_i_i (const ptrdiff_t ext_0, const bool owns_data, int*const data)
{
	struct Vector_i* dest = calloc(1,sizeof *dest); // returned

	dest->ext_0     = ext_0;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct Vector_d* constructor_move_Vector_d_d (const ptrdiff_t ext_0, const bool owns_data, double*const data)
{
	struct Vector_d* dest = calloc(1,sizeof *dest); // returned

	dest->ext_0     = ext_0;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct const_Vector_i* constructor_move_const_Vector_i_i
	(const ptrdiff_t ext_0, const bool owns_data, const int*const data)
{
	return (struct const_Vector_i*) constructor_move_Vector_i_i(ext_0,owns_data,(int*)data);
}

struct Vector_d* constructor_move_Vector_d_Matrix_d (struct Matrix_d* src)
{
	src->owns_data = false;
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	return constructor_move_Vector_d_d(size,true,src->data);
}

void const_constructor_move_Vector_d (const struct const_Vector_d*const* dest, struct Vector_d* src)
{
	*(struct const_Vector_d**) dest = (struct const_Vector_d*) src;
}

void const_constructor_move_Vector_i (const struct const_Vector_i*const* dest, struct Vector_i* src)
{
	*(struct const_Vector_i**) dest = (struct const_Vector_i*) src;
}

// Set constructors ************************************************************************************************* //

struct Vector_d* constructor_set_Vector_d_Multiarray_d (struct Multiarray_d* src, const ptrdiff_t* sub_indices)
{
	struct Vector_d* dest = constructor_default_Vector_d(); // returned
	set_Vector_from_Multiarray_d(dest,src,sub_indices);

	return dest;
}

const struct const_Vector_d* constructor_set_const_Vector_d_Multiarray_d
	(const struct const_Multiarray_d* src, const ptrdiff_t* sub_indices)
{
	return (const struct const_Vector_d*)
		constructor_set_Vector_d_Multiarray_d((struct Multiarray_d*)src,sub_indices);
}

// Special constructors ********************************************************************************************* //

struct Vector_d* constructor_sum_Vector_d_const_Matrix_d (const char sum_dir, const struct const_Matrix_d*const src)
{
	if (sum_dir != src->layout) {
		transpose_Matrix_d((struct Matrix_d*)src,true);
		struct Vector_d* dest = constructor_sum_Vector_d_const_Matrix_d(sum_dir,src);
		transpose_Matrix_d((struct Matrix_d*)src,true);
		return dest;
	}

	assert((sum_dir == 'R' || sum_dir == 'C'));

	const ptrdiff_t ext_0 = ( sum_dir == 'R' ? src->ext_0 : src->ext_1 ),
	                n_val = ( sum_dir == 'R' ? src->ext_1 : src->ext_0 );

	struct Vector_d* dest = constructor_empty_Vector_d(ext_0); // returned
	set_to_value_Vector_d(dest,0.0);

	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		const double* data_m = get_slice_const_Matrix_d(i,src);
		for (ptrdiff_t j = 0; j < n_val; ++j)
			dest->data[i] += data_m[j];
	}
	return dest;
}

const struct const_Vector_d* constructor_sum_const_Vector_d_const_Matrix_d
	(const char sum_dir, const struct const_Matrix_d*const src)
{
	return (const struct const_Vector_d*) constructor_sum_Vector_d_const_Matrix_d(sum_dir,src);
}

struct Vector_d* constructor_mv_Vector_d
	(const char trans_a_i, const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Vector_d*const b)
{
	const MKL_INT m = ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 );

	struct Vector_d* c = constructor_empty_Vector_d(m); // returned

	mv_d(trans_a_i,alpha,beta,a,b,c);

	return c;
}

const struct const_Vector_d* constructor_mv_const_Vector_d
	(const char trans_a_i, const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Vector_d*const b)
{
	return (const struct const_Vector_d*) constructor_mv_Vector_d(trans_a_i,alpha,beta,a,b);
}

struct Vector_d* constructor_sgesv_Vector_d (struct Matrix_d* A_i, struct Vector_d* B_i)
{
	// The source matrix is copied as the entries would otherwise be modified while solving the linear system.
	struct Matrix_d* A = constructor_copy_Matrix_d(A_i);         // destructed;
	struct Vector_d* X = constructor_empty_Vector_d(B_i->ext_0); // returned;

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = A->ext_0,
	                 nrhs   = 1;
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

const struct const_Vector_d* constructor_sgesv_const_Vector_d
	(const struct const_Matrix_d* A_i, const struct const_Vector_d* B_i)
{
	return (const struct const_Vector_d*) constructor_sgesv_Vector_d((struct Matrix_d*)A_i,(struct Vector_d*)B_i);
}

void set_Vector_from_Matrix_d (struct Vector_d* dest, struct Matrix_d* src, const ptrdiff_t*const sub_indices)
{
	dest->owns_data = false;
	if (src->layout == 'R') {
		dest->ext_0 = src->ext_1;
		dest->data  = get_row_Matrix_d(*sub_indices,src);
	} else {
		dest->ext_0 = src->ext_0;
		dest->data  = get_col_Matrix_d(*sub_indices,src);
	}
}

void set_const_Vector_from_Matrix_d
	(const struct const_Vector_d* dest, const struct const_Matrix_d* src, const ptrdiff_t*const sub_indices)
{
	set_Vector_from_Matrix_d((struct Vector_d*)dest,(struct Matrix_d*)src,sub_indices);
}

void set_Vector_from_Multiarray_d (struct Vector_d* dest, struct Multiarray_d* src, const ptrdiff_t*const sub_indices)
{
	dest->owns_data = false;
	dest->ext_0 = src->extents[0];
	dest->data  = &src->data[compute_index_sub_vector(src->order,src->extents,sub_indices)];
}

void set_const_Vector_from_Multiarray_d
	(const struct const_Vector_d* dest, const struct const_Multiarray_d* src, const ptrdiff_t*const sub_indices)
{
	set_Vector_from_Multiarray_d((struct Vector_d*)dest,(struct Multiarray_d*)src,sub_indices);
}


// Destructors ****************************************************************************************************** //

void destructor_Vector_d (struct Vector_d* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Vector_d (const struct const_Vector_d* a)
{
	destructor_Vector_d((struct Vector_d*)a);
}

void destructor_Vector_i (struct Vector_i* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Vector_i (const struct const_Vector_i* a)
{
	destructor_Vector_i((struct Vector_i*)a);
}

void destructor_Vector_i_2 (struct Vector_i** a, const ptrdiff_t n_src, const bool owns_data)
{
	assert(a != NULL);

	if (owns_data) {
		for (ptrdiff_t n = 0; n < n_src; n++)
			destructor_Vector_i(a[n]);
	}
	free(a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
