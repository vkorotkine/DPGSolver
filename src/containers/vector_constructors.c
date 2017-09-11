// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "vector_constructors.h"

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
	struct Vector_i* dest_m = constructor_move_Vector_i_i(ext_0,owns_data,(int*)data); // free

	const struct const_Vector_i*const dest = NULL;
	const_constructor_move_Vector_i(&dest,dest_m);

	return (struct const_Vector_i*)dest;
}

void const_constructor_move_Vector_d (const struct const_Vector_d*const* dest, struct Vector_d* src)
{
	*(struct const_Vector_d**) dest = (struct const_Vector_d*) src;
}

void const_constructor_move_Vector_i (const struct const_Vector_i*const* dest, struct Vector_i* src)
{
	*(struct const_Vector_i**) dest = (struct const_Vector_i*) src;
}

// Special constructors ********************************************************************************************* //

struct Vector_d* constructor_sum_Vector_d_const_Matrix_d (const char sum_dir, const struct const_Matrix_d*const src)
{
	if (!(sum_dir == 'R' || sum_dir == 'C'))
		EXIT_UNSUPPORTED;

	const ptrdiff_t ext_0 = ( sum_dir == 'R' ? src->ext_1 : src->ext_0 );

	struct Vector_d* dest = constructor_empty_Vector_d(ext_0); // returned
	for (ptrdiff_t j = 0; j < ext_0; ++j)
		dest->data[j] = 0.0;

	if (sum_dir != src->layout) {
		EXIT_ADD_SUPPORT;
	} else {
		if (src->layout == 'R') {
			const ptrdiff_t i_max = src->ext_0;
			for (ptrdiff_t i = 0; i < i_max; ++i) {
				const double* data_m = get_row_const_Matrix_d(i,src);
				for (ptrdiff_t j = 0; j < ext_0; ++j)
					dest->data[j] += data_m[j];
			}
		} else {
			EXIT_ADD_SUPPORT;
		}
	}

	return dest;
}

struct Vector_d* constructor_mv_Vector_d
	(const char layout, const char trans_a_i, const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Vector_d*const b)
{
	const MKL_INT m = ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 );

	struct Vector_d* c = constructor_empty_Vector_d(m); // returned

	mv_d(layout,trans_a_i,alpha,beta,a,b,c);

	return c;
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
	if (a == NULL)
		EXIT_DESTRUCTOR;

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
	if (a == NULL)
		EXIT_DESTRUCTOR;

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
	if (a == NULL)
		EXIT_DESTRUCTOR;

	if (owns_data) {
		for (ptrdiff_t n = 0; n < n_src; n++)
			destructor_Vector_i(a[n]);
	}
	free(a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
