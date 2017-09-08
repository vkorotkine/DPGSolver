// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "vector_constructors.h"

#include <string.h>

#include "macros.h"

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

// Destructors ****************************************************************************************************** //

void destructor_Vector_d (struct Vector_d* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_Vector_i (struct Vector_i* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_Vector_i_2 (struct Vector_i** a, const ptrdiff_t n_src, const bool owns_data)
{
	if (owns_data) {
		for (ptrdiff_t n = 0; n < n_src; n++)
			destructor_Vector_i(a[n]);
	}
	free(a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
