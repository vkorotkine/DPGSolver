// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "matrix_constructors.h"

#include <stdlib.h>

#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

struct Matrix_d* constructor_default_Matrix_d ()
{
	return constructor_move_Matrix_d_d('R',0,0,true,NULL);
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

struct Matrix_d* constructor_copy_Matrix_d (struct Matrix_d* src)
{
	return constructor_move_Matrix_d_d(src->layout,src->ext_0,src->ext_1,true,src->data);
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

void const_constructor_move_Matrix_i (const struct const_Matrix_i*const* dest, struct Matrix_i* src)
{
	*(struct const_Matrix_i**) dest = (struct const_Matrix_i*) src;
}

// Destructors ****************************************************************************************************** //

void destructor_Matrix_d (struct Matrix_d* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_Matrix_i (struct Matrix_i* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
