// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "matrix.h"

#include <assert.h>

#include "macros.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

double* get_row_Matrix_d (const ptrdiff_t row, const struct Matrix_d* a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

double* get_col_Matrix_d (const ptrdiff_t col, const struct Matrix_d* a)
{
	assert(a->layout == 'C');
	return &a->data[col*(a->ext_0)];
}

double* get_slice_Matrix_d (const ptrdiff_t slice, const struct Matrix_d* a)
{
	return ( a->layout == 'R' ? get_row_Matrix_d(slice,a) : get_col_Matrix_d(slice,a) );
}

const double* get_row_const_Matrix_d (const ptrdiff_t row, const struct const_Matrix_d*const a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

const double* get_col_const_Matrix_d (const ptrdiff_t col, const struct const_Matrix_d*const a)
{
	assert(a->layout == 'C');
	return &a->data[col*(a->ext_0)];
}

const double* get_slice_const_Matrix_d (const ptrdiff_t slice, const struct const_Matrix_d* a)
{
	return (const double*) get_slice_Matrix_d(slice,(const struct Matrix_d*)a);
}

int* get_row_Matrix_i (const ptrdiff_t row, const struct Matrix_i* a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

int get_val_Matrix_i (const ptrdiff_t row, const ptrdiff_t col, const struct Matrix_i*const a)
{
	assert((a->layout == 'R') || (a->layout == 'C'));

	int*const data = a->data;

	return ( a->layout == 'R' ? data[row*(a->ext_1)+col] : data[col*(a->ext_0)+row]);
}

int get_val_const_Matrix_i (const ptrdiff_t row, const ptrdiff_t col, const struct const_Matrix_i*const a)
{
	return get_val_Matrix_i(row,col,(struct Matrix_i*)a);
}

void set_row_Matrix_d (const ptrdiff_t row, struct Matrix_d* dest, const double*const data_src)
{
	assert(dest->layout == 'R');

	double*const data = get_row_Matrix_d(row,dest);

	const ptrdiff_t i_max = dest->ext_1;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		data[i] = data_src[i];
}

void set_to_value_Matrix_d (struct Matrix_d*const a, const double val)
{
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] = val;
}

char compute_opposite_layout (const char layout_i)
{
	assert((layout_i == 'R') || (layout_i == 'C'));
	return ( layout_i == 'R' ? 'C' : 'R' );
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
