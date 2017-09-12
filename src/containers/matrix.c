// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "matrix.h"

#include "macros.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

double* get_row_Matrix_d (const ptrdiff_t row, const struct Matrix_d* a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;
	return &a->data[row*(a->ext_1)];
}

double* get_col_Matrix_d (const ptrdiff_t col, const struct Matrix_d* a)
{
	if (a->layout != 'C')
		EXIT_UNSUPPORTED;
	return &a->data[col*(a->ext_0)];
}

const double* get_row_const_Matrix_d (const ptrdiff_t row, const struct const_Matrix_d*const a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;
	return &a->data[row*(a->ext_1)];
}

const double* get_col_const_Matrix_d (const ptrdiff_t col, const struct const_Matrix_d*const a)
{
	if (a->layout != 'C')
		EXIT_UNSUPPORTED;
	return &a->data[col*(a->ext_0)];
}

int* get_row_Matrix_i (const ptrdiff_t row, const struct Matrix_i* a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;

	return &a->data[row*(a->ext_1)];
}

int get_val_Matrix_i (const ptrdiff_t row, const ptrdiff_t col, const struct Matrix_i*const a)
{
	int*const data = a->data;
	switch (a->layout) {
	case 'R': {
		const ptrdiff_t ext_1 = a->ext_1;
		return data[row*ext_1+col];
	} case 'C':
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

int get_val_const_Matrix_i (const ptrdiff_t row, const ptrdiff_t col, const struct const_Matrix_i*const a)
{
	const int*const data = a->data;
	switch (a->layout) {
	case 'R': {
		const ptrdiff_t ext_1 = a->ext_1;
		return data[row*ext_1+col];
	} case 'C':
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void set_row_Matrix_d (const ptrdiff_t row, const struct Matrix_d* dest, const double*const data_src)
{
	if (dest->layout != 'R')
		EXIT_ADD_SUPPORT;

	double*const data = get_row_Matrix_d(row,dest);

	const ptrdiff_t i_max = dest->ext_1;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		data[i] = data_src[i];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
