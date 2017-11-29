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

#include "matrix.h"

#include <assert.h>

#include "macros.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to value setting functions.
 *
 *  \param dest The destination.
 *  \param src  The source.
 */
typedef void (*set_value_fptr)
	(double*const dest,
	 const double src
	);

/// \brief Version of \ref set_value_fptr inserting values.
void set_value_insert
	(double*const dest, ///< See brief.
	 const double src   ///< See brief.
	);

/// \brief Version of \ref set_value_fptr adding values.
void set_value_add
	(double*const dest, ///< See brief.
	 const double src   ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void swap_layout (char*const layout)
{
	*layout = ( *layout == 'R' ? 'C' : 'R' );
}

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

int* get_col_Matrix_i (const ptrdiff_t col, const struct Matrix_i* a)
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

const int* get_row_const_Matrix_i (const ptrdiff_t row, const struct const_Matrix_i* a)
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

void set_row_Matrix_i (const ptrdiff_t row, struct Matrix_i* dest, const int*const data_src)
{
	assert(dest->layout == 'R');

	int*const data = get_row_Matrix_i(row,dest);

	const ptrdiff_t i_max = dest->ext_1;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		data[i] = data_src[i];
}

void set_col_Matrix_i (const ptrdiff_t col, struct Matrix_i* dest, const int*const data_src)
{
	const ptrdiff_t ext_0 = dest->ext_0;
	if (dest->layout == col) {
		int*const data = get_col_Matrix_i(col,dest);

		for (ptrdiff_t i = 0; i < ext_0; ++i)
			data[i] = data_src[i];
	} else {
		const ptrdiff_t ext_1 = dest->ext_1;
		for (ptrdiff_t i = 0; i < ext_0; ++i)
			dest->data[i*ext_1+col] = data_src[i];
	}
}

void set_col_to_val_Matrix_i (const ptrdiff_t col, struct Matrix_i* dest, const int data_src)
{
	const ptrdiff_t ext_0 = dest->ext_0;
	if (dest->layout == col) {
		int*const data = get_col_Matrix_i(col,dest);

		for (ptrdiff_t i = 0; i < ext_0; ++i)
			data[i] = data_src;
	} else {
		const ptrdiff_t ext_1 = dest->ext_1;
		for (ptrdiff_t i = 0; i < ext_0; ++i)
			dest->data[i*ext_1+col] = data_src;
	}
}

void set_to_value_Matrix_i (struct Matrix_i*const a, const int val)
{
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] = val;
}

void set_to_value_Matrix_d (struct Matrix_d*const a, const double val)
{
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] = val;
}

void set_block_Matrix_d
	(struct Matrix_d* a, const struct const_Matrix_d* a_sub, const ptrdiff_t row0, const ptrdiff_t col0,
	 const char set_type)
{
	assert(a->layout == a_sub->layout); // Add support if required.
	assert(a->layout == 'R');

	set_value_fptr set_value = NULL;
	switch (set_type) {
		case 'i': set_value = set_value_insert; break;
		case 'a': set_value = set_value_add;    break;
		default:  EXIT_ERROR("Unsupported: %c.\n",set_type); break;
	}

	const ptrdiff_t ext_0 = a_sub->ext_0,
	                ext_1 = a_sub->ext_1;

	assert(row0+ext_0 <= a->ext_0);
	assert(col0+ext_1 <= a->ext_1);

	for (int i = 0, row = (int)row0; i < ext_0; ++i, ++row) {
		const double*const data_as = get_row_const_Matrix_d(i,a_sub);

		double* data_a = get_row_Matrix_d(row,a);
		data_a += col0;
		for (int j = 0; j < ext_1; ++j)
			set_value(&data_a[j],data_as[j]);
	}
}

char compute_opposite_layout (const char layout_i)
{
	assert((layout_i == 'R') || (layout_i == 'C'));
	return ( layout_i == 'R' ? 'C' : 'R' );
}

ptrdiff_t compute_index_Matrix
	(const ptrdiff_t i, const ptrdiff_t j, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const char layout)
{
	assert((layout == 'R') || (layout == 'C'));
	return ( layout == 'R' ?  i*ext_1+j : j*ext_0+i );
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

void set_value_insert (double*const dest, const double src)
{
	*dest = src;
}

void set_value_add (double*const dest, const double src)
{
	*dest += src;
}
