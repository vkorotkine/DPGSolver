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

#include <assert.h>
#include <stddef.h>

#include "macros.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to value setting functions.
 *
 *  \param dest The destination.
 *  \param src  The source.
 */
typedef void (*set_value_fptr_T)
	(Type*const dest,
	 const Type src
	);

/// \brief Version of \ref set_value_fptr_T inserting values. \{
void set_value_insert_T
	(Type*const dest,
	 const Type src
	); /// \}

/// \brief Version of \ref set_value_fptr_T adding values. \{
void set_value_add_T
	(Type*const dest,
	 const Type src
	); /// \}

// Interface functions ********************************************************************************************** //

Type* get_row_Matrix_T (const ptrdiff_t row, const struct Matrix_T* a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

const Type* get_row_const_Matrix_T (const ptrdiff_t row, const struct const_Matrix_T*const a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

Type* get_col_Matrix_T (const ptrdiff_t col, const struct Matrix_T* a)
{
	assert(a->layout == 'C');
	return &a->data[col*(a->ext_0)];
}

const Type* get_col_const_Matrix_T (const ptrdiff_t col, const struct const_Matrix_T*const a)
{
	assert(a->layout == 'C');
	return &a->data[col*(a->ext_0)];
}

Type* get_slice_Matrix_T (const ptrdiff_t slice, const struct Matrix_T* a)
{
	return ( a->layout == 'R' ? get_row_Matrix_T(slice,a) : get_col_Matrix_T(slice,a) );
}

const Type* get_slice_const_Matrix_T (const ptrdiff_t slice, const struct const_Matrix_T* a)
{
	return (const Type*) get_slice_Matrix_T(slice,(const struct Matrix_T*)a);
}

Type get_val_Matrix_T (const ptrdiff_t row, const ptrdiff_t col, const struct Matrix_T*const a)
{
	assert((a->layout == 'R') || (a->layout == 'C'));

	Type*const data = a->data;
	return ( a->layout == 'R' ? data[row*(a->ext_1)+col] : data[col*(a->ext_0)+row]);
}

Type get_val_const_Matrix_T (const ptrdiff_t row, const ptrdiff_t col, const struct const_Matrix_T*const a)
{
	return get_val_Matrix_T(row,col,(struct Matrix_T*)a);
}

void set_row_Matrix_T (const ptrdiff_t row, struct Matrix_T* dest, const Type*const data_src)
{
	assert(dest->layout == 'R');

	Type*const data = get_row_Matrix_T(row,dest);

	const ptrdiff_t i_max = dest->ext_1;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		data[i] = data_src[i];
}

void set_col_Matrix_T (const ptrdiff_t col, struct Matrix_T* dest, const Type*const data_src)
{
	const ptrdiff_t ext_0 = dest->ext_0;
	if (dest->layout == col) {
		Type*const data = get_col_Matrix_T(col,dest);
		for (ptrdiff_t i = 0; i < ext_0; ++i)
			data[i] = data_src[i];
	} else {
		const ptrdiff_t ext_1 = dest->ext_1;
		for (ptrdiff_t i = 0; i < ext_0; ++i)
			dest->data[i*ext_1+col] = data_src[i];
	}
}

void set_col_to_val_Matrix_T (const ptrdiff_t col, struct Matrix_T* dest, const Type data_src)
{
	const ptrdiff_t ext_0 = dest->ext_0;
	if (dest->layout == col) {
		Type*const data = get_col_Matrix_T(col,dest);
		for (ptrdiff_t i = 0; i < ext_0; ++i)
			data[i] = data_src;
	} else {
		const ptrdiff_t ext_1 = dest->ext_1;
		for (ptrdiff_t i = 0; i < ext_0; ++i)
			dest->data[i*ext_1+col] = data_src;
	}
}

void set_to_value_Matrix_T (struct Matrix_T*const a, const Type val)
{
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] = val;
}

void set_block_Matrix_T
	(struct Matrix_T* a, const struct const_Matrix_T* a_sub, const ptrdiff_t row0, const ptrdiff_t col0,
	 const char set_type)
{
	assert(a->layout == a_sub->layout); // Add support if required.
	assert(a->layout == 'R');

	set_value_fptr_T set_value = NULL;
	switch (set_type) {
		case 'i': set_value = set_value_insert_T; break;
		case 'a': set_value = set_value_add_T;    break;
		default:  EXIT_ERROR("Unsupported: %c.\n",set_type); break;
	}

	const ptrdiff_t ext_0 = a_sub->ext_0,
	                ext_1 = a_sub->ext_1;

	assert(row0+ext_0 <= a->ext_0);
	assert(col0+ext_1 <= a->ext_1);

	for (int i = 0, row = (int)row0; i < ext_0; ++i, ++row) {
		const Type*const data_as = get_row_const_Matrix_T(i,a_sub);

		Type* data_a = get_row_Matrix_T(row,a);
		data_a += col0;
		for (int j = 0; j < ext_1; ++j)
			set_value(&data_a[j],data_as[j]);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

void set_value_insert_T (Type*const dest, const Type src)
{
	*dest = src;
}

void set_value_add_T (Type*const dest, const Type src)
{
	*dest += src;
}
