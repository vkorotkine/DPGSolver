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
#include "definitions_tol.h"

#include "def_templates_matrix.h"
// Static function declarations ************************************************************************************* //

/** \brief Pointer to value setting functions with scaling.
 *
 *  \param alpha Scaling constant.
 *  \param dest  The destination.
 *  \param src   The source.
 */
typedef void (*set_scaled_value_fptr_T)
	(const Type alpha,
	 Type*const dest,
	 const Type src
	);

/// \brief Version of \ref set_scaled_value_fptr_T inserting values.
static void set_scaled_value_insert_T
	(const Type alpha, ///< See brief.
	 Type*const dest,  ///< See brief.
	 const Type src    ///< See brief.
	);

/// \brief Version of \ref set_scaled_value_fptr_T adding values.
static void set_scaled_value_add_T
	(const Type alpha, ///< See brief.
	 Type*const dest,  ///< See brief.
	 const Type src    ///< See brief.
	);

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
	(struct Matrix_T* dest, const ptrdiff_t row0_d, const ptrdiff_t col0_d, const struct const_Matrix_T* src,
	 const ptrdiff_t row0_s, const ptrdiff_t col0_s, const ptrdiff_t ext_0, const ptrdiff_t ext_1,
	 const char set_type)
{
	assert(dest->layout == src->layout); // Add support if required.
	const char layout = dest->layout;

	set_value_fptr_T set_value = NULL;
	switch (set_type) {
		case 'i': set_value = set_value_insert_T; break;
		case 'a': set_value = set_value_add_T;    break;
		default:  EXIT_ERROR("Unsupported: %c.\n",set_type); break;
	}

	assert(row0_d+ext_0 <= dest->ext_0);
	assert(col0_d+ext_1 <= dest->ext_1);

	assert(row0_s+ext_0 <= src->ext_0);
	assert(col0_s+ext_1 <= src->ext_1);

	if (layout == 'R') {
		ptrdiff_t row_d = row0_d,
		          row_s = row0_s;
		for (int i = 0; i < ext_0; ++i, ++row_d, ++row_s) {
			const Type*const data_s = get_row_const_Matrix_T(row_s,src)+col0_s;
			Type*const data_d       = get_row_Matrix_T(row_d,dest)+col0_d;
			for (int j = 0; j < ext_1; ++j)
				set_value(&data_d[j],data_s[j]);
		}
	} else if (layout == 'C') {
		ptrdiff_t col_d = col0_d,
		          col_s = col0_s;
		for (int j = 0; j < ext_1; ++j, ++col_d, ++col_s) {
			const Type*const data_s = get_col_const_Matrix_T(col_s,src)+row0_s;
			Type*const data_d       = get_col_Matrix_T(col_d,dest)+row0_d;
			for (int i = 0; i < ext_0; ++i)
				set_value(&data_d[i],data_s[i]);
		}
	} else {
		EXIT_ERROR("Unsupported: %c\n",layout);
	}
}

#if TYPE_RC == TYPE_COMPLEX
void set_block_Matrix_T_R
	(struct Matrix_T* dest, const ptrdiff_t row0_d, const ptrdiff_t col0_d, const struct const_Matrix_R* src,
	 const ptrdiff_t row0_s, const ptrdiff_t col0_s, const ptrdiff_t ext_0, const ptrdiff_t ext_1,
	 const char set_type)
{
	const struct const_Matrix_T*const src_c = constructor_copy_const_Matrix_T_Matrix_R(src); // destructed
	set_block_Matrix_T(dest,row0_d,col0_d,src_c,row0_s,col0_s,ext_0,ext_1,set_type);
	destructor_const_Matrix_T(src_c);
}

void set_block_Matrix_R_cmplx_step
	(struct Matrix_R* dest, const struct const_Matrix_C* src, const ptrdiff_t row0_d, const ptrdiff_t col0_d,
	 const char set_type)
{
	assert(dest->layout == src->layout); // Add support if required.
	const char layout = dest->layout;

	set_value_fptr_d set_value = NULL;
	switch (set_type) {
		case 'i': set_value = set_value_insert_d; break;
		case 'a': set_value = set_value_add_d;    break;
		default:  EXIT_ERROR("Unsupported: %c.\n",set_type); break;
	}

	const ptrdiff_t ext_0 = src->ext_0,
	                ext_1 = src->ext_1;

	assert(row0_d+ext_0 <= dest->ext_0);
	assert(col0_d+ext_1 <= dest->ext_1);

	if (layout == 'R') {
		for (int i = 0, row = (int)row0_d; i < ext_0; ++i, ++row) {
			const Complex*const data_s = get_row_const_Matrix_C(i,src);
			Real*const data_d        = get_row_Matrix_R(row,dest)+col0_d;
			for (int j = 0; j < ext_1; ++j)
				set_value(&data_d[j],cimag(data_s[j]/CX_STEP));
		}
	} else if (layout == 'C') {
		for (int j = 0, col = (int)col0_d; j < ext_1; ++j, ++col) {
			const Complex*const data_s = get_col_const_Matrix_C(j,src);
			Real*const data_d        = get_col_Matrix_R(col,dest)+row0_d;
			for (int i = 0; i < ext_0; ++i)
				set_value(&data_d[i],cimag(data_s[i]/CX_STEP));
		}
	} else {
		EXIT_ERROR("Unsupported: %c\n",layout);
	}
}
#endif

void set_scaled_block_Matrix_T
	(const Type alpha, struct Matrix_T* dest, const ptrdiff_t row0_d, const ptrdiff_t col0_d,
	 const struct const_Matrix_T* src, const ptrdiff_t row0_s, const ptrdiff_t col0_s, const ptrdiff_t ext_0,
	 const ptrdiff_t ext_1, const char set_type)
{
	assert(dest->layout == src->layout); // Add support if required.
	const char layout = dest->layout;

	set_scaled_value_fptr_T set_scaled_value = NULL;
	switch (set_type) {
		case 'i': set_scaled_value = set_scaled_value_insert_T; break;
		case 'a': set_scaled_value = set_scaled_value_add_T;    break;
		default:  EXIT_ERROR("Unsupported: %c.\n",set_type); break;
	}

	assert(row0_d+ext_0 <= dest->ext_0);
	assert(col0_d+ext_1 <= dest->ext_1);

	assert(row0_s+ext_0 <= src->ext_0);
	assert(col0_s+ext_1 <= src->ext_1);

	if (layout == 'R') {
		ptrdiff_t row_d = row0_d,
		          row_s = row0_s;
		for (int i = 0; i < ext_0; ++i, ++row_d, ++row_s) {
			const Type*const data_s = get_row_const_Matrix_T(row_s,src)+col0_s;
			Type*const data_d       = get_row_Matrix_T(row_d,dest)+col0_d;
			for (int j = 0; j < ext_1; ++j)
				set_scaled_value(alpha,&data_d[j],data_s[j]);
		}
	} else if (layout == 'C') {
		ptrdiff_t col_d = col0_d,
		          col_s = col0_s;
		for (int j = 0; j < ext_1; ++j, ++col_d, ++col_s) {
			const Type*const data_s = get_col_const_Matrix_T(col_s,src)+row0_s;
			Type*const data_d       = get_col_Matrix_T(col_d,dest)+row0_d;
			for (int i = 0; i < ext_0; ++i)
				set_scaled_value(alpha,&data_d[i],data_s[i]);
		}
	} else {
		EXIT_ERROR("Unsupported: %c\n",layout);
	}
}

void set_value_insert_T (Type*const dest, const Type src)
{
	*dest = src;
}

void set_value_add_T (Type*const dest, const Type src)
{
	*dest += src;
}

void swap_rows_Matrix_T (struct Matrix_T*const src, const int r0, const int r1)
{
	if (r0 == r1)
		return;

	const ptrdiff_t ext_1 = src->ext_1;
	Type* data_r0 = get_row_Matrix_T(r0,src);
	Type* data_r1 = get_row_Matrix_T(r1,src);
	for (int i = 0; i < ext_1; ++i) {
		const Type tmp = *data_r0;
		*data_r0++ = *data_r1;
		*data_r1++ = tmp;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_scaled_value_insert_T (const Type alpha, Type*const dest, const Type src)
{
	*dest = alpha*src;
}

static void set_scaled_value_add_T (const Type alpha, Type*const dest, const Type src)
{
	*dest += alpha*src;
}

#include "undef_templates_matrix.h"
