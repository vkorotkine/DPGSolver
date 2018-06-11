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

#include <assert.h>

#include "macros.h"

#include "def_templates_vector.h"

// Static function declarations ************************************************************************************* //

/** \brief Comparison function for std::qsort between `Type*` `a` and `b`.
 *  \return Standard required return for comparator function of qsort. */
static int cmp_T
	(const void* a, ///< Variable 1.
	 const void* b  ///< Variable 2.
	);

// Interface functions ********************************************************************************************** //

void reorder_Vector_T (struct Vector_T*const a, const int*const ordering)
{
	const ptrdiff_t size = a->ext_0;

	Type b[size];
	for (ptrdiff_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (ptrdiff_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

void resize_Vector_T (struct Vector_T*const a, const ptrdiff_t ext_0)
{
	const ptrdiff_t size_i = a->ext_0;
	a->ext_0 = ext_0;
	const ptrdiff_t size_o = a->ext_0;

	if (size_o <= size_i)
		return;

	const Type* data_i = a->data;
	a->data = malloc((size_t)size_o * sizeof *(a->data)); // keep
	if (size_i != 0) {
		for (ptrdiff_t i = 0; i < size_i; i++)
			a->data[i] = data_i[i];
	}
	free((void*)data_i);
}

void set_to_zero_Vector_T (struct Vector_T*const a)
{
	const ptrdiff_t i_max = a->ext_0;
	for (ptrdiff_t i = 0; i < i_max; i++)
		a->data[i] = 0;
}

void set_to_data_Vector_T (struct Vector_T*const a, const Type*const data_src)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] = data_src[i];
}

void set_to_Vector_Vector_T (struct Vector_T*const dest, const Type alpha, const struct const_Vector_T*const src)
{
	const ptrdiff_t ext_0 = dest->ext_0;
	assert(ext_0 == src->ext_0);
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		dest->data[i] = alpha*src->data[i];
}

void set_to_value_Vector_T (struct Vector_T*const a, const Type val)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] = val;
}

void sort_Vector_T (struct Vector_T* a)
{
	const ptrdiff_t size = a->ext_0;
	qsort(a->data,(size_t)size,sizeof(a->data[0]),cmp_T);
}

Type sum_Vector_T (const struct Vector_T* a)
{
	Type sum = 0;

	const ptrdiff_t size = a->ext_0;
	for (ptrdiff_t i = 0; i < size; ++i)
		sum += a->data[i];
	return sum;
}

Type prod_Vector_T (const struct Vector_T* a)
{
	const ptrdiff_t i_max = a->ext_0;

	Type prod = 1;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		prod *= a->data[i];

	return prod;
}

Type prod_const_Vector_T (const struct const_Vector_T* a)
{
	return prod_Vector_T ((const struct Vector_T*)a);
}

bool check_equal_Vector_T (const struct Vector_T*const a, const struct Vector_T*const b)
{
	const ptrdiff_t size = a->ext_0;
	if (size != b->ext_0)
		return false;

	const Type* data_a = a->data,
	         * data_b = b->data;

	for (ptrdiff_t i = 0; i < size; i++) {
		if (*data_a++ != *data_b++)
			return false;
	}
	return true;
}

bool check_equal_Vector_T_T (const struct Vector_T*const a, const Type* data_b)
{
	const Type* data_a = a->data;

	const ptrdiff_t size = a->ext_0;
	for (ptrdiff_t i = 0; i < size; i++) {
		if (*data_a++ != *data_b++)
			return false;
	}
	return true;
}

int cmp_Vector_T (const void *a, const void *b)
{
	const struct Vector_T*const*const ia = (const struct Vector_T*const*const) a,
	                     *const*const ib = (const struct Vector_T*const*const) b;

	const ptrdiff_t size_a = (*ia)->ext_0,
	                size_b = (*ib)->ext_0;

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const Type*const data_a = (*ia)->data,
	         *const data_b = (*ib)->data;

	for (ptrdiff_t i = 0; i < size_a; ++i) {
#if TYPE_RC == TYPE_COMPLEX
		if (creal(data_a[i]) > creal(data_b[i]))
			return 1;
		else if (creal(data_a[i]) < creal(data_b[i]))
			return -1;
#else
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
#endif
	}
	return 0;
}

void copy_data_Vector_T_Vector_T (const struct Vector_T*const src, struct Vector_T*const dest)
{
	const ptrdiff_t size_src  = src->ext_0,
	                size_dest = dest->ext_0;

	if (size_src != size_dest)
		EXIT_UNSUPPORTED;

	for (ptrdiff_t i = 0; i < size_src; ++i)
		dest->data[i] = src->data[i];
}

void push_back_Vector_T (struct Vector_T*const src, const Type val, const bool sorted, const bool unique)
{
	if (sorted)
		sort_Vector_T(src);

	const bool add_val = unique ? !find_val_Vector_T((struct const_Vector_T*)src,val,sorted) : true;

	if (!add_val)
		return;

	resize_Vector_T(src,src->ext_0+1);
	src->data[src->ext_0-1] = val;

	if (sorted)
		sort_Vector_T(src);
}

void push_back_Vector_Vector_T (struct Vector_T*const src, const struct const_Vector_T*const vals)
{
	const ptrdiff_t ext_0_src = src->ext_0,
	                ext_0_add = vals->ext_0,
	                ext_0_sum = ext_0_src+ext_0_add;
	resize_Vector_T(src,ext_0_sum);
	for (ptrdiff_t i = ext_0_src, j = 0; i < ext_0_sum; ++i, ++j)
		src->data[i] = vals->data[j];
}

bool find_val_Vector_T (const struct const_Vector_T*const src, const Type val, const bool sorted)
{
	bool found = false;
	if (!sorted) {
		const ptrdiff_t i_max = src->ext_0;
		for (ptrdiff_t i = 0; i < i_max; ++i) {
			if (src->data[i] == val) {
				found = true;
				break;
			}
		}
	} else {
		const int* ind_ptr = bsearch(&val,src->data,(size_t)src->ext_0,sizeof(src->data[0]),cmp_T);
		if (ind_ptr)
			found = true;
	}
	return found;
}

void swap_vals_Vector_T (struct Vector_T*const src, const int r0, const int r1)
{
	const Type tmp = src->data[r0];
	src->data[r0] = src->data[r1];
	src->data[r1] = tmp;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static int cmp_T (const void* a, const void* b)
{
	return (int) ( *(Type*)a - *(Type*)b );
}

#include "undef_templates_vector.h"
