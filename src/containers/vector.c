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

#include "vector.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"

#include "multiarray.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

/** \brief Comparison function for std::qsort between `int*` `a` and `b`.
 *	\return a - b.  */
static int cmp_i
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

// Interface functions ********************************************************************************************** //

void reorder_Vector_i (struct Vector_i*const a, const int*const ordering)
{
	const ptrdiff_t size = a->ext_0;

	int b[size];
	for (ptrdiff_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (ptrdiff_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

void resize_Vector_i (struct Vector_i*const a, const ptrdiff_t ext_0)
{
	const ptrdiff_t size_i = a->ext_0;
	a->ext_0 = ext_0;
	const ptrdiff_t size_o = a->ext_0;

	if (size_o <= size_i)
		return;

	const int* data_i = a->data;
	a->data = malloc(size_o * sizeof *(a->data)); // keep
	if (size_i != 0) {
		for (ptrdiff_t i = 0; i < size_i; i++)
			a->data[i] = data_i[i];
	}
	free((void*)data_i);
}

void set_to_zero_Vector_i (struct Vector_i*const a)
{
	const ptrdiff_t i_max = a->ext_0;
	for (ptrdiff_t i = 0; i < i_max; i++)
		a->data[i] = 0;
}

void set_to_data_Vector_i (struct Vector_i*const a, const int*const data_src)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] = data_src[i];
}

void set_to_value_Vector_d (struct Vector_d*const a, const double val)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] = val;
}

void sort_Vector_i (struct Vector_i* a)
{
	const ptrdiff_t size = a->ext_0;
	qsort(a->data,size,sizeof(a->data[0]),cmp_i);
}

int sum_Vector_i (const struct Vector_i* a)
{
	int sum = 0;

	const ptrdiff_t size = a->ext_0;
	for (ptrdiff_t i = 0; i < size; ++i)
		sum += a->data[i];
	return sum;
}

ptrdiff_t prod_Vector_i (const struct Vector_i* a)
{
	const ptrdiff_t i_max = a->ext_0;

	ptrdiff_t prod = 1;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		prod *= a->data[i];

	return prod;
}

ptrdiff_t prod_const_Vector_i (const struct const_Vector_i* a)
{
	return prod_Vector_i ((const struct Vector_i*)a);
}

bool check_equal_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = a->ext_0;
	if (size != b->ext_0)
		return false;

	const int* data_a = a->data,
	         * data_b = b->data;

	for (ptrdiff_t i = 0; i < size; i++) {
		if (*data_a++ != *data_b++)
			return false;
	}
	return true;
}

bool check_equal_Vector_i_i (const struct Vector_i*const a, const int* data_b)
{
	const int* data_a = a->data;

	const ptrdiff_t size = a->ext_0;
	for (ptrdiff_t i = 0; i < size; i++) {
		if (*data_a++ != *data_b++)
			return false;
	}
	return true;
}

int cmp_Vector_i (const void *a, const void *b)
{
	const struct Vector_i*const*const ia = (const struct Vector_i*const*const) a,
	                     *const*const ib = (const struct Vector_i*const*const) b;

	const ptrdiff_t size_a = (*ia)->ext_0,
	                size_b = (*ib)->ext_0;

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const int*const data_a = (*ia)->data,
	         *const data_b = (*ib)->data;

	for (ptrdiff_t i = 0; i < size_a; ++i) {
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
	}
	return 0;
}

void copy_data_Vector_i_Vector_i (const struct Vector_i*const src, struct Vector_i*const dest)
{
	const ptrdiff_t size_src  = src->ext_0,
	                size_dest = dest->ext_0;

	if (size_src != size_dest)
		EXIT_UNSUPPORTED;

	for (ptrdiff_t i = 0; i < size_src; ++i)
		dest->data[i] = src->data[i];
}

void push_back_Vector_i (struct Vector_i*const src, const int val, const bool sorted, const bool unique)
{
	if (sorted)
		sort_Vector_i(src);

	const bool add_val = unique ? !find_val_Vector_i((struct const_Vector_i*)src,val,sorted) : true;

	if (!add_val)
		return;

	resize_Vector_i(src,src->ext_0+1);
	src->data[src->ext_0-1] = val;

	if (sorted)
		sort_Vector_i(src);
}

bool find_val_Vector_i (const struct const_Vector_i*const src, const int val, const bool sorted)
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
		const int* ind_ptr = bsearch(&val,src->data,src->ext_0,sizeof(src->data[0]),cmp_i);
		if (ind_ptr)
			found = true;
	}
	return found;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static int cmp_i (const void *a, const void *b)
{
	return (int) ( *(int*)a - *(int*)b );
}
