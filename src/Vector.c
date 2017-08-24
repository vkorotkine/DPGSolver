// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "Vector.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Macros.h"

#include "allocators.h"
#include "Multiarray.h"

// Static function declarations ************************************************************************************* //

/** \brief Make a local \ref Vector_i\* (dynamic memory).
 *	\return See brief.  */
static struct Vector_i* constructor_local_Vector_i_1
	(const ptrdiff_t ext_0, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

/** \brief Comparison function for std::qsort between `int*` `a` and `b`.
 *	\return a - b.  */
static int cmp_i
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

// Constructor/Destructor functions ********************************************************************************* //

struct Vector_i* constructor_default_Vector_i ()
{
	struct Vector_i* dest = malloc(sizeof *dest); // returned
	dest->extents[0] = 0;
	dest->owns_data  = true;
	dest->data       = NULL;

	return dest;
}

struct Vector_i** constructor_default_Vector_i_2 (const ptrdiff_t n_dest)
{
	struct Vector_i** dest = malloc(n_dest * sizeof *dest); // returned;

	for (ptrdiff_t n = 0; n < n_dest; n++)
		dest[n] = constructor_default_Vector_i();

	return dest;
}

struct Vector_i* constructor_empty_Vector_i (const ptrdiff_t ext_0)
{
	int* data = mallocator(UINT_T,1,ext_0); // keep

	struct Vector_i* dest = malloc(sizeof *dest); // returned

	dest->extents[0] = ext_0;
	dest->owns_data  = true;
	dest->data       = data;

	return dest;
}

struct Vector_i* constructor_copy_Vector_i (const struct Vector_i*const src)
{
	const ptrdiff_t ext_0 = src->extents[0];
	const int*const data_src = src->data;

	int* data = malloc(ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	struct Vector_i* dest = malloc(sizeof *dest); // returned

	dest->extents[0] = ext_0;
	dest->owns_data  = true;
	dest->data       = data;

	return dest;
}

struct Vector_i* constructor_copy_Vector_i_i (const ptrdiff_t ext_0, const int*const data_src)
{
	int* data = malloc(ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	struct Vector_i* dest = malloc(sizeof *dest); // returned
	dest->extents[0] = ext_0;
	dest->owns_data  = true;
	dest->data       = data;

	return dest;
}

struct Vector_i* constructor_move_Vector_i_i (const ptrdiff_t ext_0, const bool owns_data, int*const data)
{
	struct Vector_i* dest = malloc(sizeof *dest); // returned
	dest->extents[0] = ext_0;
	dest->owns_data  = owns_data;
	dest->data       = data;

	return dest;
}

struct const_Vector_i* constructor_move_const_Vector_i_i
	(const ptrdiff_t ext_0, const bool owns_data, const int*const data)
{
	struct Vector_i* local = constructor_local_Vector_i_1(ext_0,owns_data,(int*)data); // free

	struct const_Vector_i* dest = malloc(sizeof *dest); // returned
	memcpy(dest,local,sizeof *dest);
	free(local);

	return dest;
}

void const_constructor_move_Vector_i (const struct const_Vector_i*const* dest, struct Vector_i* src)
{
	*(struct const_Vector_i**) dest = (struct const_Vector_i*) src;
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

// Helper functions ************************************************************************************************* //

void reorder_Vector_i (struct Vector_i*const a, const int*const ordering)
{
	const ptrdiff_t size = compute_size(1,a->extents);

	int b[size];
	for (ptrdiff_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (ptrdiff_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

void reserve_Vector_i (struct Vector_i*const a, const ptrdiff_t ext_0)
{
	const ptrdiff_t size_i = compute_size(1,a->extents);
	a->extents[0] = ext_0;
	const ptrdiff_t size_o = compute_size(1,a->extents);

	if (size_o <= size_i)
		return;

	const int* data_i = a->data;
	a->data = malloc(size_o * sizeof *(a->data)); // keep
	if (size_i != 0) {
		for (ptrdiff_t i = 0; i < size_i; i++)
			a->data[i] = data_i[i];
	}
}

void set_to_zero_Vector_i (struct Vector_i*const a)
{
	for (ptrdiff_t i = 0; i < a->extents[0]; i++)
		a->data[i] = 0;
}

void set_to_data_Vector_i (struct Vector_i*const a, const int*const data_src)
{
	const ptrdiff_t ext_0 = compute_size(1,a->extents);
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] = data_src[i];
}

void sort_Vector_i (struct Vector_i* a)
{
	const ptrdiff_t size = compute_size(1,a->extents);
	qsort(a->data,size,sizeof(a->data[0]),cmp_i);
}

int sum_Vector_i (struct Vector_i* a)
{
	int sum = 0;

	const ptrdiff_t size = compute_size(1,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		sum += a->data[i];
	return sum;
}

bool check_equal_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = compute_size(1,a->extents);
	if (size != compute_size(1,b->extents))
		return false;

	const int* data_a = a->data,
	         * data_b = b->data;

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

	const ptrdiff_t size_a = compute_size(1,(*ia)->extents),
	                size_b = compute_size(1,(*ib)->extents);

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
	const ptrdiff_t size_src  = compute_size(1,src->extents),
	                size_dest = compute_size(1,dest->extents);

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

	reserve_Vector_i(src,src->extents[0]+1);
	src->data[src->extents[0]-1] = val;

	if (sorted)
		sort_Vector_i(src);
}

bool find_val_Vector_i (const struct const_Vector_i*const src, const int val, const bool sorted)
{
	bool found = false;
	if (!sorted) {
		const ptrdiff_t i_max = 0;
		for (ptrdiff_t i = 0; i < i_max; ++i) {
			if (src->data[i] == val) {
				found = true;
				break;
			}
		}
	} else {
		const int* ind_ptr = bsearch(&val,src->data,src->extents[0],sizeof(src->data[0]),cmp_i);
		if (ind_ptr)
			found = true;
	}
	return found;
}

// Printing functions *********************************************************************************************** //

void print_Vector_i (const struct Vector_i*const a)
{
	const ptrdiff_t ext = a->extents[0];

	const int* data = a->data;

	for (ptrdiff_t i = 0; i < ext; i++) {
		printf("% 12d ",*data++);
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

void print_const_Vector_i (const struct const_Vector_i*const a)
{
	struct Vector_i* local = constructor_local_Vector_i_1(a->extents[0],false,(int*)a->data); // free
	print_Vector_i(local);
	free(local);
}


// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Vector_i* constructor_local_Vector_i_1 (const ptrdiff_t ext_0, const bool owns_data, int*const data)
{
	struct Vector_i* dest = malloc(sizeof *dest); // returned

	dest->extents[0] = ext_0;
	dest->owns_data  = owns_data;
	dest->data       = data;

	return dest;
}

static int cmp_i (const void *a, const void *b)
{
	return (int) ( *(int*)a - *(int*)b );
}
