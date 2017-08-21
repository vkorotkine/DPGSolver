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

/** \brief Make a local copy of a \ref Vector_ui (static memory).
 *	\return The \ref Vector_ui.
 */
static struct Vector_ui make_local_Vector_ui_0
	(const size_t ext_0,     ///< Standard.
	 const bool owns_data,   ///< Standard.
	 unsigned int*const data ///< Standard.
	);

/** \brief Comparison function for std::qsort between `unsigned int*` `a` and `b`.
 *	\return a - b.
 */
static int cmp_ui
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

// Constructor/Destructor functions ********************************************************************************* //

struct Vector_ui* constructor_default_Vector_ui ()
{
	struct Vector_ui* dest = malloc(sizeof *dest); // returned
	dest->extents[0] = 0;
	dest->owns_data  = true;
	dest->data       = NULL;

	return dest;
}

struct Vector_ui** constructor_default_Vector_ui_2 (const size_t n_dest)
{
	struct Vector_ui** dest = malloc(n_dest * sizeof *dest); // returned;

	for (size_t n = 0; n < n_dest; n++)
		dest[n] = constructor_default_Vector_ui();

	return dest;
}

struct Vector_ui* constructor_empty_Vector_ui (const size_t ext_0)
{
	unsigned int* data = mallocator(UINT_T,1,ext_0); // keep

	struct Vector_ui* dest = malloc(sizeof *dest); // returned

	dest->extents[0] = ext_0;
	dest->owns_data  = true;
	dest->data       = data;

	return dest;
}

struct Vector_ui* constructor_copy_Vector_ui (const struct Vector_ui*const src)
{
	const size_t ext_0 = src->extents[0];
	const unsigned int*const data_src = src->data;

	unsigned int* data = malloc(ext_0 * sizeof *data); // keep
	for (size_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	struct Vector_ui* dest = malloc(sizeof *dest); // returned

	dest->extents[0] = ext_0;
	dest->owns_data  = true;
	dest->data       = data;

	return dest;
}

struct Vector_ui* constructor_copy_Vector_ui_ui (const size_t ext_0, const unsigned int*const data_src)
{
	unsigned int* data = malloc(ext_0 * sizeof *data); // keep
	for (size_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	struct Vector_ui* dest = malloc(sizeof *dest); // returned
	dest->extents[0] = ext_0;
	dest->owns_data  = true;
	dest->data       = data;

	return dest;
}

struct Vector_ui* constructor_move_Vector_ui_ui (const size_t ext_0, const bool owns_data, unsigned int*const data)
{
	struct Vector_ui* dest = malloc(sizeof *dest); // returned
	dest->extents[0] = ext_0;
	dest->owns_data  = owns_data;
	dest->data       = data;

	return dest;
}

struct const_Vector_ui* constructor_move_const_Vector_ui_ui
	(const size_t ext_0, const bool owns_data, const unsigned int*const data)
{
	struct Vector_ui local = make_local_Vector_ui_0(ext_0,owns_data,(unsigned int*)data);

	struct const_Vector_ui* dest = malloc(sizeof *dest); // returned
	memcpy(dest,&local,sizeof *dest);

	return dest;
}

void const_constructor_move_Vector_ui (const struct const_Vector_ui*const* dest, struct Vector_ui* src)
{
	*(struct const_Vector_ui**) dest = (struct const_Vector_ui*) src;
}

void destructor_Vector_ui (struct Vector_ui* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_Vector_ui_2 (struct Vector_ui** a, const size_t n_src, const bool owns_data)
{
	if (owns_data) {
		for (size_t n = 0; n < n_src; n++)
			destructor_Vector_ui(a[n]);
	}
	free(a);
}

// Helper functions ************************************************************************************************* //

void reorder_Vector_ui (struct Vector_ui*const a, const unsigned int*const ordering)
{
	const size_t size = compute_size(1,a->extents);

	unsigned int b[size];
	for (size_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (size_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

void reserve_Vector_ui (struct Vector_ui*const a, const size_t ext_0)
{
	const size_t size_i = compute_size(1,a->extents);
	a->extents[0] = ext_0;
	const size_t size_o = compute_size(1,a->extents);

	if (size_o <= size_i)
		return;

	const unsigned int* data_i = a->data;
	a->data = malloc(size_o * sizeof *(a->data)); // keep
	if (size_i != 0) {
		for (size_t i = 0; i < size_i; i++)
			a->data[i] = data_i[i];
	}
}

void set_to_zero_Vector_ui (struct Vector_ui*const a)
{
	for (size_t i = 0; i < a->extents[0]; i++)
		a->data[i] = 0;
}

void set_to_data_Vector_ui (struct Vector_ui*const a, const unsigned int*const data_src)
{
	const size_t ext_0 = compute_size(1,a->extents);
	for (size_t i = 0; i < ext_0; ++i)
		a->data[i] = data_src[i];
}

void sort_Vector_ui (struct Vector_ui* a)
{
	const size_t size = compute_size(1,a->extents);
	qsort(a->data,size,sizeof(a->data[0]),cmp_ui);
}

unsigned int sum_Vector_ui (struct Vector_ui* a)
{
	unsigned int sum = 0;

	const size_t size = compute_size(1,a->extents);
	for (size_t i = 0; i < size; ++i)
		sum += a->data[i];
	return sum;
}

bool check_equal_Vector_ui (const struct Vector_ui*const a, const struct Vector_ui*const b)
{
	const size_t size = compute_size(1,a->extents);
	if (size != compute_size(1,b->extents))
		return false;

	const unsigned int* data_a = a->data,
	                  * data_b = b->data;

	for (size_t i = 0; i < size; i++) {
		if (*data_a++ != *data_b++)
			return false;
	}
	return true;
}

int cmp_Vector_ui (const void *a, const void *b)
{
	const struct Vector_ui*const*const ia = (const struct Vector_ui*const*const) a,
	                      *const*const ib = (const struct Vector_ui*const*const) b;

	const size_t size_a = compute_size(1,(*ia)->extents),
	             size_b = compute_size(1,(*ib)->extents);

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const unsigned int*const data_a = (*ia)->data,
	                  *const data_b = (*ib)->data;

	for (size_t i = 0; i < size_a; ++i) {
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
	}
	return 0;
}

void copy_data_Vector_ui_Vector_ui (const struct Vector_ui*const src, struct Vector_ui*const dest)
{
	const size_t size_src  = compute_size(1,src->extents),
	             size_dest = compute_size(1,dest->extents);

	if (size_src != size_dest)
		EXIT_UNSUPPORTED;

	for (size_t i = 0; i < size_src; ++i)
		dest->data[i] = src->data[i];
}

// Printing functions *********************************************************************************************** //

void print_Vector_ui (const struct Vector_ui*const a)
{
	const size_t ext = a->extents[0];

	const unsigned int* data = a->data;

	for (size_t i = 0; i < ext; i++) {
		printf("% 12d ",*data++);
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

void print_const_Vector_ui (const struct const_Vector_ui*const a)
{
	struct Vector_ui local = make_local_Vector_ui_0(a->extents[0],false,(unsigned int*)a->data);
	print_Vector_ui(&local);
}


// Static functions ************************************************************************************************* //

static struct Vector_ui make_local_Vector_ui_0 (const size_t ext_0, const bool owns_data, unsigned int*const data)
{
	struct Vector_ui local =
		{ .extents[0] = ext_0,
		  .owns_data  = owns_data,
		  .data       = data,
		};
	return local;
}

static int cmp_ui (const void *a, const void *b)
{
	return (int) ( *(unsigned int*)a - *(unsigned int*)b );
}
