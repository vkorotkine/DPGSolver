// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "Vector.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "allocators.h"

// Add number indicating level of dereferencing for constructors here and for Matrix.

static struct Vector_ui make_local_Vector_ui_0 (const size_t ext_0, const bool owns_data, unsigned int*const data);

struct Vector_ui* constructor_default_Vector_ui_1 ()
{
	struct Vector_ui* dest = malloc(sizeof *dest); // returned
	dest->owns_data = false;

	return dest;
}

struct Vector_ui** constructor_default_Vector_ui_2 (const size_t n_dest)
{
	struct Vector_ui** dest = malloc(n_dest * sizeof *dest); // returned;

	for (size_t n = 0; n < n_dest; n++)
		dest[n] = constructor_default_Vector_ui_1();

	return dest;
}

//	struct Vector_ui** data = constructor_default_Vector_ui_2(compute_size(order,extents)); // keep

struct Vector_ui* constructor_empty_Vector_ui_1 (const size_t ext_0)
{
	unsigned int* data = mallocator(UINT_T,1,ext_0); // keep

	struct Vector_ui local = { .extents[0] = ext_0,
	                           .owns_data  = true,
	                           .data       = data,
	                         };

	struct Vector_ui* dest = malloc(sizeof *dest); // returned
	memcpy(dest,&local,sizeof *dest);

	return dest;
}

struct const_Vector_ui* constructor_move_const_Vector_ui_1_Vector_ui (struct Vector_ui*const src)
{
	src->owns_data = false;

	struct Vector_ui local = make_local_Vector_ui_0(src->extents[0],true,src->data);

	struct const_Vector_ui* dest = malloc(sizeof *dest); // returned
	memcpy(dest,&local,sizeof *dest);

	return dest;
}

void const_constructor_const_Vector_ui_1_Vector_ui (const struct const_Vector_ui*const* dest, struct Vector_ui* src)
{
	struct const_Vector_ui* local = constructor_move_const_Vector_ui_1_Vector_ui(src); // keep
	*(struct const_Vector_ui**) dest = local;
}

void reorder_Vector_ui (struct Vector_ui*const a, unsigned int*const ordering)
{
	const size_t size = a->extents[0];

	unsigned int b[size];
	for (size_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (size_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

void print_Vector_ui (const struct Vector_ui*const a)
{
	const size_t ext = a->extents[0];

	const unsigned int* data = a->data;

	for (size_t i = 0; i < ext; i++)
		printf("% 12d ",*data++);
	printf("\n");
}

// Static functions ************************************************************************************************* //

/// \brief Make a local copy of a \ref Vector_ui (static memory).
static struct Vector_ui make_local_Vector_ui_0
	(const size_t ext_0,     ///< Standard.
	 const bool owns_data,   ///< Standard.
	 unsigned int*const data ///< Standard.
	)
{
	struct Vector_ui local = { .extents[0] = ext_0,
	                           .owns_data  = owns_data,
	                           .data       = data,
	                         };
	return local;
}
