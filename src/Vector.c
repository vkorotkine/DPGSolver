// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "Vector.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "allocators.h"
#include "Multiarray.h"

// Add number indicating level of dereferencing for constructors here and for Matrix.

static struct Vector_ui make_local_Vector_ui_0 (const size_t ext_0, const bool owns_data, unsigned int*const data);

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

	struct Vector_ui local = { .extents[0] = ext_0,
	                           .owns_data  = true,
	                           .data       = data,
	                         };

	struct Vector_ui* dest = malloc(sizeof *dest); // returned
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

void destructor_Vector_ui_2 (struct Vector_ui** a, const size_t n_src)
{
	for (size_t n = 0; n < n_src; n++)
		destructor_Vector_ui(a[n]);
	free(a);
}

// Helper functions ************************************************************************************************* //

void reorder_Vector_ui (struct Vector_ui*const a, unsigned int*const ordering)
{
	const size_t size = a->extents[0];

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
