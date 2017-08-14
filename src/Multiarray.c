// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "Multiarray.h"

#include <string.h>

#include "Macros.h"

#include "allocators.h"
#include "Vector.h"

static struct Multiarray_Vector_ui make_local_Multiarray_Vector_ui_0
	(const char layout, const size_t order, size_t*const extents, const bool owns_data, struct Vector_ui**const data);

// Constructor/Destructor functions ********************************************************************************* //

struct Multiarray_d* constructor_move_Multiarray_d_d (const char layout, double*const data, const size_t order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	size_t* extents = set_extents(order,ap); // keep
	va_end(ap);

	struct Multiarray_d local = { .layout    = layout,
	                              .order     = order,
	                              .extents   = extents,
	                              .owns_data = false,
	                              .data      = data,
	                            };

	struct Multiarray_d* A = malloc(sizeof *A); // returned
	memcpy(A,&local,sizeof *A);

	return A;
}

struct Multiarray_Vector_ui* constructor_empty_Multiarray_Vector_ui (const size_t order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	size_t* extents = set_extents(order,ap); // keep
	va_end(ap);

	struct Vector_ui** data = constructor_default_Vector_ui_2(compute_size(order,extents)); // keep

	struct Multiarray_Vector_ui local = make_local_Multiarray_Vector_ui_0('R',order,extents,true,data);

	struct Multiarray_Vector_ui* dest = malloc(sizeof *dest); // returned
	memcpy(dest,&local,sizeof *dest);

	return dest;
}

void const_constructor_move_Multiarray_Vector_ui
	(const struct const_Multiarray_Vector_ui*const* dest, struct Multiarray_Vector_ui* src)
{
	*(struct const_Multiarray_Vector_ui**) dest = (struct const_Multiarray_Vector_ui*) src;
}

void destructor_Multiarray_d (struct Multiarray_d* a)
{
	free(a->extents);
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_Multiarray_Vector_ui (struct Multiarray_Vector_ui* a)
{
	if (a->owns_data)
		destructor_Vector_ui_2(a->data,compute_size(a->order,a->extents));
	free(a->extents);
	free(a);
}

// Helper functions ************************************************************************************************* //

size_t* set_extents (const size_t order, va_list ap)
{
	size_t* extents = mallocator(SIZE_T_T,1,order); // returned

	for (size_t i = 0; i < order; i++)
		extents[i] = va_arg(ap,size_t);

	return extents;
}

size_t compute_size (const size_t order, const size_t *const extents)
{
	size_t size = 1;
	for (size_t i = 0; i < order; i++)
		size *= extents[i];

	return size;
}

void set_Multiarray_Vector_ui_ui
	(struct Multiarray_Vector_ui* a, const unsigned int*data_V, const unsigned int*const ext_V)
{
	const size_t size = compute_size(a->order,a->extents);
	for (size_t i = 0; i < size; i++) {
		reserve_Vector_ui(a->data[i],ext_V[i]);
		for (size_t j = 0; j < ext_V[i]; j++)
			a->data[i]->data[j] = *data_V++;
	}
}

// Printing functions *********************************************************************************************** //

void print_Multiarray_Vector_ui (const struct Multiarray_Vector_ui*const a)
{
	const size_t order          = a->order;
	const size_t *const extents = a->extents;

	printf("Multi-array extents:");
	for (size_t i = 0; i < order; i++)
		printf(" %zu",extents[i]);
	printf("\n\n");

	switch (order) {
	case 1:
		for (size_t i = 0; i < extents[0]; i++)
			print_Vector_ui(a->data[i]);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void print_const_Multiarray_Vector_ui (const struct const_Multiarray_Vector_ui*const a)
{
	struct Multiarray_Vector_ui local =
		make_local_Multiarray_Vector_ui_0('R',a->order,(size_t*)a->extents,false,(struct Vector_ui**)a->data);
	print_Multiarray_Vector_ui(&local);
}

// Static functions ************************************************************************************************* //

/// \brief Make a local copy of a \ref Multiarray_Vector_ui (static memory).
static struct Multiarray_Vector_ui make_local_Multiarray_Vector_ui_0
	(const char layout,           ///< Standard.
	 const size_t order,          ///< Standard.
	 size_t*const extents,        ///< Standard.
	 const bool owns_data,        ///< Standard.
	 struct Vector_ui**const data ///< Standard.
	)
{
	struct Multiarray_Vector_ui local = { .layout    = layout,
	                                      .order     = order,
	                                      .extents   = extents,
	                                      .owns_data = owns_data,
	                                      .data      = data,
	                                    };
	return local;
}
