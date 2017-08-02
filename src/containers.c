// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "containers.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "Macros.h"

struct Multiarray_d* constructor_move_Multiarray_d_1_d (const char layout, double*const data, const size_t order, ...)
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

void destructor_Multiarray_d_1 (struct Multiarray_d* A)
{
	free(A->extents);
	if (A->owns_data)
		free(A->data);
	FREE_NULL(A);
}

size_t* set_extents (const size_t order, va_list ap)
{
	size_t* extents = malloc(order * sizeof *extents); // returned

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
