// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "containers_c.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>

#include "Macros.h"

#include "Multiarray.h"

struct Multiarray_c* constructor_empty_Multiarray_c_1 (const char layout, const size_t order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	size_t* extents = set_extents(order,ap); // keep
	va_end(ap);

	struct Multiarray_c* A = malloc(sizeof *A); // returned

	A->layout    = layout;
	A->order     = order;
	A->extents   = extents;
	A->owns_data = true;
	A->data      = calloc(compute_size(order,extents) , sizeof *(A->data)); // keep

	return A;
}

void destructor_Multiarray_c_1 (struct Multiarray_c* A)
{
	free(A->extents);
	if (A->owns_data)
		free(A->data);
	FREE_NULL(A);
}
