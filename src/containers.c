// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "containers.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

/*	Purpose:
 *		Provide functions relating to standard containers.
 *
 *	Notation: (ToBeModified)
 *		Function naming convention:
 *
 *		constructor_[0]_(1)_[2]_[3]_(4):
 *			Elements in square brackets [] are required, those in round brackets () are optional.
 *			[0] : Type of constructor
 *				Options: move
 *			(1) : Optional const specifier
 *			[2] : Type of container to be returned
 *				Options: Multiarray_d
 *			[3] : Level of dereferencing of the returned container object
 *			(4) : Type of input from which the container is constructed
 *
 *	Comments:
 *		The constructors for the const_Multiarray_* structs must initialize const members on the heap. This is done so
 *		that it is not possible to modify values of const moved objects accidentally.
 *
 *		The procedure used here for this purpose was taken from the [SO answer] referenced below. See the discussion
 *		following the accepted answer for the explanation of why this does not violate the standard regarding modifying
 *		an object defined with a const-qualified type.
 *
 *		Several functions from <stdarg.h> are used such that a variable number of arguments can be passed. This is
 *		similar to the implementation of printf. An example with a more detailed explanation can be found on the
 *		referenced Wikipedia article on [stdarg.h].
 *
 *	References:
 *		[SO answer]: https://stackoverflow.com/questions/2219001/how-to-initialize-const-members-of-structs-on-the-heap
 *		[stdarg.h]: https://en.wikipedia.org/wiki/Stdarg.h
 */

/*
static struct Multiarray_d* constructor_default_Multiarray_d_1 (const size_t order)
{
	struct Multiarray_d* A = calloc(1 , sizeof *A); // returned

	A->extents = calloc(order , sizeof *(A->extents)); // keep
	A->order   = order;

	return A;
}
*/

struct Multiarray_d* constructor_move_Multiarray_d_1_d (const char layout, double*const data, const size_t order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	size_t* extents = set_extents(order,ap); // keep
	va_end(ap);

	struct Multiarray_d local = { .layout  = layout,
	                              .order   = order,
	                              .extents = extents,
	                              .data    = data,
	                            };

	struct Multiarray_d* A = malloc(sizeof *A); // returned
	memcpy(A,&local,sizeof *A);

	return A;
}

void destructor_moved_Multiarray_d_1 (struct Multiarray_d* A)
{
	free(A->extents);
	free(A);
}

size_t* set_extents (const size_t order, va_list ap)
{
	size_t* extents = malloc(order * sizeof *extents); // returned

	for (size_t i = 0; i < order; i++)
		extents[i] = va_arg(ap,size_t);

	return extents;
}

size_t compute_size (const size_t n_ext, const size_t *const extents)
{
	size_t size = 1;
	for (size_t i = 0; i < n_ext; i++)
		size *= extents[i];

	return size;
}
