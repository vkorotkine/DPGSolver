// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "Multiarray.h"

#include <string.h>

#include "Macros.h"

#include "allocators.h"
#include "Vector.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for a vector with an associated index.
struct Vector_ui_indexed {
	size_t index;             ///< The index.
	struct Vector_ui* vector; ///< The \ref Vector_ui\*.
};

/// Move constructs for a \ref Vector_ui_indexed\*\* from a \ref Vector_ui\*\*.
static struct Vector_ui_indexed** constructor_move_Vector_ui_indexed
	(const size_t size,      ///< The number of elements.
	 struct Vector_ui** data ///< The data to be moved.
	);

/// Destructor for a \ref Vector_ui_indexed\*\*.
static void destructor_Vector_ui_indexed
	(struct Vector_ui_indexed** src, ///< Standard.
	 const size_t size               ///< The number of elements.
	);

/** \brief Make a local copy of a \ref Multiarray_Vector_ui (static memory).
 *	\return The \ref Multiarray_Vector_ui.
 */
static struct Multiarray_Vector_ui make_local_Multiarray_Vector_ui_0
	(const char layout,           ///< Standard.
	 const size_t order,          ///< Standard.
	 size_t*const extents,        ///< Standard.
	 const bool owns_data,        ///< Standard.
	 struct Vector_ui**const data ///< Standard.
	);

/** \brief Comparison function for std::qsort between `struct Vector_ui_indexed**` `a` and `b`.
 *	\return The lexicographical comparison of `a` and `b`.
 *
 *	\note Input Vectors must be have sorted data.
 */
static int cmp_Vector_ui_indexed
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/** \brief Reorder a \ref Multiarray_Vector_ui based on the provided ordering.
 *	\warning This is not currently done in place.
 */
static void reorder_Multiarray_Vector_ui
	(struct Multiarray_Vector_ui*const a, ///< Standard.
	 const unsigned int*const ordering    ///< The ordering.
	);

/// \brief Compute the total number of entries in a \ref Multiarray_Vector_ui\*.
static size_t compute_total_entries
	(const struct Multiarray_Vector_ui*const src ///< Standard.
	);

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

struct Multiarray_Vector_ui* constructor_copy_Multiarray_Vector_ui
	(const unsigned int* data_V, const unsigned int*const ext_V, const size_t order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	size_t* extents = set_extents(order,ap); // free (reallocated in constructor_empty)
	va_end(ap);

	size_t extents_l[order];
	for (size_t i = 0; i < order; ++i)
		extents_l[i] = extents[i];
	free(extents);

	struct Multiarray_Vector_ui* dest;
	switch (order) {
	case 1:
		dest = constructor_empty_Multiarray_Vector_ui(order,extents_l[0]);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	set_Multiarray_Vector_ui_ui(dest,data_V,ext_V);

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
	destructor_Vector_ui_2(a->data,compute_size(a->order,a->extents),a->owns_data);
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
	(struct Multiarray_Vector_ui* a, const unsigned int* data_V, const unsigned int*const ext_V)
{
	const size_t size = compute_size(a->order,a->extents);
	for (size_t i = 0; i < size; i++) {
		reserve_Vector_ui(a->data[i],ext_V[i]);
		for (size_t j = 0; j < ext_V[i]; j++)
			a->data[i]->data[j] = *data_V++;
	}
}

struct Vector_ui* sort_Multiarray_Vector_ui (struct Multiarray_Vector_ui* a, const bool return_indices)
{
	// sort the individual Vector entries
	const size_t size = compute_size(a->order,a->extents);
	for (size_t i = 0; i < size; ++i)
		sort_Vector_ui(a->data[i]);

	// sort the Vectors
	unsigned int* ordering_ui = NULL;
	if (!return_indices) {
		qsort(a->data,size,sizeof(a->data[0]),cmp_Vector_ui);
	} else {
		struct Vector_ui_indexed** a_indexed = constructor_move_Vector_ui_indexed(size,a->data); // destructed

		qsort(a_indexed,size,sizeof(a_indexed[0]),cmp_Vector_ui_indexed);

		ordering_ui = malloc(size * sizeof *ordering_ui); // keep
		for (size_t i = 0; i < size; ++i)
			ordering_ui[i] = a_indexed[i]->index;

		reorder_Multiarray_Vector_ui(a,ordering_ui);

		destructor_Vector_ui_indexed(a_indexed,size);
	}

	struct Vector_ui* ordering = constructor_move_Vector_ui_ui(size,true,ordering_ui); // returned

	return ordering;
}

struct Vector_ui* collapse_Multiarray_Vector_ui (const struct Multiarray_Vector_ui*const src)
{
	const size_t n_entries = compute_total_entries(src);

	unsigned int*const data = malloc(n_entries * sizeof *data); // keep

	size_t ind_d = 0;
	const size_t size = compute_size(src->order,src->extents);
	for (size_t i = 0; i < size; ++i) {
		struct Vector_ui*const src_curr = src->data[i];
		const size_t size_V = compute_size(1,src_curr->extents);
		for (size_t j = 0; j < size_V; ++j) {
			data[ind_d] = src_curr->data[j];
			++ind_d;
		}
	}

	struct Vector_ui*const dest = malloc(sizeof *dest); // returned
	dest->extents[0] = n_entries;
	dest->owns_data  = true;
	dest->data       = data;

	return dest;
}

// Printing functions *********************************************************************************************** //

void print_Multiarray_Vector_ui (const struct Multiarray_Vector_ui*const a)
{
	const size_t order          = a->order;
	const size_t *const extents = a->extents;

	printf("Multi-array extents: {");
	for (size_t i = 0; i < order; i++)
		printf(" %zu,",extents[i]);
	printf(" }\n\n");

	switch (order) {
	case 1:
		for (size_t i = 0; i < extents[0]; i++)
			print_Vector_ui(a->data[i]);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	printf("\n");
}

void print_const_Multiarray_Vector_ui (const struct const_Multiarray_Vector_ui*const a)
{
	struct Multiarray_Vector_ui local =
		make_local_Multiarray_Vector_ui_0('R',a->order,(size_t*)a->extents,false,(struct Vector_ui**)a->data);
	print_Multiarray_Vector_ui(&local);
}

// Static functions ************************************************************************************************* //

static struct Vector_ui_indexed** constructor_move_Vector_ui_indexed
	(const size_t size, struct Vector_ui** data)
{
	struct Vector_ui_indexed** dest = malloc(size * sizeof *dest); // returned

	for (size_t i = 0; i < size; i++) {
		dest[i] = malloc(sizeof *dest[i]); // keep
		dest[i]->index  = i;
		dest[i]->vector = data[i];
	}
	return dest;
}

static void destructor_Vector_ui_indexed (struct Vector_ui_indexed** src, const size_t size)
{
	for (size_t i = 0; i < size; ++i)
		free(src[i]);
	free(src);
}

static struct Multiarray_Vector_ui make_local_Multiarray_Vector_ui_0
	(const char layout, const size_t order, size_t*const extents, const bool owns_data, struct Vector_ui**const data)
{
	struct Multiarray_Vector_ui local =
		{ .layout    = layout,
		  .order     = order,
		  .extents   = extents,
		  .owns_data = owns_data,
		  .data      = data,
		};
	return local;
}

static void reorder_Multiarray_Vector_ui (struct Multiarray_Vector_ui*const a, const unsigned int*const ordering)
{
	const size_t size = compute_size(1,a->extents);

	struct Vector_ui* b[size];
	for (size_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (size_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

static int cmp_Vector_ui_indexed (const void *a, const void *b)
{
	const struct Vector_ui_indexed*const*const ia = (const struct Vector_ui_indexed*const*const) a,
	                              *const*const ib = (const struct Vector_ui_indexed*const*const) b;

	const size_t size_a = compute_size(1,(*ia)->vector->extents),
	             size_b = compute_size(1,(*ib)->vector->extents);

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const unsigned int*const data_a = (*ia)->vector->data,
	                  *const data_b = (*ib)->vector->data;

	for (size_t i = 0; i < size_a; ++i) {
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
	}
	return 0;
}

static size_t compute_total_entries (const struct Multiarray_Vector_ui*const src)
{
	const size_t size = compute_size(src->order,src->extents);

	size_t count = 0;
	for (size_t i = 0; i < size; ++i)
		count += compute_size(1,src->data[i]->extents);

	return count;
}

