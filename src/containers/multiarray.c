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
/** \file
 */

#include "multiarray.h"

#include "macros.h"

#include "vector.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for a vector with an associated index.
struct Vector_i_indexed {
	ptrdiff_t index;         ///< The index.
	struct Vector_i* vector; ///< The \ref Vector_i\*.
};

/** \brief Move constructor for a \ref Vector_i_indexed\*\* from a \ref Vector_i\*\*.
 *  \return Standard. */
static struct Vector_i_indexed** constructor_move_Vector_i_indexed
	(const ptrdiff_t size,  ///< The number of elements.
	 struct Vector_i** data ///< The data to be moved.
	);

/// \brief Destructor for a \ref Vector_i_indexed\*\*.
static void destructor_Vector_i_indexed
	(struct Vector_i_indexed** src, ///< Standard.
	 const ptrdiff_t size           ///< The number of elements.
	);

/** \brief Comparison function for std::qsort between `struct Vector_i_indexed**` `a` and `b`.
 *  \return The lexicographical comparison of `a` and `b`.
 *
 *  \note Input Vectors must be have sorted data.
 */
static int cmp_Vector_i_indexed
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/** \brief Reorder a \ref Multiarray_Vector_i based on the provided ordering.
 *  \warning This is not currently done in place.
 */
static void reorder_Multiarray_Vector_i
	(struct Multiarray_Vector_i*const a, ///< Standard.
	 const int*const ordering    ///< The ordering.
	);

/** \brief Compute the total number of entries in a \ref Multiarray_Vector_i\*.
 *  \return See brief. */
static ptrdiff_t compute_total_entries
	(const struct Multiarray_Vector_i*const src ///< Standard.
	);

// Interface functions ********************************************************************************************** //

ptrdiff_t compute_size (const int order, const ptrdiff_t*const extents)
{
	ptrdiff_t size = 1;
	for (ptrdiff_t i = 0; i < order; i++)
		size *= extents[i];

	return size;
}

void set_Multiarray_Vector_i_i
	(struct Multiarray_Vector_i* a, const int* data_V, const int*const ext_V)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; i++) {
		resize_Vector_i(a->data[i],ext_V[i]);
		for (ptrdiff_t j = 0; j < ext_V[i]; j++)
			a->data[i]->data[j] = *data_V++;
	}
}

struct Vector_i* sort_Multiarray_Vector_i (struct Multiarray_Vector_i* a, const bool return_indices)
{
	// sort the individual Vector entries
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		sort_Vector_i(a->data[i]);

	// sort the Vectors
	int* ordering_i = NULL;
	if (!return_indices) {
		qsort(a->data,size,sizeof(a->data[0]),cmp_Vector_i);
	} else {
		struct Vector_i_indexed** a_indexed = constructor_move_Vector_i_indexed(size,a->data); // destructed

		qsort(a_indexed,size,sizeof(a_indexed[0]),cmp_Vector_i_indexed);

		ordering_i = malloc(size * sizeof *ordering_i); // keep
		for (ptrdiff_t i = 0; i < size; ++i)
			ordering_i[i] = a_indexed[i]->index;

		reorder_Multiarray_Vector_i(a,ordering_i);

		destructor_Vector_i_indexed(a_indexed,size);
	}

	struct Vector_i* ordering = constructor_move_Vector_i_i(size,true,ordering_i); // returned

	return ordering;
}

struct Vector_i* collapse_Multiarray_Vector_i (const struct Multiarray_Vector_i*const src)
{
	const ptrdiff_t n_entries = compute_total_entries(src);

	int*const data = malloc(n_entries * sizeof *data); // keep

	ptrdiff_t ind_d = 0;
	const ptrdiff_t size = compute_size(src->order,src->extents);
	for (ptrdiff_t i = 0; i < size; ++i) {
		struct Vector_i*const src_curr = src->data[i];
		const ptrdiff_t size_V = src_curr->ext_0;
		for (ptrdiff_t j = 0; j < size_V; ++j) {
			data[ind_d] = src_curr->data[j];
			++ind_d;
		}
	}

	struct Vector_i*const dest = calloc(1,sizeof *dest); // returned
	dest->ext_0     = n_entries;
	dest->owns_data = true;
	dest->data      = data;

	return dest;
}

ptrdiff_t compute_index_sub_vector (const int order, const ptrdiff_t*const extents, const ptrdiff_t*const sub_indices)
{
	return compute_index_sub_container(order,1,extents,sub_indices);
}

ptrdiff_t compute_index_sub_matrix (const int order, const ptrdiff_t*const extents, const ptrdiff_t*const sub_indices)
{
	return compute_index_sub_container(order,2,extents,sub_indices);
}

ptrdiff_t compute_index_sub_container
	(const int order_i, const int order_o, const ptrdiff_t*const extents, const ptrdiff_t*const sub_indices)
{
	const ptrdiff_t*const extents_tail = &extents[order_o];

	ptrdiff_t base = 1;
	for (int i = 0; i < order_o; ++i)
		base *= extents[i];
	ptrdiff_t ind_sub = 0;
	for (int i = 0; i < order_i-order_o; ++i) {
		ind_sub += base*sub_indices[i];
		base *= extents_tail[i];
	}
	return ind_sub;
}

ptrdiff_t compute_index_sub_container_pi
	(const int order_i, const int order_o, const ptrdiff_t*const extents, const int*const sub_indices)
{
	const int n_sub = order_i-order_o;
	ptrdiff_t sub_indices_p[n_sub];
	for (int i = 0; i < n_sub; ++i)
		sub_indices_p[i] = sub_indices[i];

	return compute_index_sub_container(order_i,order_o,extents,sub_indices_p);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Vector_i_indexed** constructor_move_Vector_i_indexed
	(const ptrdiff_t size, struct Vector_i** data)
{
	struct Vector_i_indexed** dest = malloc(size * sizeof *dest); // returned

	for (ptrdiff_t i = 0; i < size; i++) {
		dest[i] = calloc(1,sizeof *dest[i]); // keep
		dest[i]->index  = i;
		dest[i]->vector = data[i];
	}
	return dest;
}

static void destructor_Vector_i_indexed (struct Vector_i_indexed** src, const ptrdiff_t size)
{
	for (ptrdiff_t i = 0; i < size; ++i)
		free(src[i]);
	free(src);
}

static void reorder_Multiarray_Vector_i (struct Multiarray_Vector_i*const a, const int*const ordering)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	struct Vector_i* b[size];
	for (ptrdiff_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (ptrdiff_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

static int cmp_Vector_i_indexed (const void *a, const void *b)
{
	const struct Vector_i_indexed*const*const ia = (const struct Vector_i_indexed*const*const) a,
	                             *const*const ib = (const struct Vector_i_indexed*const*const) b;

	const ptrdiff_t size_a = (*ia)->vector->ext_0,
	                size_b = (*ib)->vector->ext_0;

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const int*const data_a = (*ia)->vector->data,
	         *const data_b = (*ib)->vector->data;

	for (ptrdiff_t i = 0; i < size_a; ++i) {
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
	}
	return 0;
}

static ptrdiff_t compute_total_entries (const struct Multiarray_Vector_i*const src)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);

	ptrdiff_t count = 0;
	for (ptrdiff_t i = 0; i < size; ++i)
		count += src->data[i]->ext_0;

	return count;
}
