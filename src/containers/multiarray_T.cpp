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

#include <assert.h>
#include <math.h>

#include "macros.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "def_templates_math_functions.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for a vector with an associated index.
struct Vector_T_indexed {
	ptrdiff_t index;         ///< The index.
	struct Vector_T* vector; ///< The \ref Vector_T\*.
};

/** \brief Move constructor for a \ref Vector_T_indexed\*\* from a \ref Vector_T\*\*.
 *  \return Standard. */
static struct Vector_T_indexed** constructor_move_Vector_T_indexed
	(const ptrdiff_t size,  ///< The number of elements.
	 struct Vector_T** data ///< The data to be moved.
	);

/// \brief Destructor for a \ref Vector_T_indexed\*\*.
static void destructor_Vector_T_indexed
	(struct Vector_T_indexed** src, ///< Standard.
	 const ptrdiff_t size           ///< The number of elements.
	);

/** \brief Comparison function for std::qsort between `struct Vector_T_indexed**` `a` and `b`.
 *  \return The lexicographical comparison of `a` and `b`.
 *
 *  \note Input Vectors must be have sorted data.
 */
static int cmp_Vector_T_indexed
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/** \brief Comparison function for std::qsort between `struct Vector_T*` `a` and `b`.
 *  \return The lexicographical comparison of `a` and `b`.
 *
 *  \note The value of CMP_TOL must be set to the desired comparison tolerance before calling this function.
 */
static int cmp_Vector_tol_T
	(const void* a, ///< Variable 1.
	 const void* b  ///< Variable 2.
	 );

/** \brief Compute the total number of entries in a \ref Multiarray_Vector_T\*.
 *  \return See brief. */
static ptrdiff_t compute_total_entries
	(const struct Multiarray_Vector_T*const src ///< Standard.
	);

// Interface functions ********************************************************************************************** //

Type* get_row_Multiarray_T (const ptrdiff_t row, const struct Multiarray_T* a)
{
	assert(a->order == 2);
	assert(a->layout == 'R');
	return &a->data[row*(a->extents[1])];
}

const Type* get_row_const_Multiarray_T (const ptrdiff_t row, const struct const_Multiarray_T* a)
{
	assert(a->order >  0);
	assert(a->order <= 2);
	assert(a->layout == 'R');

	const ptrdiff_t ext_1 = ( a->order == 1 ? 1 : a->extents[1] );
	return &a->data[row*ext_1];
}

Type* get_col_Multiarray_T (const ptrdiff_t col, struct Multiarray_T* a)
{
	assert(a->layout == 'C');

	const ptrdiff_t ext_0 = a->extents[0];
	return &a->data[col*ext_0];
}

const Type* get_col_const_Multiarray_T (const ptrdiff_t col, const struct const_Multiarray_T* a)
{
	return (const Type*) get_col_Multiarray_T(col,(struct Multiarray_T*)a);
}

void remove_col_Multiarray_T (const ptrdiff_t col, struct Multiarray_T* a)
{
	assert(a->layout == 'C');

	Type*const data_col   = get_col_Multiarray_T(col,a),
	    *const data_colp1 = get_col_Multiarray_T(col+1,a);

	const ptrdiff_t size    = compute_size(a->order,a->extents),
	                ind_cp1 = (a->extents[0])*(col+1);

	for (ptrdiff_t i = 0; i < size-ind_cp1; ++i)
		data_col[i] = data_colp1[i];
	a->extents[1] -= 1;
}

void set_to_value_Multiarray_T (struct Multiarray_T*const a, const Type val)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] = val;
}

void set_Multiarray_Vector_T_T
	(struct Multiarray_Vector_T* a, const Type* data_V, const int*const ext_V)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; i++) {
		resize_Vector_T(a->data[i],ext_V[i]);
		for (ptrdiff_t j = 0; j < ext_V[i]; j++)
			a->data[i]->data[j] = *data_V++;
	}
}

void set_Multiarray_T (struct Multiarray_T* a_o, const struct const_Multiarray_T* a_i)
{
	assert(a_o->owns_data);

	resize_Multiarray_T(a_o,a_i->order,a_i->extents);
	a_o->layout = a_i->layout;

	const ptrdiff_t size = compute_size(a_i->order,a_i->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a_o->data[i] = a_i->data[i];
}

void set_Multiarray_T_Multiarray_R (struct Multiarray_T* a, const struct const_Multiarray_R* b)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	assert(size == compute_size(b->order,b->extents));

	for (int i = 0; i < size; ++i)
		a->data[i] = (Type)b->data[i];
}

struct Vector_i* sort_Multiarray_Vector_T (struct Multiarray_Vector_T* a, const bool return_indices)
{
	// sort the individual Vector entries
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		sort_Vector_T(a->data[i]);

	// sort the Vectors
	if (!return_indices) {
		qsort(a->data,(size_t)size,sizeof(a->data[0]),cmp_Vector_T);
		return NULL;
	} else {
		struct Vector_T_indexed** a_indexed = constructor_move_Vector_T_indexed(size,a->data); // destructed

		qsort(a_indexed,(size_t)size,sizeof(a_indexed[0]),cmp_Vector_T_indexed);

		int*const ordering_i = malloc((size_t)size * sizeof *ordering_i); // keep
		for (ptrdiff_t i = 0; i < size; ++i)
			ordering_i[i] = (int)a_indexed[i]->index;

		reorder_Multiarray_Vector_T(a,ordering_i);

		destructor_Vector_T_indexed(a_indexed,size);

		return constructor_move_Vector_i_i(size,true,ordering_i);
	}
	EXIT_ERROR("Should not have made it here.\n");
}

struct Vector_i* sort_Multiarray_Vector_tol_T (struct Multiarray_Vector_T* a, const bool return_indices, const double tol)
{
	CMP_TOL = tol;
	// Note: Does not sort the individual Vector entries

	const ptrdiff_t size = compute_size(a->order,a->extents);
	if (!return_indices) {
		qsort(a->data,(size_t)size,sizeof(a->data[0]),cmp_Vector_tol_T);
		return NULL;
	} else {
		EXIT_ADD_SUPPORT;
	}
}

void reorder_Multiarray_Vector_T (struct Multiarray_Vector_T*const a, const int*const ordering)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	struct Vector_T* b[size];
	for (ptrdiff_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (ptrdiff_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

struct Vector_T* collapse_Multiarray_Vector_T (const struct Multiarray_Vector_T*const src)
{
	const ptrdiff_t n_entries = compute_total_entries(src);

	Type*const data = malloc((size_t)n_entries * sizeof *data); // keep

	ptrdiff_t ind_d = 0;
	const ptrdiff_t size = compute_size(src->order,src->extents);
	for (ptrdiff_t i = 0; i < size; ++i) {
		struct Vector_T*const src_curr = src->data[i];
		const ptrdiff_t size_V = src_curr->ext_0;
		for (ptrdiff_t j = 0; j < size_V; ++j) {
			data[ind_d] = src_curr->data[j];
			++ind_d;
		}
	}

	struct Vector_T*const dest = calloc(1,sizeof *dest); // returned
	dest->ext_0     = n_entries;
	dest->owns_data = true;
	dest->data      = data;

	return dest;
}

void resize_Multiarray_T (struct Multiarray_T* a, const int order, const ptrdiff_t* extents)
{
#if 1
	assert(a->order == order); // Make flexible if needed.
#else
	a->extents = realloc(a->extents,order * sizeof *a->extents);
#endif
	a->order   = order;
	for (int i = 0; i < order; ++i)
		a->extents[i] = extents[i];

	a->data = realloc(a->data,(size_t)compute_size(order,extents) * sizeof *a->data);
}

const struct const_Vector_T* get_const_Multiarray_Vector_T
	(const struct const_Multiarray_Vector_T* src, const ptrdiff_t*const sub_indices)
{
	assert(src != NULL);
	const struct const_Vector_T*const vec =
		src->data[compute_index_sub_container(src->order,0,src->extents,sub_indices)];
	assert(vec != NULL);
	return vec;
}

struct Vector_T interpret_Multiarray_as_Vector_T (struct Multiarray_T*const a_Ma)
{
	struct Vector_T a =
		{ .ext_0     = compute_size(a_Ma->order,a_Ma->extents),
		  .owns_data = false,
		  .data      = a_Ma->data, };
	return a;
}

struct const_Vector_T interpret_const_Multiarray_as_Vector_T (const struct const_Multiarray_T* a_Ma)
{
	assert(a_Ma->order == 1);
	struct const_Vector_T a =
		{ .ext_0     = a_Ma->extents[0],
		  .owns_data = false,
		  .data      = a_Ma->data, };
	return a;
}

struct Matrix_T interpret_Multiarray_as_Matrix_T (const struct Multiarray_T* a_Ma)
{
	assert(a_Ma->order == 2);
	struct Matrix_T a =
		{ .layout    = a_Ma->layout,
		  .ext_0     = a_Ma->extents[0],
		  .ext_1     = a_Ma->extents[1],
		  .owns_data = false,
		  .data      = a_Ma->data, };
	return a;
}

struct const_Matrix_T interpret_const_Multiarray_as_Matrix_T (const struct const_Multiarray_T* a_Ma)
{
	assert(a_Ma->order == 2);
	struct const_Matrix_T a =
		{ .layout    = a_Ma->layout,
		  .ext_0     = a_Ma->extents[0],
		  .ext_1     = a_Ma->extents[1],
		  .owns_data = false,
		  .data      = a_Ma->data, };
	return a;
}

struct Multiarray_T interpret_Multiarray_as_slice_T
	(const struct Multiarray_T* src, const int order_o, const ptrdiff_t*const sub_indices)
{
	const ptrdiff_t ind_data = compute_index_sub_container(src->order,order_o,src->extents,sub_indices);

	struct Multiarray_T dest =
		{ .layout    = src->layout,
		  .order     = order_o,
		  .extents   = src->extents,
		  .owns_data = false,
		  .data      = &src->data[ind_data], };
	return dest;
}

struct const_Multiarray_T interpret_const_Multiarray_as_slice_T
	(const struct const_Multiarray_T*const src, const int order_o, const ptrdiff_t*const sub_indices)
{
	const ptrdiff_t ind_data = compute_index_sub_container(src->order,order_o,src->extents,sub_indices);

	struct const_Multiarray_T dest =
		{ .layout    = src->layout,
		  .order     = order_o,
		  .extents   = src->extents,
		  .owns_data = false,
		  .data      = &src->data[ind_data], };
	return dest;
}

struct Vector_T interpret_Multiarray_slice_as_Vector_T
	(const struct Multiarray_T*const src, const ptrdiff_t*const sub_indices)
{
	struct Multiarray_T tmp = interpret_Multiarray_as_slice_T(src,1,sub_indices);
	return interpret_Multiarray_as_Vector_T(&tmp);
}

void copy_into_Multiarray_T (struct Multiarray_T*const dest, const struct const_Multiarray_T*const src)
{
	assert(dest->order == src->order);
	for (int i = 0; i < src->order; ++i)
		assert(dest->extents[i] == src->extents[i]);
	assert(dest->layout == src->layout);
	assert(dest->owns_data == true);

	const ptrdiff_t size = compute_size(src->order,src->extents);
	for (int i = 0; i < size; ++i)
		dest->data[i] = src->data[i];
}

#ifdef TYPE_RC
void copy_into_Multiarray_T_from_R (struct Multiarray_T*const dest, const struct const_Multiarray_R*const src)
{
	assert(dest->order == src->order);
	for (int i = 0; i < src->order; ++i)
		assert(dest->extents[i] == src->extents[i]);
	assert(dest->layout == src->layout);
	assert(dest->owns_data == true);

	const ptrdiff_t size = compute_size(src->order,src->extents);
	for (int i = 0; i < size; ++i)
		dest->data[i] = src->data[i];
}
#endif

void push_back_Multiarray_T (struct Multiarray_T*const dest, const struct const_Multiarray_T*const src)
{
	assert(dest->order == src->order);
	assert(dest->order == 2); // Can possibly be made flexible if needed.
	assert(dest->layout == src->layout);

	const int order = dest->order;
	const char layout = dest->layout;
	const ptrdiff_t*const extents_d_i = dest->extents,
	               *const extents_s_i = src->extents;
	ptrdiff_t extents_o[order];
	ptrdiff_t rc_start = -1;
	for (int i = 0; i < order; ++i)
		extents_o[i] = extents_d_i[i];

	int ind_mod[2];
	if (layout == 'R') {
		ind_mod[0] = 0;
		ind_mod[1] = 1;
	} else {
		ind_mod[0] = 1;
		ind_mod[1] = 0;
	}
	extents_o[ind_mod[0]] += extents_s_i[ind_mod[0]];
	assert(extents_o[ind_mod[1]] == extents_s_i[ind_mod[1]]);
	rc_start = extents_d_i[ind_mod[0]];

	resize_Multiarray_T(dest,order,extents_o);
	const Type* data_src = src->data;
	Type* data_dest = ( layout == 'R' ? get_row_Multiarray_T(rc_start,dest) : get_col_Multiarray_T(rc_start,dest) );

	const ptrdiff_t size = compute_size(src->order,src->extents);
	for (int i = 0; i < size; ++i)
		*data_dest++ = *data_src++;
}

#ifdef TYPE_RC

const struct const_Vector_i* make_unique_row_Multiarray_T
	(struct Multiarray_T*const src, const Real tol, const bool return_indices)
{
	const ptrdiff_t ext_0 = src->extents[0],
	                ext_1 = compute_size(src->order,src->extents)/ext_0;

	int ind_u[ext_0];
	ind_u[0] = 0;

	int n_u = 1;
	for (int i = 1; i < ext_0; ++i) {
		Type*const data_u = get_row_Multiarray_T(ind_u[n_u-1],src);
		Type*const data_i = get_row_Multiarray_T(i,src);

		bool unique = false;
		for (int j = 0; j < ext_1; ++j) {
			if (!equal_T(data_u[j],data_i[j],tol)) {
				unique = true;
				break;
			}
		}
		if (unique) {
			ind_u[n_u] = i;
			++n_u;
		}
	}

	for (int i = 1; i < n_u; ++i) {
		if (i == ind_u[i])
			continue;

		Type*const data_u = get_row_Multiarray_T(ind_u[i],src);
		Type*const data_i = get_row_Multiarray_T(i,src);
		for (int j = 0; j < ext_1; ++j)
			data_i[j] = data_u[j];
	}

	const int order = src->order;
	ptrdiff_t extents_new[order];
	extents_new[0] = n_u;
	for (int i = 1; i < order; ++i)
		extents_new[i] = src->extents[i];

	resize_Multiarray_T(src,order,extents_new);
	if (!return_indices) {
		return NULL;
	} else {
		int* indices = malloc((size_t)ext_0 * sizeof *indices); // keep
		for (int i = 0, u = 0; i < ext_0; ++i) {
			indices[i] = u;
			if (u < n_u-1 && ind_u[u+1] == i+1)
				++u;
		}
		return constructor_move_const_Vector_i_i(ext_0,true,indices);
	}
}

void update_rows_Multiarray_T
	(struct Multiarray_T*const dest, const struct const_Multiarray_T*const src,
	 const struct const_Vector_i*const row_inds)
{
	assert(dest->order == 2); // Can be made flexible if desired.
	assert(src->order == 2);
	assert(dest->extents[1] == src->extents[1]);
	assert(row_inds->ext_0 == src->extents[0]);

	const bool requires_transpose = ( dest->layout == 'C' ? true : false );
	if (requires_transpose)
		transpose_Multiarray_T(dest,true);

	for (int i = 0; i < row_inds->ext_0; ++i) {
		Type*const dest_row      = get_row_Multiarray_T(row_inds->data[i],dest);
		const Type*const src_row = get_row_const_Multiarray_T(i,src);
		for (int j = 0; j < src->extents[1]; ++j)
			dest_row[j] = src_row[j];
	}

	if (requires_transpose)
		transpose_Multiarray_T(dest,true);
}

#endif

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Vector_T_indexed** constructor_move_Vector_T_indexed (const ptrdiff_t size, struct Vector_T** data)
{
	struct Vector_T_indexed** dest = malloc((size_t)size * sizeof *dest); // returned

	for (ptrdiff_t i = 0; i < size; i++) {
		dest[i] = calloc(1,sizeof *dest[i]); // keep
		dest[i]->index  = i;
		dest[i]->vector = data[i];
	}
	return dest;
}

static void destructor_Vector_T_indexed (struct Vector_T_indexed** src, const ptrdiff_t size)
{
	for (ptrdiff_t i = 0; i < size; ++i)
		free(src[i]);
	free(src);
}

static int cmp_Vector_T_indexed (const void *a, const void *b)
{
	const struct Vector_T_indexed*const*const ia = (const struct Vector_T_indexed*const*const) a;
	const struct Vector_T_indexed*const*const ib = (const struct Vector_T_indexed*const*const) b;

	const ptrdiff_t size_a = (*ia)->vector->ext_0;
	const ptrdiff_t size_b = (*ib)->vector->ext_0;

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const Type*const data_a = (*ia)->vector->data;
	const Type*const data_b = (*ib)->vector->data;

	for (ptrdiff_t i = 0; i < size_a; ++i) {
#if TYPE_RC == TYPE_COMPLEX
		if (creal(data_a[i]) > creal(data_b[i]))
			return 1;
		else if (creal(data_a[i]) < creal(data_b[i]))
			return -1;
#else
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
#endif
	}
	return 0;
}

static int cmp_Vector_tol_T (const void* a, const void* b)
{
	const struct Vector_T*const*const ia = (const struct Vector_T*const*const) a;
	const struct Vector_T*const*const ib = (const struct Vector_T*const*const) b;

	const ptrdiff_t size_a = (*ia)->ext_0;
	const ptrdiff_t size_b = (*ib)->ext_0;

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const Type*const data_a = (*ia)->data;
	const Type*const data_b = (*ib)->data;

	for (ptrdiff_t i = 0; i < size_a; ++i) {
#if TYPE_RC == TYPE_COMPLEX
		const Real da = creal(data_a[i]);
		const Real db = creal(data_b[i]);
		if (abs_R(da-db) < CMP_TOL)
			continue;
#else
		const Type da = data_a[i];
		const Type db = data_b[i];
		if (abs_T(da-db) < CMP_TOL)
			continue;
#endif
		if (da > db)
			return 1;
		else if (da < db)
			return -1;
	}
	return 0;
}

static ptrdiff_t compute_total_entries (const struct Multiarray_Vector_T*const src)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);

	ptrdiff_t count = 0;
	for (ptrdiff_t i = 0; i < size; ++i)
		count += src->data[i]->ext_0;

	return count;
}

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
#include "undef_templates_math_functions.h"
