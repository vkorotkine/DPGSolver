// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "multiarray_constructors.h"

#include <stdlib.h>

#include "macros.h"

#include "multiarray.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

/** \brief Allocated and set the `extents` for a `Multiarray_*`.
\todo Change function name.
 *	\return See brief. */
static ptrdiff_t* allocate_and_set_extents
	(const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

struct Multiarray_d* constructor_default_Multiarray_d ()
{
	const int order = 1;
/// \todo Likely remove dynamic memory allocation here.
	ptrdiff_t* extents = calloc(order , sizeof *extents); // keep

	return constructor_move_Multiarray_d_dyn_extents('C',order,extents,true,NULL);
}

// Empty constructors *********************************************************************************************** //

struct Multiarray_d* constructor_empty_Multiarray_d (const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep

	double* data = malloc(compute_size(order,extents) * sizeof *data); // keep

	return constructor_move_Multiarray_d_dyn_extents('C',order,extents,true,data);
}

struct Multiarray_Vector_i* constructor_empty_Multiarray_Vector_i
	(const bool alloc_V, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep

	struct Vector_i** data = NULL;
	if (alloc_V)
		data = constructor_default_Vector_i_2(compute_size(order,extents)); // keep
	else
		data = malloc(compute_size(order,extents) * sizeof *data); // keep

	return constructor_move_Multiarray_Vector_i_dyn_extents(order,extents,true,data);
}

// Copy constructors ************************************************************************************************ //

struct Multiarray_Vector_i* constructor_copy_Multiarray_Vector_i_i
	(const int* data_V, const int*const ext_V, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // free

	struct Multiarray_Vector_i* dest = constructor_empty_Multiarray_Vector_i(true,order,extents); // returned
	free(extents);

	set_Multiarray_Vector_i_i(dest,data_V,ext_V);

	return dest;
}

// Move constructors ************************************************************************************************ //

struct Multiarray_d* constructor_move_Multiarray_d_d
	(const char layout, double*const data, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // free (reallocated in constructor_empty)
	return constructor_move_Multiarray_d_dyn_extents(layout,order,extents,false,data);
}

struct Multiarray_d* constructor_move_Multiarray_d_dyn_extents
	(const char layout, const int order, ptrdiff_t*const extents, const bool owns_data, double*const data)
{
	struct Multiarray_d* dest = calloc(1,sizeof *dest); // returned

	dest->layout    = layout;
	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

struct Multiarray_Vector_i* constructor_move_Multiarray_Vector_i_dyn_extents
	(const int order, ptrdiff_t*const extents, const bool owns_data, struct Vector_i**const data)
{
	struct Multiarray_Vector_i* dest = calloc(1,sizeof *dest); // returned

	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

void const_constructor_move_Multiarray_d (const struct const_Multiarray_d*const* dest, struct Multiarray_d* src)
{
	*(struct const_Multiarray_d**) dest = (struct const_Multiarray_d*) src;
}

void const_constructor_move_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const* dest, struct Multiarray_Vector_i* src)
{
	*(struct const_Multiarray_Vector_i**) dest = (struct const_Multiarray_Vector_i*) src;
}

// Special constructors ********************************************************************************************* //

void set_Multiarray_Matrix_d_from_Multiarray_Matrix_d
	(struct Multiarray_Matrix_d* dest, struct Multiarray_Matrix_d* src, const int order_o,
	 const ptrdiff_t*const sub_indices)
{
	dest->owns_data = false;
	dest->order     = order_o;
	dest->extents   = src->extents;
	dest->data      = &src->data[compute_index_sub_container(src->order,dest->order,src->extents,sub_indices)];
}

void set_const_Multiarray_Matrix_d_from_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d* dest, const struct const_Multiarray_Matrix_d* src, const int order_o,
	 const ptrdiff_t*const sub_indices)
{
	set_Multiarray_Matrix_d_from_Multiarray_Matrix_d(
		(struct Multiarray_Matrix_d*)dest,(struct Multiarray_Matrix_d*)src,order_o,sub_indices);
}


// Destructors ****************************************************************************************************** //

void destructor_Multiarray_d (struct Multiarray_d* a)
{
	if (a == NULL)
		EXIT_DESTRUCTOR;

	free(a->extents);
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_Multiarray_Vector_i (struct Multiarray_Vector_i* a)
{
	if (a == NULL)
		EXIT_DESTRUCTOR;

	destructor_Vector_i_2(a->data,compute_size(a->order,a->extents),a->owns_data);
	free(a->extents);
	free(a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static ptrdiff_t* allocate_and_set_extents (const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t* extents = malloc(order * sizeof *extents); // returned

	for (ptrdiff_t i = 0; i < order; i++)
		extents[i] = extents_i[i];

	return extents;
}
