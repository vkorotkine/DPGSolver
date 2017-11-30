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

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

/** \brief Empty constructor for a Multiarray_\*.
 *  \return See brief. */
struct Multiarray_T* constructor_empty_Multiarray_T
	(const char layout,              ///< Standard.
	 const int order,                ///< Standard.
	 const ptrdiff_t*const extents_i ///< The input extents.
	)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_empty_Multiarray_T_dyn_extents(layout,order,extents);
}

/** \brief Empty constructor for a Multiarray_\* with extents having been previously dynamically allocated.
 *  \return See brief. */
struct Multiarray_T* constructor_empty_Multiarray_T_dyn_extents
	(const char layout,            ///< Standard.
	 const int order,              ///< Standard.
	 const ptrdiff_t*const extents ///< The input extents.
	)
{
	Type* data = malloc((size_t)compute_size(order,extents) * sizeof *data); // keep

	return constructor_move_Multiarray_T_dyn_extents(layout,order,(ptrdiff_t*)extents,true,data);
}

// Zero constructors ************************************************************************************************ //

// Copy constructors ************************************************************************************************ //

/** \brief Copy constructor for a Multiarray_\* from a Multiarray_\*.
 *  \return See brief. */
struct Multiarray_T* constructor_copy_Multiarray_T
	(struct Multiarray_T* src ///< Source.
	)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);
	Type* data = malloc((size_t)size * sizeof *data); // moved
	for (int i = 0; i < size; ++i)
		data[i] = src->data[i];

	return constructor_move_Multiarray_T_T(src->layout,src->order,src->extents,true,data);
}

/** \brief `const` version of \ref constructor_copy_Multiarray_T.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_copy_const_Multiarray_T
	(const struct const_Multiarray_T*const src ///< Source.
	)
{
	return (const struct const_Multiarray_T*) constructor_copy_Multiarray_T((struct Multiarray_T*)src);
}

// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a Multiarray_\* from a `Type*`.
 *  \return See brief.. */
struct Multiarray_T* constructor_move_Multiarray_T_T
	(const char layout,               ///< Standard.
	 const int order,                 ///< Standard.
	 const ptrdiff_t*const extents_i, ///< The input extents.
	 const bool owns_data,            ///< Standard.
	 Type*const data                  ///< Standard.
	)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_move_Multiarray_T_dyn_extents(layout,order,extents,owns_data,data);
}

/** \brief `const` version of \ref constructor_move_Multiarray_T_T.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_move_const_Multiarray_T_T
	(const char layout,               ///< See brief.
	 const int order,                 ///< See brief.
	 const ptrdiff_t*const extents_i, ///< See brief.
	 const bool owns_data,            ///< See brief.
	 const Type*const data            ///< See brief.
	)
{
	return (const struct const_Multiarray_T*)
		constructor_move_Multiarray_T_T(layout,order,extents_i,owns_data,(Type*)data);
}

/** \brief Move constructor for a Multiarray_\* with extents having been previously dynamically allocated.
 *  \return See brief. */
struct Multiarray_T* constructor_move_Multiarray_T_dyn_extents
	(const char layout,       ///< Standard.
	 const int order,         ///< Standard.
	 ptrdiff_t*const extents, ///< Standard.
	 const bool owns_data,    ///< Standard.
	 Type*const data          ///< Standard.
	)
{
	struct Multiarray_T* dest = calloc(1,sizeof *dest); // returned

	dest->layout    = layout;
	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

/// \brief Destructs a Multiarray_\*.
void destructor_Multiarray_T
	(struct Multiarray_T* a ///< Standard.
	)
{
	assert(a != NULL);

	free(a->extents);
	if (a->owns_data)
		free(a->data);
	free(a);
}

/// \brief `const` version of \ref destructor_Multiarray_T.
void destructor_const_Multiarray_T
	(const struct const_Multiarray_T* a ///< Standard.
	)
{
	destructor_Multiarray_T((struct Multiarray_T*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
