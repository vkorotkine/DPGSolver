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

#include "complex_multiarray_constructors.h"

#include <assert.h>
#include <stdlib.h>

#include "complex_multiarray.h"
#include "multiarray.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for an empty \ref Multiarray_c\* with extents having been previously dynamically allocated.
 *  \return Standard. */
static struct Multiarray_c* constructor_empty_Multiarray_c_dyn_extents
	(const char layout,            ///< Defined in \ref Multiarray_c.
	 const int order,              ///< Defined in \ref Multiarray_c.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_c.
	);

// Interface functions ********************************************************************************************** //

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

struct Multiarray_c* constructor_empty_Multiarray_c
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_empty_Multiarray_c_dyn_extents(layout,order,extents);
}

// Zero constructors ************************************************************************************************ //

// Copy constructors ************************************************************************************************ //

// Move constructors ************************************************************************************************ //

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

void destructor_Multiarray_c (struct Multiarray_c* a)
{
	assert(a != NULL);

	free(a->extents);
	if (a->owns_data)
		free(a->data);
	free(a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Move constructor for a \ref Multiarray_c\* with the input extents having been previously dynamically
 *         allocated.
 *  \return See brief. */
static struct Multiarray_c* constructor_move_Multiarray_c_dyn_extents
	(const char layout,        ///< Standard.
	 const int order,          ///< Standard.
	 ptrdiff_t*const extents,  ///< Standard.
	 const bool owns_data,     ///< Standard.
	 double complex*const data ///< Standard.
	);

static struct Multiarray_c* constructor_empty_Multiarray_c_dyn_extents
	(const char layout, const int order, const ptrdiff_t*const extents)
{
	double complex* data = malloc(compute_size(order,extents) * sizeof *data); // keep
	return constructor_move_Multiarray_c_dyn_extents(layout,order,(ptrdiff_t*)extents,true,data);
}

// Level 1 ********************************************************************************************************** //

static struct Multiarray_c* constructor_move_Multiarray_c_dyn_extents
	(const char layout, const int order, ptrdiff_t*const extents, const bool owns_data, double complex*const data)
{
	struct Multiarray_c* dest = calloc(1,sizeof *dest); // returned

	dest->layout    = layout;
	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}
