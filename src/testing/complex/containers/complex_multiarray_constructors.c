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

#include <stdlib.h>
#include <assert.h>

#include "macros.h"

#include "complex_multiarray_minimal.h"
#include "multiarray.h"

// Static function declarations ************************************************************************************* //

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

/** \brief Same as \ref constructor_empty_Multiarray_c_dyn_extents but with data calloc'ed.
 *  \return Standard. */
struct Multiarray_c* constructor_zero_Multiarray_c_dyn_extents
	(const char layout,            ///< Defined in \ref Multiarray_d.
	 const int order,              ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_d.
	);

// Interface functions ********************************************************************************************** //

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

struct Multiarray_c* constructor_empty_Multiarray_c_dyn_extents
	(const char layout, const int order, const ptrdiff_t*const extents)
{
	double complex* data = malloc(compute_size(order,extents) * sizeof *data); // keep
	return constructor_move_Multiarray_c_dyn_extents(layout,order,(ptrdiff_t*)extents,true,data);
}

// Zero constructors ************************************************************************************************ //

struct Multiarray_c* constructor_zero_Multiarray_c
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_zero_Multiarray_c_dyn_extents(layout,order,extents);
}

// Copy constructors ************************************************************************************************ //

// Move constructors ************************************************************************************************ //

struct Multiarray_c* constructor_move_Multiarray_c_c
	(const char layout, const int order, const ptrdiff_t*const extents_i, const bool owns_data,
	 double complex*const data)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_move_Multiarray_c_dyn_extents(layout,order,extents,owns_data,data);
}

const struct const_Multiarray_c* constructor_move_const_Multiarray_c_c
	(const char layout, const int order, const ptrdiff_t*const extents_i, const bool owns_data,
	 const double complex*const data)
{
	return (const struct const_Multiarray_c*)
		constructor_move_Multiarray_c_c(layout,order,extents_i,owns_data,(double complex*)data);
}

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

void destructor_const_Multiarray_c (const struct const_Multiarray_c* a)
{
	destructor_Multiarray_c((struct Multiarray_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

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

struct Multiarray_c* constructor_zero_Multiarray_c_dyn_extents
	(const char layout, const int order, const ptrdiff_t*const extents)
{
	double complex* data = calloc(compute_size(order,extents) , sizeof *data); // keep
	return constructor_move_Multiarray_c_dyn_extents(layout,order,(ptrdiff_t*)extents,true,data);
}
