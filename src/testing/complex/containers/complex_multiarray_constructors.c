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
static struct Multiarray_c* constructor_zero_Multiarray_c_dyn_extents
	(const char layout,            ///< Defined in \ref Multiarray_d.
	 const int order,              ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_d.
	);

/** \brief Copy constructor for a \ref Multiarray_c\* from a \ref Multiarray_d\*.
 *  \return Standard. */
static struct Multiarray_c* constructor_copy_Multiarray_c_Multiarray_d
	(struct Multiarray_d* src ///< Source.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "def_templates_multiarray_c.h"
#include "def_templates_multiarray_constructors_c.h"
#include "multiarray_constructors_T.c"

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

// Zero constructors ************************************************************************************************ //

struct Multiarray_c* constructor_zero_Multiarray_c
	(const char layout, const int order, const ptrdiff_t*const extents_i)
{
	ptrdiff_t*const extents = allocate_and_set_extents(order,extents_i); // keep
	return constructor_zero_Multiarray_c_dyn_extents(layout,order,extents);
}

// Copy constructors ************************************************************************************************ //

const struct const_Multiarray_c* constructor_copy_const_Multiarray_c_Multiarray_d
	(const struct const_Multiarray_d* src)
{
	return (struct const_Multiarray_c*) constructor_copy_Multiarray_c_Multiarray_d((struct Multiarray_d*)src);
}

// Move constructors ************************************************************************************************ //

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Multiarray_c* constructor_zero_Multiarray_c_dyn_extents
	(const char layout, const int order, const ptrdiff_t*const extents)
{
	double complex* data = calloc((size_t)compute_size(order,extents) , sizeof *data); // keep
	return constructor_move_Multiarray_c_dyn_extents(layout,order,(ptrdiff_t*)extents,true,data);
}

static struct Multiarray_c* constructor_copy_Multiarray_c_Multiarray_d (struct Multiarray_d* src)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);
	double complex* data = malloc((size_t)size * sizeof *data); // moved
	for (int i = 0; i < size; ++i)
		data[i] = src->data[i];

	return constructor_move_Multiarray_c_c(src->layout,src->order,src->extents,true,data);
}
