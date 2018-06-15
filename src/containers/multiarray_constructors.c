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

#include "multiarray_constructors.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "const_cast.h"
#include "operator.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "multiarray_constructors_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "multiarray_constructors_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_i.h"
#include "multiarray_constructors_T.c"
#include "undef_templates_type.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

ptrdiff_t* allocate_and_set_extents (const int order, const ptrdiff_t*const extents_i)
{
	// MSB The extents_i structure is an array of ptrdiff_t that holds the size of each 
	// dimension of the multiarray. Here, I think the issue is that we need to move 
	// a const ptrdiff_t values into a non const array.
	ptrdiff_t* extents = malloc((size_t)order * sizeof *extents); // returned

	for (ptrdiff_t i = 0; i < order; i++)
		extents[i] = extents_i[i];

	return extents;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
