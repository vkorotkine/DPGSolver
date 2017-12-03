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

#include "complex_multiarray.h"

#include <assert.h>

#include "macros.h"

#include "multiarray.h"
#include "complex_matrix.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "def_templates_multiarray_c.h"
#include "multiarray_T.c"

void set_Multiarray_c_Multiarray_d (struct Multiarray_c* a, const struct const_Multiarray_d* b)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	assert(size == compute_size(b->order,b->extents));

	for (int i = 0; i < size; ++i)
		a->data[i] = b->data[i];
}

struct Matrix_c interpret_Multiarray_as_Matrix_c (const struct Multiarray_c* a_Ma)
{
	assert(a_Ma->order == 2);
	struct Matrix_c a =
		{ .layout    = a_Ma->layout,
		  .ext_0     = a_Ma->extents[0],
		  .ext_1     = a_Ma->extents[1],
		  .owns_data = false,
		  .data      = a_Ma->data, };
	return a;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
