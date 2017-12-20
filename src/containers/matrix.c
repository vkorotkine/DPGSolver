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

#include "matrix.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_matrix.h"
#include "matrix_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

#include "def_templates_type_i.h"
#include "def_templates_matrix_i.h"
#include "matrix_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void swap_layout (char*const layout)
{
	*layout = ( *layout == 'R' ? 'C' : 'R' );
}

void swap_layout_and_extents (char*const layout, ptrdiff_t*const ext_0, ptrdiff_t*const ext_1)
{
	swap_layout(layout);

	const ptrdiff_t tmp = *ext_0;
	*ext_0 = *ext_1;
	*ext_1 = tmp;
}

char compute_opposite_layout (const char layout_i)
{
	assert((layout_i == 'R') || (layout_i == 'C'));
	return ( layout_i == 'R' ? 'C' : 'R' );
}

ptrdiff_t compute_index_Matrix
	(const ptrdiff_t i, const ptrdiff_t j, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const char layout)
{
	assert((layout == 'R') || (layout == 'C'));
	return ( layout == 'R' ?  i*ext_1+j : j*ext_0+i );
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
