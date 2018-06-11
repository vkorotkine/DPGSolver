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

#include "multiarray_print.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "multiarray_print_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "multiarray_print_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_i.h"
#include "multiarray_print_T.c"
#include "undef_templates_type.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void print_Multiarray_counter (const int order, const ptrdiff_t*const counter)
{
	printf("{:,:");
	for (int i = 0; i < order-2; ++i)
		printf(",%td",counter[i]);
	printf("}\n");
}

void increment_counter (const int order, const ptrdiff_t*const extents, ptrdiff_t*const counter)
{
	const ptrdiff_t*const extents_tail = &extents[2];

	for (int i = 0; i < order-2; ++i) {
		++counter[i];
		if (counter[i] == extents_tail[i]) {
			counter[i] = 0;
			continue;
		}
		return;
	}
}

void print_Multiarray_extents (const int order, const ptrdiff_t*const extents)
{
	printf("Multi-array extents: {");
	for (ptrdiff_t i = 0; i < order; i++)
		printf(" %td,",extents[i]);
	printf(" }\n\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
