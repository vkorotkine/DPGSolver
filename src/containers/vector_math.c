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

#include "vector_math.h"

#include <assert.h>

#include "vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void invert_Vector_d (struct Vector_d* a)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		assert(a->data[i] != 0.0);
		a->data[i] = 1.0/(a->data[i]);
	}
}

void add_to_Vector_d_d (struct Vector_d* a, const double* b)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] += b[i];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
