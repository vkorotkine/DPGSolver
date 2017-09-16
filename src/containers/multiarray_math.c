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

#include "multiarray_math.h"

#include "macros.h"

#include "multiarray.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void transpose_Multiarray_d (struct Multiarray_d* a, const bool mem_only)
{
	if (a->order != 2)
		EXIT_UNSUPPORTED;

	struct Matrix_d* dest = constructor_default_Matrix_d(); // destructed
	set_Matrix_from_Multiarray_d(dest,a,(ptrdiff_t[]){});

	transpose_Matrix_d(dest,mem_only);
	destructor_Matrix_d(dest);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
