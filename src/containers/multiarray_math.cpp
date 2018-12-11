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

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "math_functions.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "multiarray_math_T.cpp"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "multiarray_math_T.cpp"
#include "undef_templates_type.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

ptrdiff_t* compute_extents_mm_MMa (const ptrdiff_t ext_0, const int order, const ptrdiff_t* extents_i)
{
	ptrdiff_t* extents = malloc((size_t)order * sizeof *extents); // returned

	extents[0] = ext_0;
	for (int i = 1; i < order; ++i)
		extents[i] = extents_i[i];

	return extents;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
