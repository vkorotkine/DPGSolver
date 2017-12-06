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

#include "geometry.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "computational_elements.h"
#include "volume_solver.h"
#include "face_solver.h"

#include "const_cast.h"
#include "element_geometry.h"
#include "intrusive.h"
#include "operator.h"
#include "multiarray_operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief See return.
 *  \return The permutation required for conversion to the standard Jacobian ordering from the transposed ordering. */
static const ptrdiff_t* set_jacobian_permutation
	(const int d ///< The dimension
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "geometry_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const ptrdiff_t* set_jacobian_permutation (const int d)
{
	switch (d) {
	case 1:
		return NULL;
		break;
	case 2: {
		static const ptrdiff_t perm_2d[] = { 0, 2, 1, 3, };
		return perm_2d;
		break;
	} case 3: {
		static const ptrdiff_t perm_3d[] = { 0, 3, 6, 1, 4, 7, 2, 5, 8, };
		return perm_3d;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",d);
		break;
	}
	return NULL;
}
