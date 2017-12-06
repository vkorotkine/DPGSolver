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

#ifndef DPG__multiarray_math_h__INCLUDED
#define DPG__multiarray_math_h__INCLUDED
/** \file
 *  \brief Provides real Multiarray_\* math functions.
 */

#include <stddef.h>

#include "def_templates_type_d.h"
#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "multiarray_math_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

/** \brief Compute the extents of the output multiarray from a matrix-multiarray multiplication.
 *  \return Dynamically allocated extents. */
ptrdiff_t* compute_extents_mm_MMa
	(const ptrdiff_t ext_0,     ///< The value of `extents[0]`.
	 const int order,           ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t* extents_i ///< The input extents. Used to set all but the first entry.
	);

#endif // DPG__multiarray_math_h__INCLUDED
