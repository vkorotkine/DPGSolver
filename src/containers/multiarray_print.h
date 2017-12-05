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

#ifndef DPG__multiarray_print_h__INCLUDED
#define DPG__multiarray_print_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* printing functions.
 */

#include <stddef.h>

#include "def_templates_type_d.h"
#include "def_templates_matrix_d.h"
#include "def_templates_multiarray_d.h"
#include "def_templates_vector_d.h"
#include "multiarray_print_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "def_templates_type_i.h"
#include "def_templates_matrix_i.h"
#include "def_templates_multiarray_i.h"
#include "def_templates_vector_i.h"
#include "multiarray_print_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

/// \brief Print the counter for the indices of order > 2 when printing sub-Matrices of the Multiarray.
void print_Multiarray_counter
	(const int order,              ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const counter ///< The counter for the indices of order > 2.
	);

/// \brief Increment the counter for the indicies of order > 2 by 1.
void increment_counter
	(const int order,               ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents, ///< Defined in \ref Multiarray_T.
	 ptrdiff_t*const counter        ///< The counter for the indices of order > 2.
	);

/// \brief Print the extents of the Multiarray.
void print_Multiarray_extents
	(const int order,              ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_T.
	);

#endif // DPG__multiarray_print_h__INCLUDED
