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

#ifndef DPG__vector_constructors_h__INCLUDED
#define DPG__vector_constructors_h__INCLUDED
/** \file
 *  \brief Provides Vector_\* constructors and destructors.
 */

#include "def_templates_type_d.h"
#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "vector_constructors_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "def_templates_type_i.h"
#include "def_templates_matrix_i.h"
#include "def_templates_multiarray_i.h"
#include "def_templates_vector_i.h"
#include "vector_constructors_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

/** \brief Constructor the \ref Vector_i of column indices where each of the rows of the input matrix equals 1.0.
 *  \return See brief. */
const struct const_Vector_i* constructor_const_Vector_i_inds_eq_1_const_Matrix_d
	(const struct const_Matrix_d*const a ///< Input matrix.
	);

#endif // DPG__vector_constructors_h__INCLUDED
