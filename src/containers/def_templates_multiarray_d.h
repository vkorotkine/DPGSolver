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
 *  \brief Provides the macro definitions used for c-style templating related to the `double` multiarray
 *         containers/functions.
 */

#include "def_templates_multiarray_constructors_d.h"
#include "def_templates_multiarray_math_d.h"

///\{ \name Data types
#define Multiarray_T       Multiarray_d
#define Multiarray_R       Multiarray_d       
#define const_Multiarray_T const_Multiarray_d
#define const_Multiarray_R const_Multiarray_d
///\}

///\{ \name Function names
#define set_to_value_Multiarray_T  set_to_value_Multiarray_d

#define get_col_const_Multiarray_T get_col_const_Multiarray_d
#define get_col_Multiarray_T       get_col_Multiarray_d

#define interpret_Multiarray_as_Matrix_T interpret_Multiarray_as_Matrix_d
#define interpret_const_Multiarray_as_Vector_R interpret_const_Multiarray_as_Vector_d
///\}
