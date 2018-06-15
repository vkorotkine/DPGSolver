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
 *  \brief Undefine macro definitions for c-style templating relating to def_templates_multiarray_\*.h.
 */

#include "undef_templates_multiarray_constructors.h"
#include "undef_templates_multiarray_math.h"
#include "undef_templates_multiarray_print.h"


#undef Multiarray_T
#undef const_Multiarray_T

#undef Multiarray_Vector_T
#undef const_Multiarray_Vector_T
#undef Multiarray_Matrix_T
#undef const_Multiarray_Matrix_T


#undef Vector_T_indexed

#undef constructor_move_Vector_T_indexed
#undef destructor_Vector_T_indexed
#undef cmp_Vector_T_indexed
#undef reorder_Multiarray_Vector_T
#undef compute_total_entries

#undef get_row_Multiarray_T
#undef get_row_const_Multiarray_T
#undef get_col_Multiarray_T
#undef get_col_const_Multiarray_T
#undef remove_col_Multiarray_T
#undef set_to_value_Multiarray_T
#undef set_Multiarray_Vector_T_T
#undef set_Multiarray_T
#undef set_Multiarray_T_Multiarray_R

#undef sort_Multiarray_Vector_T
#undef collapse_Multiarray_Vector_T
#undef resize_Multiarray_T
#undef get_const_Multiarray_Vector_T
#undef interpret_Multiarray_as_Vector_T
#undef interpret_const_Multiarray_as_Vector_T
#undef interpret_Multiarray_as_Matrix_T
#undef interpret_const_Multiarray_as_Matrix_T
#undef interpret_Multiarray_as_slice_T
#undef interpret_const_Multiarray_as_slice_T
#undef interpret_Multiarray_slice_as_Vector_T
#undef copy_into_Multiarray_T
#undef copy_into_Multiarray_T_from_R
#undef update_rows_Multiarray_T
#undef push_back_Multiarray_T
#undef make_unique_row_Multiarray_T


#undef Multiarray_R
#undef const_Multiarray_R

#undef get_row_Multiarray_R
#undef get_row_const_Multiarray_R
#undef get_col_Multiarray_R
#undef get_col_const_Multiarray_R
#undef set_to_value_Multiarray_R
#undef set_Multiarray_R

#undef resize_Multiarray_R
#undef interpret_const_Multiarray_as_Vector_R
#undef interpret_Multiarray_as_Matrix_R
#undef interpret_const_Multiarray_as_Matrix_R
#undef update_rows_Multiarray_R
