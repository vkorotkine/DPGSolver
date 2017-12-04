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
 *  \brief Provides the macro definitions used for c-style templating related to the `double complex` multiarray
 *         containers/functions.
 */

#include "def_templates_multiarray_constructors_c.h"
#include "def_templates_multiarray_math_c.h"
#include "def_templates_multiarray_print_c.h"

///\{ \name Data types
#define Multiarray_T       Multiarray_c
#define Multiarray_R       Multiarray_d
#define const_Multiarray_T const_Multiarray_c
#define const_Multiarray_R const_Multiarray_d

#define Multiarray_Vector_T       Multiarray_Vector_c
#define const_Multiarray_Vector_T const_Multiarray_Vector_c
#define Multiarray_Matrix_T       Multiarray_Matrix_c
#define const_Multiarray_Matrix_T const_Multiarray_Matrix_c
///\}

///\{ \name Function names
#define Vector_T_indexed Vector_c_indexed

#define constructor_move_Vector_T_indexed constructor_move_Vector_c_indexed
#define destructor_Vector_T_indexed       destructor_Vector_c_indexed
#define cmp_Vector_T_indexed              cmp_Vector_c_indexed
#define reorder_Multiarray_Vector_T       reorder_Multiarray_Vector_c
#define compute_total_entries             compute_total_entries_c

#define get_row_Multiarray_T          get_row_Multiarray_c
#define get_row_const_Multiarray_T    get_row_const_Multiarray_c
#define get_col_Multiarray_T          get_col_Multiarray_c
#define get_col_const_Multiarray_T    get_col_const_Multiarray_c
#define set_to_value_Multiarray_T     set_to_value_Multiarray_c
#define set_Multiarray_Vector_T_T     set_Multiarray_Vector_c_c
#define set_Multiarray_T              set_Multiarray_c
#define set_Multiarray_T_Multiarray_R set_Multiarray_c_Multiarray_d

#define sort_Multiarray_Vector_T               sort_Multiarray_Vector_c
#define collapse_Multiarray_Vector_T           collapse_Multiarray_Vector_c
#define resize_Multiarray_T                    resize_Multiarray_c
#define get_const_Multiarray_Vector_T          get_const_Multiarray_Vector_c
#define interpret_const_Multiarray_as_Vector_T interpret_const_Multiarray_as_Vector_c
#define interpret_const_Multiarray_as_Vector_R interpret_const_Multiarray_as_Vector_d
#define interpret_Multiarray_as_Matrix_T       interpret_Multiarray_as_Matrix_c
///\}
