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
 *  \brief Provides the macro definitions used for c-style templating related to the `int` vector
 *         containers/functions.
 */

#include "def_templates_vector_constructors_i.h"
#include "def_templates_vector_print_i.h"

///\{ \name Data types
#define Vector_T       Vector_i
#define const_Vector_T const_Vector_i
///\}

///\{ \name Function names
#define cmp_T cmp_i

#define reorder_Vector_T            reorder_Vector_i
#define resize_Vector_T             resize_Vector_i
#define set_to_zero_Vector_T        set_to_zero_Vector_i
#define set_to_data_Vector_T        set_to_data_Vector_i
#define set_to_value_Vector_T       set_to_value_Vector_i
#define sort_Vector_T               sort_Vector_i
#define sum_Vector_T                sum_Vector_i
#define prod_Vector_T               prod_Vector_i
#define prod_const_Vector_T         prod_const_Vector_i
#define check_equal_Vector_T        check_equal_Vector_i
#define check_equal_Vector_T_T      check_equal_Vector_i_i
#define cmp_Vector_T                cmp_Vector_i
#define copy_data_Vector_T_Vector_T copy_data_Vector_i_Vector_i
#define push_back_Vector_T          push_back_Vector_i
#define find_val_Vector_T           find_val_Vector_i
///\}
