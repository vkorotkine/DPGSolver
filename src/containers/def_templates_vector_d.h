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
 *  \brief Provides the macro definitions used for c-style templating related to the `double` vector
 *         containers/functions.
 */

#include "def_templates_vector_constructors_d.h"
#include "def_templates_vector_math_d.h"
#include "def_templates_vector_print_d.h"

///\{ \name Data types
#define Vector_T       Vector_d
#define const_Vector_T const_Vector_d
#define Vector_R       Vector_d
#define const_Vector_R const_Vector_d
///\}

///\{ \name Function names
#define cmp_T cmp_d

#define reorder_Vector_T            reorder_Vector_d
#define resize_Vector_T             resize_Vector_d
#define set_to_zero_Vector_T        set_to_zero_Vector_d
#define set_to_data_Vector_T        set_to_data_Vector_d
#define set_to_value_Vector_T       set_to_value_Vector_d
#define sort_Vector_T               sort_Vector_d
#define sum_Vector_T                sum_Vector_d
#define prod_Vector_T               prod_Vector_d
#define prod_const_Vector_T         prod_const_Vector_d
#define check_equal_Vector_T        check_equal_Vector_d
#define check_equal_Vector_T_T      check_equal_Vector_d_d
#define cmp_Vector_T                cmp_Vector_d
#define copy_data_Vector_T_Vector_T copy_data_Vector_d_Vector_d
#define push_back_Vector_T          push_back_Vector_d
#define find_val_Vector_T           find_val_Vector_d
///\}
