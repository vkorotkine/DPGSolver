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
 *  \brief Undefine macro definitions for c-style templated relating to templates_vector_\*.h.
 */

#include "undef_templates_vector_constructors.h"
#include "undef_templates_vector_math.h"
#include "undef_templates_vector_print.h"

#undef Vector_T
#undef const_Vector_T


#undef cmp_T

#undef reorder_Vector_T
#undef resize_Vector_T
#undef set_to_zero_Vector_T
#undef set_to_data_Vector_T
#undef set_to_value_Vector_T
#undef sort_Vector_T
#undef sum_Vector_T
#undef prod_Vector_T
#undef prod_const_Vector_T
#undef check_equal_Vector_T
#undef check_equal_Vector_T_T
#undef cmp_Vector_T
#undef copy_data_Vector_T_Vector_T
#undef push_back_Vector_T
#undef find_val_Vector_T
