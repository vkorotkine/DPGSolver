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
/// \file

#include "element_operators_tp.h"

#include <assert.h>

#include "macros.h"

#include "multiarray.h"
#include "vector.h"

#include "simulation.h"
#include "element_operators.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_operators_tp
	(struct Operators_TP* ops_tp, const struct const_Multiarray_Matrix_d* op_00,
	 const struct const_Multiarray_Matrix_d* op_01, const struct const_Multiarray_Matrix_d* op_10,
	 const struct const_Multiarray_Matrix_d* op_11)
{
	ops_tp->op[0][0] = op_00;
	ops_tp->op[0][1] = op_01;
	ops_tp->op[1][0] = op_10;
	ops_tp->op[1][1] = op_11;
}

const struct const_Multiarray_Matrix_d* constructor_operators_tp
	(const char*const name_type, const char*const name_in, const char*const name_out, const char*const name_range,
	 const struct const_Element* element, const struct Simulation* sim, const struct Operators_TP* ops_tp)
{
	const int p_ref[2] = { -1, -1 };
	struct Operator_Info* op_info =
		constructor_Operator_Info(name_type,name_in,name_out,name_range,p_ref,element); // destructed

printf("ex_op_tp\n");
print_const_Vector_i(op_info->extents_op);
	const struct const_Multiarray_Matrix_d* op =
		constructor_empty_const_Multiarray_Matrix_d_V(false,op_info->extents_op); // returned
UNUSED(sim);
UNUSED(ops_tp);

	destructor_Operator_Info(op_info);

	return op;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

