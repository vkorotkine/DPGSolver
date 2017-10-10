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
 */

#include "test_unit_operators_tp.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "test_base.h"
#include "test_support.h"
#include "test_support_multiarray.h"
#include "test_support_matrix.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_tol.h"

#include "multiarray.h"
#include "matrix.h"

#include "element_operators_tp.h"
#include "operator.h"

// Static function declarations ************************************************************************************* //

/** \brief Provides unit tests for the construction of the standard operator from the tensor-product of the
 *         sub-operators. */
static void test_unit_construct_std_from_tp
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const e_type           ///< \ref Element::type.
	);

/// \brief Provides unit tests for the application of tensor-product sub-operators.
static void test_unit_apply_tp
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const e_type           ///< \ref Element::type.
	);

// Interface functions ********************************************************************************************** //

void test_unit_operators_tp (struct Test_Info*const test_info)
{
	test_unit_construct_std_from_tp(test_info,"quad");
	test_unit_construct_std_from_tp(test_info,"hex");
	test_unit_construct_std_from_tp(test_info,"wedge");

	test_unit_apply_tp(test_info,"quad");
	test_unit_apply_tp(test_info,"hex");
	test_unit_apply_tp(test_info,"wedge");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void test_unit_construct_std_from_tp (struct Test_Info*const test_info, const char*const e_type)
{
	sprintf(test_info->name,"%s%s%s","Operators - construct std from tp (",e_type,")");

	bool pass         = false;
	double* tol       = NULL;
	bool* differences = NULL;

	const char*const file_name_full = constructor_file_name_unit("operators/operators_tp");

	char name_ops_tp[STRLEN_MIN] = { 0, },
	     name_op_std[STRLEN_MIN] = { 0, };

	sprintf(name_ops_tp,"%s%s",e_type,"_ops_tp");
	sprintf(name_op_std,"%s%s",e_type,"_op_std");

	const struct const_Multiarray_Matrix_d* ops_tp_r =
		constructor_file_name_const_Multiarray_Matrix_d(name_ops_tp,file_name_full); // destructed
	const struct const_Matrix_d* op_std_r =
		constructor_file_name_const_Matrix_d(name_op_std,file_name_full); // destructed
	const struct const_Matrix_d* op_std_c = constructor_op_std(ops_tp_r);   // destructed

	destructor_const_Multiarray_Matrix_d(ops_tp_r);

	tol = (double[]) { EPS, };
	differences = (bool[])
		{ diff_const_Matrix_d(op_std_r,op_std_c,tol[0]),
		};
	if (check_diff(1,differences,&pass)) {
		if (differences[0])
			print_diff_const_Matrix_d(op_std_r,op_std_c,tol[0]);
	}
	test_increment_and_print(test_info,pass);

	destructor_const_Matrix_d(op_std_r);
	destructor_const_Matrix_d(op_std_c);
}

static void test_unit_apply_tp(struct Test_Info*const test_info, const char*const e_type)
{
	sprintf(test_info->name,"%s%s%s","Operators - apply tp (",e_type,")");

	bool pass         = false;
	double* tol       = NULL;
	bool* differences = NULL;
	char var_name[STRLEN_MIN] = { 0, };

	const char*const file_name_full = constructor_file_name_unit("operators/operators_tp");

	sprintf(var_name,"%s%s",e_type,"_ops_tp");
	const struct const_Multiarray_Matrix_d* ops_tp =
		constructor_file_name_const_Multiarray_Matrix_d(var_name,file_name_full); // destructed
	struct Operator op =
		{ .ops_tp = ops_tp,
		  .op_std = constructor_op_std(ops_tp), // destructed
		  .op_csr = NULL, };

	sprintf(var_name,"%s%s",e_type,"_b");
	const struct const_Multiarray_d* b_r =
		constructor_file_name_const_Multiarray_d(var_name,file_name_full); // destructed

	sprintf(var_name,"%s%s",e_type,"_c");
	const struct const_Multiarray_d* c_r =
		constructor_file_name_const_Multiarray_d(var_name,file_name_full); // destructed

	const struct const_Multiarray_d* c_std = constructor_mm_NN1C_const_Multiarray_d(op.op_std,b_r); // destructed
	const struct const_Multiarray_d* c_tp  = constructor_mm_tp_NN1C_const_Multiarray_d(op.ops_tp,b_r); // destructed

	destructor_const_Matrix_d(op.op_std);
	destructor_const_Multiarray_Matrix_d(op.ops_tp);
	destructor_const_Multiarray_d(b_r);

	tol = (double[]) { EPS, EPS, };
	differences = (bool[])
		{ diff_const_Multiarray_d(c_r,c_std,tol[0]),
		  diff_const_Multiarray_d(c_r,c_tp,tol[1]),
		};
	if (check_diff(2,differences,&pass)) {
		if (differences[0])
			print_diff_const_Multiarray_d(c_r,c_std,tol[0]);
		if (differences[1])
			print_diff_const_Multiarray_d(c_r,c_tp,tol[1]);
	}
	test_increment_and_print(test_info,pass);

	destructor_const_Multiarray_d(c_r);
	destructor_const_Multiarray_d(c_std);
	destructor_const_Multiarray_d(c_tp);
}
