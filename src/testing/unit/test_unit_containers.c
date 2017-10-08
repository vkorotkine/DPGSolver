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

#include "test_unit_containers.h"

#include <stdlib.h>
#include <stdio.h>

#include "test_base.h"
#include "test_support.h"
#include "test_support_matrix.h"
#include "test_support_vector.h"

#include "macros.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the matrix-matrix multiplication functions.
static void test_unit_matrix_mm
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the matrix-vector multiplication functions.
static void test_unit_matrix_mv
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_containers (struct Test_Info*const test_info)
{
	test_unit_matrix_mm(test_info);
	test_unit_matrix_mv(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void test_unit_matrix_mm (struct Test_Info*const test_info)
{
	sprintf(test_info->name,"%s","Containers - Matrix mm");
	bool pass = true;

	const char*const file_name_full = constructor_file_name_unit("containers/matrix"); // free

	const struct const_Matrix_d* a = constructor_file_name_const_Matrix_d("a_Matrix",file_name_full), // destructed
	                           * b = constructor_file_name_const_Matrix_d("b_Matrix",file_name_full); // destructed
	struct Matrix_d* c_R = constructor_file_name_Matrix_d("c_Matrix",file_name_full); // destructed
	free((void*)file_name_full);

	const struct const_Matrix_d* a_t = constructor_copy_transpose_const_Matrix_d(a,false), // destructed
	                           * b_t = constructor_copy_transpose_const_Matrix_d(b,false); // destructed

	// row major
	struct Matrix_d* c_NNR = constructor_mm_Matrix_d('N','N',1.0,0.0,a,  b,  'R'), // destructed
	               * c_TNR = constructor_mm_Matrix_d('T','N',1.0,0.0,a_t,b,  'R'), // destructed
	               * c_NTR = constructor_mm_Matrix_d('N','T',1.0,0.0,a  ,b_t,'R'), // destructed
	               * c_TTR = constructor_mm_Matrix_d('T','T',1.0,0.0,a_t,b_t,'R'); // destructed

	transpose_const_Matrix_d(a,true);
	transpose_const_Matrix_d(b,true);
	transpose_const_Matrix_d(a_t,true);
	transpose_const_Matrix_d(b_t,true);

	// col major
	struct Matrix_d* c_C   = constructor_copy_transpose_Matrix_d(c_R,true),        // destructed
	               * c_NNC = constructor_mm_Matrix_d('N','N',1.0,0.0,a,  b,  'C'), // destructed
	               * c_TNC = constructor_mm_Matrix_d('T','N',1.0,0.0,a_t,b,  'C'), // destructed
	               * c_NTC = constructor_mm_Matrix_d('N','T',1.0,0.0,a  ,b_t,'C'), // destructed
	               * c_TTC = constructor_mm_Matrix_d('T','T',1.0,0.0,a_t,b_t,'C'); // destructed

	destructor_const_Matrix_d(a);
	destructor_const_Matrix_d(b);
	destructor_const_Matrix_d(a_t);
	destructor_const_Matrix_d(b_t);

	if (diff_Matrix_d(c_R,c_NNR,EPS) || diff_Matrix_d(c_R,c_TNR,EPS) ||
	    diff_Matrix_d(c_R,c_NTR,EPS) || diff_Matrix_d(c_R,c_TTR,EPS) ||
	    diff_Matrix_d(c_C,c_NNC,EPS) || diff_Matrix_d(c_C,c_TNC,EPS) ||
	    diff_Matrix_d(c_C,c_NTC,EPS) || diff_Matrix_d(c_C,c_TTC,EPS))
	{
		pass = false;

		if (diff_Matrix_d(c_R,c_NNR,EPS)) print_diff_Matrix_d(c_R,c_NNR,EPS);
		if (diff_Matrix_d(c_R,c_TNR,EPS)) print_diff_Matrix_d(c_R,c_TNR,EPS);
		if (diff_Matrix_d(c_R,c_NTR,EPS)) print_diff_Matrix_d(c_R,c_NTR,EPS);
		if (diff_Matrix_d(c_R,c_TTR,EPS)) print_diff_Matrix_d(c_R,c_TTR,EPS);
		if (diff_Matrix_d(c_C,c_NNC,EPS)) print_diff_Matrix_d(c_C,c_NNC,EPS);
		if (diff_Matrix_d(c_C,c_TNC,EPS)) print_diff_Matrix_d(c_C,c_TNC,EPS);
		if (diff_Matrix_d(c_C,c_NTC,EPS)) print_diff_Matrix_d(c_C,c_NTC,EPS);
		if (diff_Matrix_d(c_C,c_TTC,EPS)) print_diff_Matrix_d(c_C,c_TTC,EPS);
	}

	destructor_Matrix_d(c_R);
	destructor_Matrix_d(c_NNR);
	destructor_Matrix_d(c_TNR);
	destructor_Matrix_d(c_NTR);
	destructor_Matrix_d(c_TTR);

	destructor_Matrix_d(c_C);
	destructor_Matrix_d(c_NNC);
	destructor_Matrix_d(c_TNC);
	destructor_Matrix_d(c_NTC);
	destructor_Matrix_d(c_TTC);

	test_increment_and_print(test_info,pass);
}

static void test_unit_matrix_mv (struct Test_Info*const test_info)
{
	sprintf(test_info->name,"%s","Containers - Matrix mv");
	bool pass = true;

	const char*const file_name_full = constructor_file_name_unit("containers/matrix"); // free

	const struct const_Matrix_d* a = constructor_file_name_const_Matrix_d("a_Matrix",file_name_full); // destructed
	const struct const_Vector_d* b = constructor_file_name_const_Vector_d("b_Vector",file_name_full); // destructed
	struct Vector_d* c = constructor_file_name_Vector_d("c_Vector",file_name_full); // destructed
	free((void*)file_name_full);

	const struct const_Matrix_d* a_t = constructor_copy_transpose_const_Matrix_d(a,false); // destructed

	struct Vector_d* c_N = constructor_mv_Vector_d('N',1.0,0.0,a,  b), // destructed
	               * c_T = constructor_mv_Vector_d('T',1.0,0.0,a_t,b); // destructed
	destructor_const_Matrix_d(a);
	destructor_const_Matrix_d(a_t);
	destructor_const_Vector_d(b);

	if (diff_Vector_d(c,c_N,EPS) ||
	    diff_Vector_d(c,c_T,EPS))
	{
		pass = false;

		if (diff_Vector_d(c,c_N,EPS)) print_diff_Vector_d(c,c_N,EPS);
		if (diff_Vector_d(c,c_T,EPS)) print_diff_Vector_d(c,c_T,EPS);
	}

	destructor_Vector_d(c);
	destructor_Vector_d(c_N);
	destructor_Vector_d(c_T);

	test_increment_and_print(test_info,pass);
}
