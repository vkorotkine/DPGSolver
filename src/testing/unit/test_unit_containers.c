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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

/// \brief Provides unit tests for the matrix-diagonal matrix multiplication functions.
static void test_unit_matrix_mm_diag
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs unit testing for the containers (\ref test_unit_containers.c).
 *  \return 0 on success. */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 2,"Invalid number of input arguments");
	const char* test_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };
	if (strcmp(test_name,"matrix_mm") == 0)
		test_unit_matrix_mm(&test_info);
	else if (strcmp(test_name,"matrix_mv") == 0)
		test_unit_matrix_mv(&test_info);
	else if (strcmp(test_name,"matrix_mm_diag") == 0)
		test_unit_matrix_mm_diag(&test_info);
	else
		EXIT_ERROR("Invalid test name: %s\n",test_name);

	output_warning_count(&test_info);

	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void test_unit_matrix_mm (struct Test_Info*const test_info)
{
	UNUSED(test_info);
	bool pass = true;

	const char*const file_name_full = set_data_file_name_unit("containers/matrix");

	const struct const_Matrix_d* a = constructor_file_name_const_Matrix_d("a_Matrix",file_name_full), // destructed
	                           * b = constructor_file_name_const_Matrix_d("b_Matrix",file_name_full); // destructed
	struct Matrix_d* c_R = constructor_file_name_Matrix_d("c_Matrix",file_name_full); // destructed

	const struct const_Matrix_d* a_t = constructor_copy_transpose_const_Matrix_d(a,false), // destructed
	                           * b_t = constructor_copy_transpose_const_Matrix_d(b,false); // destructed

	// row major
	struct Matrix_d* c_NNR = constructor_mm_Matrix_d('N','N',1.0,a,  b,  'R'), // destructed
	               * c_TNR = constructor_mm_Matrix_d('T','N',1.0,a_t,b,  'R'), // destructed
	               * c_NTR = constructor_mm_Matrix_d('N','T',1.0,a  ,b_t,'R'), // destructed
	               * c_TTR = constructor_mm_Matrix_d('T','T',1.0,a_t,b_t,'R'); // destructed

	transpose_const_Matrix_d(a,true);
	transpose_const_Matrix_d(b,true);
	transpose_const_Matrix_d(a_t,true);
	transpose_const_Matrix_d(b_t,true);

	// col major
	struct Matrix_d* c_C   = constructor_copy_transpose_Matrix_d(c_R,true),        // destructed
	               * c_NNC = constructor_mm_Matrix_d('N','N',1.0,a,  b,  'C'), // destructed
	               * c_TNC = constructor_mm_Matrix_d('T','N',1.0,a_t,b,  'C'), // destructed
	               * c_NTC = constructor_mm_Matrix_d('N','T',1.0,a  ,b_t,'C'), // destructed
	               * c_TTC = constructor_mm_Matrix_d('T','T',1.0,a_t,b_t,'C'); // destructed

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

	assert_condition(pass);
}

static void test_unit_matrix_mv (struct Test_Info*const test_info)
{
	UNUSED(test_info);
	bool pass = true;

	const char*const file_name_full = set_data_file_name_unit("containers/matrix");

	const struct const_Matrix_d* a = constructor_file_name_const_Matrix_d("a_Matrix",file_name_full); // destructed
	const struct const_Vector_d* b = constructor_file_name_const_Vector_d("b_Vector",file_name_full); // destructed
	struct Vector_d* c = constructor_file_name_Vector_d("c_Vector",file_name_full); // destructed

	const struct const_Matrix_d* a_t = constructor_copy_transpose_const_Matrix_d(a,false); // destructed

	struct Vector_d* c_N = constructor_mv_Vector_d('N',1.0,a,  b), // destructed
	               * c_T = constructor_mv_Vector_d('T',1.0,a_t,b); // destructed
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

	assert_condition(pass);
}

static void test_unit_matrix_mm_diag (struct Test_Info*const test_info)
{
	UNUSED(test_info);
	bool pass = true;

	const char*const file_name_full = set_data_file_name_unit("containers/matrix");

	const struct const_Matrix_d* a  = constructor_file_name_const_Matrix_d("a_Matrix", file_name_full); // destructed
	const struct const_Vector_d* bl = constructor_file_name_const_Vector_d("bl_Vector",file_name_full); // destructed
	const struct const_Vector_d* br = constructor_file_name_const_Vector_d("br_Vector",file_name_full); // destructed
	const struct const_Matrix_d* cl = constructor_file_name_const_Matrix_d("cl_Matrix",file_name_full); // destructed
	const struct const_Matrix_d* cr = constructor_file_name_const_Matrix_d("cr_Matrix",file_name_full); // destructed

	const struct const_Matrix_d* a_t  = constructor_copy_transpose_const_Matrix_d(a,true);  // destructed
	const struct const_Matrix_d* cl_t = constructor_copy_transpose_const_Matrix_d(cl,true); // destructed
	const struct const_Matrix_d* cr_t = constructor_copy_transpose_const_Matrix_d(cr,true); // destructed

	const struct const_Matrix_d* cl_R = constructor_mm_diag_const_Matrix_d(1.0,a,bl,'L',false);   // destructed
	const struct const_Matrix_d* cr_R = constructor_mm_diag_const_Matrix_d(1.0,a,br,'R',false);   // destructed
	const struct const_Matrix_d* cl_C = constructor_mm_diag_const_Matrix_d(1.0,a_t,bl,'L',false); // destructed
	const struct const_Matrix_d* cr_C = constructor_mm_diag_const_Matrix_d(1.0,a_t,br,'R',false); // destructed

	struct Matrix_d* cl_R_alt = constructor_empty_Matrix_d(cl_R->layout,cl_R->ext_0,cl_R->ext_1); // destructed
	struct Matrix_d* cr_R_alt = constructor_empty_Matrix_d(cr_R->layout,cr_R->ext_0,cr_R->ext_1); // destructed
	struct Matrix_d* cl_C_alt = constructor_empty_Matrix_d(cl_C->layout,cl_C->ext_0,cl_C->ext_1); // destructed
	struct Matrix_d* cr_C_alt = constructor_empty_Matrix_d(cr_C->layout,cr_C->ext_0,cr_C->ext_1); // destructed

	mm_diag_d('L',1.0,0.0,a,  bl,cl_R_alt,false);
	mm_diag_d('R',1.0,0.0,a,  br,cr_R_alt,false);
	mm_diag_d('L',1.0,0.0,a_t,bl,cl_C_alt,false);
	mm_diag_d('R',1.0,0.0,a_t,br,cr_C_alt,false);

	destructor_const_Matrix_d(a);
	destructor_const_Matrix_d(a_t);
	destructor_const_Vector_d(bl);
	destructor_const_Vector_d(br);

	const bool* differences = (bool[])
		{ diff_const_Matrix_d(cl,cl_R,EPS),
		  diff_const_Matrix_d(cr,cr_R,EPS),
		  diff_const_Matrix_d(cl_t,cl_C,EPS),
		  diff_const_Matrix_d(cr_t,cr_C,EPS),
		  diff_const_Matrix_d(cl,(struct const_Matrix_d*)cl_R_alt,EPS),
		  diff_const_Matrix_d(cr,(struct const_Matrix_d*)cr_R_alt,EPS),
		  diff_const_Matrix_d(cl_t,(struct const_Matrix_d*)cl_C_alt,EPS),
		  diff_const_Matrix_d(cr_t,(struct const_Matrix_d*)cr_C_alt,EPS),
		};
	if (check_diff(8,differences,&pass)) {
		if (differences[0])
			print_diff_const_Matrix_d(cl,cl_R,EPS);
		if (differences[1])
			print_diff_const_Matrix_d(cr,cr_R,EPS);
		if (differences[2])
			print_diff_const_Matrix_d(cl_t,cl_C,EPS);
		if (differences[3])
			print_diff_const_Matrix_d(cr_t,cr_C,EPS);

		int ind = 4;
		print_diff_cond_const_Matrix_d(cl,  (struct const_Matrix_d*)cl_R_alt,EPS,differences[ind++]);
		print_diff_cond_const_Matrix_d(cr,  (struct const_Matrix_d*)cr_R_alt,EPS,differences[ind++]);
		print_diff_cond_const_Matrix_d(cl_t,(struct const_Matrix_d*)cl_C_alt,EPS,differences[ind++]);
		print_diff_cond_const_Matrix_d(cr_t,(struct const_Matrix_d*)cr_C_alt,EPS,differences[ind++]);
	}

	destructor_const_Matrix_d(cl);
	destructor_const_Matrix_d(cr);
	destructor_const_Matrix_d(cl_t);
	destructor_const_Matrix_d(cr_t);

	destructor_const_Matrix_d(cl_R);
	destructor_const_Matrix_d(cr_R);
	destructor_const_Matrix_d(cl_C);
	destructor_const_Matrix_d(cr_C);

	destructor_Matrix_d(cl_R_alt);
	destructor_Matrix_d(cr_R_alt);
	destructor_Matrix_d(cl_C_alt);
	destructor_Matrix_d(cr_C_alt);

	assert_condition(pass);
}
