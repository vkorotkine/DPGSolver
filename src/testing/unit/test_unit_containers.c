// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_unit_containers.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_base.h"
#include "test_support.h"
#include "test_support_matrix.h"
#include "test_support_vector.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

///	\brief Provides unit tests for the matrix-matrix multiplication functions.
static void test_unit_matrix_mm
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

///	\brief Provides unit tests for the matrix-vector multiplication functions.
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

/**	\brief Allocates memory and sets the name of the data file.
 *	\return See brief. */
static char* constructor_file_name
	(const char*const file_name_spec ///< The specific name of the data file (`file_name` without path or extension).
	);

static void test_unit_matrix_mm (struct Test_Info*const test_info)
{
	char* test_name = "Containers - Matrix mm";
	bool pass = true;

	const char*const file_name_full = constructor_file_name("matrix"); // free

	const struct const_Matrix_d* a = constructor_file_name_const_Matrix_d("a_Matrix",file_name_full), // destructed
	                           * b = constructor_file_name_const_Matrix_d("b_Matrix",file_name_full); // destructed
	struct Matrix_d* c_R = constructor_file_name_Matrix_d("c_Matrix",file_name_full); // destructed
	free((void*)file_name_full);

	const struct const_Matrix_d* a_t = constructor_transpose_const_Matrix_d(a,false), // destructed
	                           * b_t = constructor_transpose_const_Matrix_d(b,false); // destructed

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
	struct Matrix_d* c_C   = constructor_transpose_Matrix_d(c_R,true),             // destructed
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

	test_increment_and_print(test_info,pass,test_name);
}

static void test_unit_matrix_mv (struct Test_Info*const test_info)
{
	char* test_name = "Containers - Matrix mv";
	bool pass = true;

	const char*const file_name_full = constructor_file_name("matrix"); // free

	const struct const_Matrix_d* a = constructor_file_name_const_Matrix_d("a_Matrix",file_name_full); // destructed
	const struct const_Vector_d* b = constructor_file_name_const_Vector_d("b_Vector",file_name_full); // destructed
	struct Vector_d* c = constructor_file_name_Vector_d("c_Vector",file_name_full); // destructed
	free((void*)file_name_full);

	const struct const_Matrix_d* a_t = constructor_transpose_const_Matrix_d(a,false); // destructed

	// row major
	struct Vector_d* c_NR = constructor_mv_Vector_d('R','N',1.0,0.0,a,  b), // destructed
	               * c_TR = constructor_mv_Vector_d('R','T',1.0,0.0,a_t,b); // destructed

	transpose_const_Matrix_d(a,true);
	transpose_const_Matrix_d(a_t,true);

	// col major
	struct Vector_d* c_NC = constructor_mv_Vector_d('C','N',1.0,0.0,a,  b), // destructed
	               * c_TC = constructor_mv_Vector_d('C','T',1.0,0.0,a_t,b); // destructed

	destructor_const_Matrix_d(a);
	destructor_const_Matrix_d(a_t);
	destructor_const_Vector_d(b);

	if (diff_Vector_d(c,c_NR,EPS) ||
	    diff_Vector_d(c,c_TR,EPS) ||
	    diff_Vector_d(c,c_NC,EPS) ||
	    diff_Vector_d(c,c_TC,EPS))
	{
		pass = false;

		if (diff_Vector_d(c,c_NR,EPS)) print_diff_Vector_d(c,c_NR,EPS);
		if (diff_Vector_d(c,c_TR,EPS)) print_diff_Vector_d(c,c_TR,EPS);
		if (diff_Vector_d(c,c_NC,EPS)) print_diff_Vector_d(c,c_NC,EPS);
		if (diff_Vector_d(c,c_TC,EPS)) print_diff_Vector_d(c,c_TC,EPS);
	}

	destructor_Vector_d(c);
	destructor_Vector_d(c_NR);
	destructor_Vector_d(c_TR);
	destructor_Vector_d(c_NC);
	destructor_Vector_d(c_TC);

	test_increment_and_print(test_info,pass,test_name);
}

// Level 1 ********************************************************************************************************** //

/**	\brief Allocates memory and sets the path to the data file.
 *	\return See brief. */
static char* constructor_file_name_base ();

static char* constructor_file_name (const char*const file_name_spec)
{
	char*const file_name = constructor_file_name_base();

	strcat(file_name,"/");
	strcat(file_name,file_name_spec);
	strcat(file_name,".data");

	return file_name;
}

// Level 2 ********************************************************************************************************** //

static char* constructor_file_name_base ()
{
	char*const file_name_base = malloc(STRLEN_MAX * sizeof *file_name_base); // returned

	strcpy(file_name_base,"../testing/unit/containers");

	return file_name_base;
}
