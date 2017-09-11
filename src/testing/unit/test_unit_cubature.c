// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "test_unit_cubature.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "test_base.h"
#include "test_support_matrix.h"

#include "matrix.h"

#include "cubature.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the tensor-product basis functions.
static void test_unit_cubature_tensor_product
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_cubature (struct Test_Info*const test_info)
{
	test_unit_cubature_tensor_product(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void test_unit_cubature_tensor_product (struct Test_Info*const test_info)
{
	char* test_name = "Cubature - TP";
	bool pass = true;

	const char*const file_name_full = constructor_file_name_unit("cubature/cubature_tp"); // free

	struct Matrix_d* a = constructor_file_name_Matrix_d("EQ_d1_p3_Cubature",file_name_full); // destructed
/// \todo Change this to read Cubature struct.
print_Matrix_d(a,1e-10);

	pass = false;

	test_increment_and_print(test_info,pass,test_name);
}
