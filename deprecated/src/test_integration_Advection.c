// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_integration_Advection.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Macros.h"

#include "test_code_integration_conv_order.h"
#include "test_code_integration_linearization.h"
#include "test_integration_linearization.h"
#include "test_support.h"

#include "allocators.h"

void test_integration_Advection(int nargc, char **argv)
{
	UNUSED(argv);
bool const TestHDG = 0;

	const bool run_tests_linearization = 1,
	           run_tests_conv_order    = 1;

	const size_t n_argv_new = 2;
	char** argv_new = mallocator(CHAR_T,2,STRLEN_MAX,n_argv_new); // free
	char* test_name = mallocator(CHAR_T,1,STRLEN_MAX);            // free

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (run_tests_linearization) {
		struct Test_Linearization data_l = { .nargc     = nargc,
		                                     .argv_new  = argv_new,
		                                     .test_name = test_name, };

if (TestHDG) {
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightTRI");
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightQUAD");
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightMIXED2D");
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_CurvedTRI");
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_CurvedQUAD");
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_CurvedMIXED2D");
}

		// 2D (Mixed TRI/QUAD mesh)
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_n-Cube_StraightTRI");
if (0) {
		test_integration_linearization("test/Advection/Test_Advection_Default_n-Cube_StraightTRI");
}
		test_linearization(&data_l,"test/Advection/Test_Advection_Default_n-Cube_StraightQUAD");

		test_print_warning("Advection curved element testing not yet implemented");
//		test_linearization(&data_l,"Advection_CurvedTRI");
//		test_linearization(&data_l,"Advection_ToBeCurvedTRI");

		test_print_warning("Advection 3D testing needs to be implemented");
	} else {
		test_print_warning("Advection linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
	if (run_tests_conv_order) {
		struct S_convorder data_c = { .nargc     = nargc,
		                              .argvNew   = argv_new,
		                              .PrintName = test_name, };
if (!TestHDG) {
		test_conv_order(&data_c,"test/Advection/Test_Advection_Default_n-Cube_StraightTRI");
		test_conv_order(&data_c,"test/Advection/Test_Advection_Default_n-Cube_StraightQUAD");

		test_conv_order(&data_c,"test/Advection/Test_Advection_Peterson_n-Cube_StraightTRI");
} else {
		test_conv_order(&data_c,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightTRI");
}
	} else {
		test_print_warning("Advection convergence order testing currently disabled");
	}

	deallocator(argv_new,CHAR_T,2,STRLEN_MAX,n_argv_new);
	deallocator(test_name,CHAR_T,1,STRLEN_MAX);
}
