// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
// \file

#include "test_integration_Poisson.h"

#include <stdlib.h>
#include <stdio.h>

#include "Macros.h"
#include "Parameters.h"

#include "test_code_integration_conv_order.h"
#include "test_code_integration_linearization.h"
#include "test_support.h"

#include "allocators.h"

void test_integration_Poisson(int nargc, char **argv)
{
	UNUSED(argv);

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

		// 1D
		test_linearization(&data_l,"test/Poisson/Test_Poisson_n-Cube_StraightLINE_Col");
		test_linearization(&data_l,"test/Poisson/Test_Poisson_n-Cube_StraightLINE");

		// 2D (Mixed TRI/QUAD mesh)
		test_linearization(&data_l,"test/Poisson/Test_Poisson_n-Ball_HollowSection_CurvedMIXED2D");

		// 3D (TET mesh)
//		test_linearization(&data_l,"Poisson_TET");
		test_print_warning("Poisson 3D testing needs to be updated");
		// Revisit when parametrized 3D geometry is available. (ToBeDeleted)
	} else {
		test_print_warning("Poisson linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
	if (run_tests_conv_order) {
		struct S_convorder data_c = { .nargc     = nargc,
		                              .argvNew   = argv_new,
		                              .PrintName = test_name, };

		test_conv_order(&data_c,"Poisson_n-Cube_LINE");

		test_conv_order(&data_c,"Poisson_n-Ellipsoid_HollowSection_TRI_extended");
		test_conv_order(&data_c,"Poisson_n-Ellipsoid_HollowSection_TRI");
		test_conv_order(&data_c,"Poisson_n-Ellipsoid_HollowSection_ToBeCurvedTRI");
		test_conv_order(&data_c,"Poisson_n-Ellipsoid_HollowSection_ToBeCurvedQUAD");
		test_conv_order(&data_c,"Poisson_n-Ellipsoid_HollowSection_QUAD");
		test_conv_order(&data_c,"Poisson_n-Ellipsoid_HollowSection_MIXED2D");
	} else {
		test_print_warning("Poisson convergence order testing currently disabled");
	}

	deallocator(argv_new,CHAR_T,2,STRLEN_MAX,n_argv_new);
	deallocator(test_name,CHAR_T,1,STRLEN_MAX);
}
