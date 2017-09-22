// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Macros.h"

#include "test_code_integration_equivalence_real_complex.h"
#include "test_code_integration_equivalence_algorithms.h"
#include "test_code_integration_linearization.h"
#include "test_code_integration_conv_order.h"
#include "test_integration_linearization.h"
#include "test_support.h"

#include "allocators.h"

void test_integration_Euler(int nargc, char **argv)
{
	UNUSED(argv);

	const bool run_tests_equivalence_real_complex = 1,
	           run_tests_equivalence_algorithms   = 1,
	           run_tests_linearization            = 1,
	           run_tests_conv_order               = 1;

	const size_t n_argv_new = 2;
	char** argv_new = mallocator(CHAR_T,2,STRLEN_MAX,n_argv_new); // free
	char* test_name = mallocator(CHAR_T,1,STRLEN_MAX);            // free

	// **************************************************************************************************** //
	// Real/Complex Equivalence
	// **************************************************************************************************** //
	if (run_tests_equivalence_real_complex) {
		struct S_equivalence_rc data_rc =  { .nargc     = nargc,
		                                     .argvNew   = argv_new,
		                                     .PrintName = test_name, };

		test_equivalence_real_complex(&data_rc,"Euler_n-Cylinder_HollowSection_CurvedMIXED2D");
	} else {
		test_print_warning("Euler equivalence real/complex testing currently disabled");
	}

	// **************************************************************************************************** //
	// Algorithm equivalence
	// **************************************************************************************************** //
	if (run_tests_equivalence_algorithms) {
		struct S_equivalence_algs data_alg = { .nargc     = nargc,
		                                       .argvNew   = argv_new,
		                                       .PrintName = test_name, };

		test_equivalence_algorithms(&data_alg,"Euler_n-Cylinder_HollowSection_CurvedMIXED3D_HW");
	} else {
		test_print_warning("Euler equivalence algorithms testing currently disabled");
	}

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (run_tests_linearization) {
		struct Test_Linearization data_l = { .nargc     = nargc,
		                                     .argv_new  = argv_new,
		                                     .test_name = test_name, };

#if 0
	#if 0
		#if 1
test_integration_linearization("test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED2D");
		#else
test_integration_linearization("test/Euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
		#endif
	#else
		test_integration_linearization("test/Euler/Test_Euler_PeriodicVortex_Stationary_QUAD");
	#endif
#endif
		test_linearization(&data_l,"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED2D");
		test_linearization(&data_l,"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED3D_TP");
		test_linearization(&data_l,"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED3D_HW");
	} else {
		test_print_warning("Euler linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order
	// **************************************************************************************************** //
	if (run_tests_conv_order) {
		struct S_convorder data_c = { .nargc     = nargc,
		                              .argvNew   = argv_new,
		                              .PrintName = test_name, };

		test_conv_order(&data_c,"Euler_n-Parabolic_Pipe_ToBeCurvedTRI");
		test_conv_order(&data_c,"Euler_n-Elliptic_Pipe_ToBeCurvedTRI");
		test_conv_order(&data_c,"Euler_n-Elliptic_Pipe_ToBeCurvedQUAD");

		test_conv_order(&data_c,"Euler_n-GaussianBump_CurvedQUAD");
		test_conv_order(&data_c,"Euler_n-Cylinder_HollowSection_CurvedMIXED2D");
		test_conv_order(&data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedMIXED2D");

		test_conv_order(&data_c,"Euler_PeriodicVortex_Stationary_n-Cube_QUAD");
		test_conv_order(&data_c,"Euler_PeriodicVortex_n-Cube_TRI");
		test_conv_order(&data_c,"Euler_PeriodicVortex_n-Cube_QUAD");

bool const test_3D = 0;
if (test_3D) {
		test_conv_order(&data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedTET");   // Revisit with DPG (ToBeDeleted)
		test_conv_order(&data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedHEX");   // ~Optimal
		test_conv_order(&data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedWEDGE"); // ~Optimal
} else {
		test_print_warning("3D SupersonicVortex testing is currently disabled");
}

	} else {
		test_print_warning("Euler convergence order testing currently disabled");
	}

	deallocator(argv_new,CHAR_T,2,STRLEN_MAX,n_argv_new);
	deallocator(test_name,CHAR_T,1,STRLEN_MAX);
}
