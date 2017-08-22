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
#include "test_support.h"

#include "allocators.h"

void test_integration_NavierStokes(int nargc, char **argv)
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
		struct S_equivalence_rc data_rc = { .nargc     = nargc,
		                                    .argvNew   = argv_new,
		                                    .PrintName = test_name, };

		test_equivalence_real_complex(&data_rc,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D_Collocated");
		test_equivalence_real_complex(&data_rc,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D");
	} else {
		test_print_warning("Navier-Stokes equivalence real/complex testing currently disabled");
	}

	// **************************************************************************************************** //
	// Algorithm equivalence
	// **************************************************************************************************** //
	if (run_tests_equivalence_algorithms) {
		struct S_equivalence_algs data_alg = { .nargc     = nargc,
		                                       .argvNew   = argv_new,
		                                       .PrintName = test_name, };

		test_equivalence_algorithms(&data_alg,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D");
		test_print_warning("WEDGE equivalence algorithms not being tested for Navier-Stokes");
	} else {
		test_print_warning("Navier-Stokes equivalence algorithms testing currently disabled");
	}

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (run_tests_linearization) {
		struct Test_Linearization data_l = { .nargc     = nargc,
		                                     .argv_new  = argv_new,
		                                     .test_name = test_name, };

		test_linearization(&data_l,"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedMIXED2D_Col");
		test_linearization(&data_l,"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedMIXED2D");
	} else {
		test_print_warning("Navier-Stokes linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order
	// **************************************************************************************************** //
	if (run_tests_conv_order) {
		struct S_convorder data_c = { .nargc     = nargc,
		                              .argvNew   = argv_new,
		                              .PrintName = test_name, };

		test_conv_order(&data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedTRI");
		test_conv_order(&data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedQUAD");
//		test_conv_order(&data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D");

		/**	\todo Add tests for:
		 *		1. PlaneCouette (Illingworth(1950), problem 3)
		 *		2. 3D Curved manufactured solution
		 */

//		test_conv_order(&data_c,"NavierStokes_n-Cube_StraightQUAD"); // PlaneCouette (Possibly not yet working)
	} else {
		test_print_warning("Navier-Stokes convergence order testing currently disabled");
	}

	deallocator(argv_new,CHAR_T,2,STRLEN_MAX,n_argv_new);
	deallocator(test_name,CHAR_T,1,STRLEN_MAX);
}
