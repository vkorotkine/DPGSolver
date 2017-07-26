// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"

#include "test_code_integration_equivalence_real_complex.h"
#include "test_code_integration_equivalence_algorithms.h"
#include "test_code_integration_linearization.h"
#include "test_code_integration_conv_order.h"
#include "test_support.h"

#include "array_free.h"

/*
 *	Purpose:
 *		Test various aspects of the Navier-Stokes solver implementation:
 *			- Equivalence between real and complex versions of functions;
 *			- Equivalence between running using different algorithms (with different flop counts);
 *			- Linearization;
 *			- Optimal convergence orders.
 *
 *	Comments:
 *		Both the collocated and non-collocated version of the code are checked below as the operators are different in
 *		the two cases.
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_NavierStokes(int nargc, char **argv)
{
	bool const RunTests_equivalence_real_complex = 1,
	           RunTests_equivalence_algorithms   = 1,
	           RunTests_linearization            = 1,
	           RunTests_conv_order               = 1;

	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	PrintName  = malloc(STRLEN_MAX * sizeof *PrintName); // free

	// silence
	strcpy(argvNew[0],argv[0]);

	// **************************************************************************************************** //
	// Real/Complex Equivalence
	// **************************************************************************************************** //
	if (RunTests_equivalence_real_complex) {
		struct S_equivalence_rc *const data_rc = calloc(1 , sizeof *data_rc); // free
		data_rc->nargc     = nargc;
		data_rc->argvNew   = argvNew;
		data_rc->PrintName = PrintName;

		test_equivalence_real_complex(data_rc,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D_Collocated");
		test_equivalence_real_complex(data_rc,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D");

		free(data_rc);
	} else {
		test_print_warning("Navier-Stokes equivalence real/complex testing currently disabled");
	}

	// **************************************************************************************************** //
	// Algorithm equivalence
	// **************************************************************************************************** //
	if (RunTests_equivalence_algorithms) {
		struct S_equivalence_algs *const data_alg = calloc(1 , sizeof *data_alg); // free

		data_alg->nargc     = nargc;
		data_alg->argvNew   = argvNew;
		data_alg->PrintName = PrintName;

		test_equivalence_algorithms(data_alg,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D");
		test_print_warning("WEDGE equivalence algorithms not being tested for Navier-Stokes");

		free(data_alg);
	} else {
		test_print_warning("Navier-Stokes equivalence algorithms testing currently disabled");
	}

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (RunTests_linearization) {
		struct S_linearization *const data_l = calloc(1 , sizeof *data_l); // free

		data_l->nargc     = nargc;
		data_l->argvNew   = argvNew;
		data_l->PrintName = PrintName;

		test_linearization(data_l,"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedMIXED2D_Col");
		test_linearization(data_l,"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedMIXED2D");

		free(data_l);
	} else {
		test_print_warning("Navier-Stokes linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order
	// **************************************************************************************************** //
	if (RunTests_conv_order) {
		struct S_convorder *const data_c = calloc(1 , sizeof *data_c); // free

		data_c->nargc     = nargc;
		data_c->argvNew   = argvNew;
		data_c->PrintName = PrintName;

		test_conv_order(data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedTRI");
		test_conv_order(data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedQUAD");
//		test_conv_order(data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D");

		// Add tests for:
		// 1) PlaneCouette (Illingworth(1950), problem 3)
		// 2) 3D Curved manufactured solution

//		test_conv_order(data_c,"NavierStokes_n-Cube_StraightQUAD"); // PlaneCouette (Possibly not yet working)

		free(data_c);
	} else {
		test_print_warning("Navier-Stokes convergence order testing currently disabled");
	}

	array_free2_c(2,argvNew);
	free(PrintName);
}
