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
 *		Test various aspects of the Euler solver implementation:
 *			- Equivalence between real and complex versions of functions;
 *			- Equivalence between running using different algorithms (with different flop counts);
 *			- Linearization;
 *			- Optimal convergence orders.
 *
 *	Comments:
 *		Optimal convergence orders for 3D curved meshes which are not associated with extruded 2D meshes has so far not
 *		been obtained. This is potentially a result of mesh regularity issues or because memory constraints do not allow
 *		for the attainment of the asymptotic regime. (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_Euler(int nargc, char **argv)
{
//	bool const (ToBeModified)
	bool RunTests_equivalence_real_complex = 0,
	     RunTests_equivalence_algorithms   = 0,
	     RunTests_linearization            = 0,
	     RunTests_conv_order               = 1;

	// ToBeDeleted after Manmeet has finished his initial verification.
	bool const PeriodicVortexOnly = 0;
	if (PeriodicVortexOnly) {
		RunTests_equivalence_real_complex = 0;
		RunTests_equivalence_algorithms   = 0;
		RunTests_linearization            = 0;
	}

	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);   // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew);  // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew);  // free
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

		test_equivalence_real_complex(data_rc,"Euler_n-Cylinder_HollowSection_CurvedMIXED2D");

		free(data_rc);
	} else {
		test_print_warning("Euler equivalence real/complex testing currently disabled");
	}

	// **************************************************************************************************** //
	// Algorithm equivalence
	// **************************************************************************************************** //
	if (RunTests_equivalence_algorithms) {
		struct S_equivalence_algs *const data_alg = calloc(1 , sizeof *data_alg); // free

		data_alg->nargc     = nargc;
		data_alg->argvNew   = argvNew;
		data_alg->PrintName = PrintName;

		test_equivalence_algorithms(data_alg,"Euler_n-Cylinder_HollowSection_CurvedMIXED3D_HW");

		free(data_alg);
	} else {
		test_print_warning("Euler equivalence algorithms testing currently disabled");
	}

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (RunTests_linearization) {
		struct S_linearization *const data_l = calloc(1 , sizeof *data_l); // free

		data_l->nargc     = nargc;
		data_l->argvNew   = argvNew;
		data_l->PrintName = PrintName;

		test_linearization(data_l,"Euler_MIXED2D");
		test_linearization(data_l,"Euler_MIXED_TET_PYR");
		test_linearization(data_l,"Euler_MIXED_HEX_WEDGE");

		free(data_l);
	} else {
		test_print_warning("Euler linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order
	// **************************************************************************************************** //
	if (RunTests_conv_order) {
		struct S_convorder *const data_c = calloc(1 , sizeof *data_c); // free

		data_c->nargc     = nargc;
		data_c->argvNew   = argvNew;
		data_c->PrintName = PrintName;

if (!PeriodicVortexOnly) {
		test_conv_order(data_c,"Euler_n-GaussianBump_CurvedQUAD");
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_CurvedMIXED2D");
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedMIXED2D");
}

bool const test_3D = 0;
if (test_3D) {
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedTET");   // Revisit with DPG (ToBeDeleted)
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedHEX");   // ~Optimal
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedWEDGE"); // ~Optimal
} else {
		test_print_warning("3D SupersonicVortex testing is currently disabled");
}

if (PeriodicVortexOnly) {
//		test_conv_order(data_c,"Euler_PeriodicVortex_Stationary_n-Cube_QUAD");
//		test_conv_order(data_c,"Euler_PeriodicVortex_n-Cube_TRI");
		test_conv_order(data_c,"Euler_PeriodicVortex_n-Cube_QUAD");
}
		free(data_c);
	} else {
		test_print_warning("Euler convergence order testing currently disabled");
	}

	array_free2_c(2,argvNew);
	free(PrintName);
}
