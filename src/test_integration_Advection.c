// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_integration_Advection.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Macros.h"

#include "test_code_integration_conv_order.h"
#include "test_code_integration_linearization.h"
#include "test_support.h"

#include "array_free.h"

/*
 *	Purpose:
 *		Test various aspects of the Advection solver implementation:
 *			- Linearization;
 *			- Optimal convergence orders.
 *
 *	Comments:
 *		Convergence order tests failed (TRI/QUAD) meshes when using a GLL nodal basis (but succeeded for WSH/GL nodal
 *		and all modal). The failure was due to KSPConvergedReason = -11. Aditya tried with a pseudotimestepping code and
 *		obtained optimal orders for P2 EQ nodes.
 *		INVESTIGATE (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_Advection(int nargc, char **argv)
{
bool const TestHDG = 0;
	bool const RunTests_linearization = 1,
	           RunTests_conv_order    = 1;

	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	PrintName  = malloc(STRLEN_MAX * sizeof *PrintName); // free

	// silence
	strcpy(argvNew[0],argv[0]);

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (RunTests_linearization) {
		struct Test_Linearization *data_l = calloc(1 , sizeof *data_l); // free

		data_l->nargc     = nargc;
		data_l->argvNew   = argvNew;
		data_l->PrintName = PrintName;

if (TestHDG) {
		test_linearization(data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightTRI");
		test_linearization(data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightQUAD");
		test_linearization(data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightMIXED2D");
		test_linearization(data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_CurvedTRI");
		test_linearization(data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_CurvedQUAD");
		test_linearization(data_l,"test/Advection/Test_Advection_Default_HDG_n-Cube_CurvedMIXED2D");
}

		// 2D (Mixed TRI/QUAD mesh)
		test_linearization(data_l,"test/Advection/Test_Advection_Default_n-Cube_StraightTRI");
		test_linearization(data_l,"test/Advection/Test_Advection_Default_n-Cube_StraightQUAD");

		test_print_warning("Advection curved element testing not yet implemented");
//		test_linearization(data_l,"Advection_CurvedTRI");
//		test_linearization(data_l,"Advection_ToBeCurvedTRI");

		test_print_warning("Advection 3D testing needs to be implemented");

		free(data_l);
	} else {
		test_print_warning("Advection linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
	if (RunTests_conv_order) {
		struct S_convorder *data_c = calloc(1 , sizeof *data_c); // free

		data_c->nargc     = nargc;
		data_c->argvNew   = argvNew;
		data_c->PrintName = PrintName;
if (!TestHDG) {
		test_conv_order(data_c,"test/Advection/Test_Advection_Default_n-Cube_StraightTRI");
		test_conv_order(data_c,"test/Advection/Test_Advection_Default_n-Cube_StraightQUAD");

		test_conv_order(data_c,"test/Advection/Test_Advection_Peterson_n-Cube_StraightTRI");
} else {
		test_conv_order(data_c,"test/Advection/Test_Advection_Default_HDG_n-Cube_StraightTRI");
}
		free(data_c);
	} else {
		test_print_warning("Advection convergence order testing currently disabled");
	}

	array_free2_c(2,argvNew);
	free(PrintName);
}
