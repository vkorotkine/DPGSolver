// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_code_integration_linearization.h"
#include "test_code_integration_conv_order.h"
#include "test_support.h"

#include "adaptation.h"
#include "output_to_paraview.h"
#include "solver_explicit.h"
#include "solver_implicit.h"
#include "compute_errors.h"
#include "array_free.h"
#include "test_integration_Poisson.h"
#include "element_functions.h"
#include "array_norm.h"
#include "array_print.h"
#include "test_integration_Euler.h"
#include "update_VOLUMEs.h"
#include "update_FACEs.h"

/*
 *	Purpose:
 *		Test various aspects of the Navier-Stokes solver implementation:
 *			1) Linearization
 *			2) Optimal convergence orders
 *
 *	Comments:
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
		test_print_warning("Equivalence real vs complex not yet implemented for Navier-Stokes");
	}

	// **************************************************************************************************** //
	// Algorithm equivalence
	// **************************************************************************************************** //
	if (RunTests_equivalence_algorithms) {
		test_print_warning("Equivalence algorithms not yet implemented for Navier-Stokes");
	}

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (RunTests_linearization) {
struct S_convorder *const data_c = calloc(1 , sizeof *data_c); // free
data_c->nargc     = nargc;
data_c->argvNew   = argvNew;
data_c->PrintName = PrintName;
//if (0)
test_conv_order(data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedTRI");
free(data_c);
		struct S_linearization *const data_l = calloc(1 , sizeof *data_l); // free

		data_l->nargc     = nargc;
		data_l->argvNew   = argvNew;
		data_l->PrintName = PrintName;

		test_linearization(data_l,"NavierStokes_TRI");

		free(data_l);
	}

	// **************************************************************************************************** //
	// Convergence Order
	// **************************************************************************************************** //
	if (RunTests_conv_order) {
		struct S_convorder *const data_c = calloc(1 , sizeof *data_c); // free

		data_c->nargc     = nargc;
		data_c->argvNew   = argvNew;
		data_c->PrintName = PrintName;

		// Note: The case is converging fastest with non-collocated GL-WSH Nodal (ToBeDeleted)
		test_conv_order(data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedTRI");
//		test_conv_order(data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedQUAD");
//		test_conv_order(data_c,"NavierStokes_n-Cylinder_Hollow_ToBeCurvedMIXED2D");

		// Add tests for:
		// 1) Analytical Solution for simple shearing motion between parallel flat plates (Illingworth(1950), problem 3)
		// 2) Stokes flow over a sphere (analytical solution for compressible?)

		free(data_c);
	}

	array_free2_c(2,argvNew);
	free(PrintName);
}
