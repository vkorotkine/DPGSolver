// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "update_VOLUMEs.h"
#include "solver.h"

#include "implicit_GradW.h"

#include "finalize_LHS.h"
#include "solver_implicit.h"
#include "output_to_paraview.h"
#include "test_code_output_to_paraview.h"

/*
 *	Purpose:
 *		Solve the Poisson equation.
 *
 *	Comments:
 *		The Navier-Stokes solver functions are called here where the inviscid flux is ignored and the viscous flux is
 *		assumed only to have dependence on Q (setting Fv_func_of_W = false).
 *
 *	Notation:
 *
 *	References:
 */

void solver_Poisson(bool PrintEnabled)
{
	if (!strstr(DB.SolverType,"Implicit"))
		EXIT_UNSUPPORTED;

	update_VOLUME_Ops();
	update_memory_VOLUMEs();

	struct S_RLHS_info RLHS_info = constructor_RLHS_info_2(PrintEnabled,DB.imex_type);
	compute_RLHS(&RLHS_info);

	Mat A = NULL;
	Vec b = NULL, x = NULL;
	KSP ksp = NULL;

	solver_implicit_linear_system(&A,&b,&x,&ksp,0,PrintEnabled);
	solver_implicit_update_What(x);

	KSPDestroy(&ksp);
	finalize_ksp(&A,&b,&x,2);

	// Update Qhat based on computed solution
	implicit_GradW(false);
// Try with explicit here instead, should be sufficient
//	explicit_GradW();

	char *const fNameOut = get_fNameOut("SolFinal_"); // free
	output_to_paraview(fNameOut);
	free(fNameOut);
}
