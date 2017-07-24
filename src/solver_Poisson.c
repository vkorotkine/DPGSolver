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

#include "compute_GradW_DG.h"

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

	struct S_solver_info solver_info = constructor_solver_info(PrintEnabled,true,false,DB.imex_type,DB.Method);
	initialize_petsc_structs(&solver_info);
	compute_RLHS(&solver_info);

//	Mat A = NULL;
//	Vec b = NULL, x = NULL;
	KSP ksp = NULL;

	Mat A = solver_info.A;
	Vec b = solver_info.b;
	Vec x = solver_info.x;

	solver_implicit_linear_system(&A,&b,&x,&ksp,0,PrintEnabled);
	solver_implicit_update_What(x);

	KSPDestroy(&ksp);
	finalize_ksp(&A,&b,&x,2);

	// Update Qhat based on computed solution
	solver_info.display   = false;
	solver_info.imex_type = 'E';
	compute_GradW_DG(&solver_info);

	output_to_paraview("SolFinal_");
}
