// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Advection.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "solver.h"
#include "solver_Advection_functions.h"
#include "solver_implicit.h"
#include "finalize_LHS.h"

#include "output_to_paraview.h"
#include "test_code_output_to_paraview.h"

/*
 *	Purpose:
 *		Solve the linear Advection equation.
 *
 *	Comments:
 *		It is currently assumed that div (dot) b = 0.
 *
 *	Notation:
 *
 *	References:
 */

void solver_Advection(bool const PrintEnabled)
{
	struct S_RLHS_info RLHS_info = constructor_RLHS_info_2(PrintEnabled,DB.imex_type);
	compute_RLHS(&RLHS_info);

	if (DB.Method == METHOD_DG) {
		if (strstr(DB.SolverType,"Explicit")) {
			EXIT_UNSUPPORTED;
		} else if (strstr(DB.SolverType,"Implicit")) {
			Mat A = NULL;
			Vec b = NULL, x = NULL;
			KSP ksp = NULL;

			solver_implicit_linear_system(&A,&b,&x,&ksp,0,PrintEnabled);
			solver_implicit_update_What(x);

			KSPDestroy(&ksp);
			finalize_ksp(&A,&b,&x,2);

			char *const fNameOut = get_fNameOut("SolFinal_"); // free
			output_to_paraview(fNameOut);
			free(fNameOut);
		}
	} else if (DB.Method == METHOD_HDG) {

		EXIT_UNSUPPORTED;
	} else {
		EXIT_UNSUPPORTED;
	}
}
