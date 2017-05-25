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
#include "S_VOLUME.h"

#include "solver_functions.h"
#include "solver_Advection_functions.h"
#include "explicit_VOLUME_info.h"
#include "explicit_FACE_info.h"
#include "implicit_VOLUME_info.h"
#include "implicit_FACE_info.h"
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
	if (strstr(DB.SolverType,"Explicit")) {
		explicit_VOLUME_info();
		explicit_FACE_info();
		EXIT_UNSUPPORTED;
	} else if (strstr(DB.SolverType,"Implicit")) {
		if (PrintEnabled) { printf("V");  } implicit_VOLUME_info();
		if (PrintEnabled) { printf("F");  } implicit_FACE_info();

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
}
