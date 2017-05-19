// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Advection.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

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

#include "output_to_paraview.h"
#include "test_code_output_to_paraview.h"

#include "matrix_functions.h"
/*
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "S_FACE.h"
#include "Test.h"

#include "update_VOLUMEs.h"

#include "finalize_LHS.h"
#include "solver_implicit.h"

#include "array_print.h"
*/
/*
 *	Purpose:
 *		Provide functions for the Advection solver.
 *
 *	Comments:
 *		Many of the RHS terms computed are 0; they are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *
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
		printf("%d\n",PrintEnabled);
	} else if (strstr(DB.SolverType,"Implicit")) {
		implicit_VOLUME_info();
		implicit_FACE_info();

		char *const fNameOut = get_fNameOut("SolFinal_");
		output_to_paraview(fNameOut);
		free(fNameOut);

		EXIT_UNSUPPORTED;
	}
}
