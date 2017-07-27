// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_c.h"
#include "solver.h"

#include <stdlib.h>
#include <stdio.h>

#include "Macros.h"
#include "Parameters.h"

#include "compute_GradW_DG_c.h"

/*
 *	Purpose:
 *		Provide complex solver related functions.
 */

void compute_GradW_c (const struct S_solver_info*const solver_info, const char stage)
{
	if (!(stage == 'A' || stage == 'F'))
		EXIT_UNSUPPORTED;

	switch (solver_info->method) {
	case METHOD_DG:
		if (stage == 'A')
			compute_GradW_DG_c(solver_info);
		else if (stage == 'F')
			free_GradW_DG_c(solver_info);

		break;
	case METHOD_HDG:
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}
