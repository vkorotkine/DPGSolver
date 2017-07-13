// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "compute_VOLUME_info_HDG.h"
#include "compute_FACE_info_HDG.h"

#include "explicit_FACE_info.h"
#include "compute_GradW_DG.h"
#include "compute_VOLUME_RLHS_DG.h"
#include "compute_FACE_RLHS_DG.h"

/*
 *	Purpose:
 *		Provide solver related functions.
 *
 *	Comments:
 *		This is the interface to all functions related to the solving stage of the computation:
 *			- compute_RLHS;
 *			- update_solution;
 *			- update_space.
 */

/*bool evaluate_exit_condition (const struct S_solver_info*const solver_info)
{
	if (solver_info->steady) {
		if (solver_info->linear)
			return true;

		double const EXIT_RHS_RATIO = 1e10;
//		double const maxRHS = compute_maxRHS();
double maxRHS = 0.0, maxRHS0 = 0.0;

		if ((maxRHS0/maxRHS > EXIT_RHS_RATIO) || (maxRHS < 1e1*EPS))
			return true;
	} else {
		if (time == final_time)
			return true;
	}
	return false;
}*/

void compute_final_solution (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Compute the final solution.
	 */

	for (bool finished = false; !finished; ) {
EXIT_UNSUPPORTED;
finished = false;
		compute_RLHS(solver_info);
//		finalize_RLHS(solver_info);
//		update_solution(solver_info);
//		update_space(solver_info); // For adaptation (ToBeDeleted: remove this comment)

//		finished = evaluate_exit_condition(solver_info);
	}

// Move to postprocessing function: (ToBeDeleted)
//	update_gradients();
//	if (solver_info->output)
//		output_to_paraview("SolFinal_");
}

void compute_RLHS (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Compute the (R)ight and optionally the (L)eft (H)and (S)ide residual contributions.
	 *
	 *	Comments:
	 *		The RHS includes all terms of the discretization except for the time-stepping term (i.e. d/dt What = RHS).
	 *		The LHS represents the linearization of the RHS (i.e. LHS = d(RHS)/d(What)). The LHS is then used to update
	 *		the solution using Newton's method:
	 *			1) approximate: RHS(What_exact) = 0 = RHS(What) + LHS(What)*dWhat + O(dWhat^2);
	 *			2) solve:       LHS*dWhat = -RHS(What);
	 *			3) update:      What += dWhat
	 *			4) repeat:      if LHS is a function of What, repeat until convergence.
	 */

	switch (solver_info->method) {
	case METHOD_DG: {
		compute_GradW_DG(solver_info);
		compute_VOLUME_RLHS_DG(solver_info);
		compute_FACE_RLHS_DG(solver_info);
// Delete explicit_FACE_info file (ToBeDeleted)

		break;
	} case METHOD_HDG:
// Change info to RLHS (ToBeDeleted)
		compute_VOLUME_info_HDG(solver_info);
		compute_FACE_info_HDG(solver_info);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void finalize_RLHS (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Finalize the (R)ight and (L) (H)and (S)ide terms.
	 *
	 *	Comments:
	 *		The original implementation uses redundant storage of RHS and LHS terms separated based on whether they are
	 *		computed using VOLUME or FACE integrals. This function combines these contributions if applicable.
	 *		(ToBeModified)
	 */

	switch (solver_info->method) {
	case METHOD_DG: {
// Refactor such that the RHS and LHS terms are directly store in the same memory location. (ToBeDeleted)
		EXIT_UNSUPPORTED;
		break;
	} case METHOD_HDG:
// Change info to RLHS (ToBeDeleted)
		EXIT_UNSUPPORTED;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void update_solution (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Update the solution using the appropriate explicit/implicit updating procedure.
	 *
	 *	Comments:
	 *		For the implicit solver, the LHS is then used to update the solution using Newton's method:
	 *			Taylor-Expansion: RHS(What_converged) := 0 = RHS(What) + LHS(What)*dWhat + O(dWhat^2).
	 *
	 *			1) solve:   LHS*dWhat = -RHS(What);
	 *			2) update:  What += dWhat;
	 *			3) iterate: if the LHS is a function of What (i.e. the equation is nonlinear), repeat until convergence.
	 */

	if (solver_info->imex_type == 'E') {
	} else if (solver_info->imex_type == 'I') {
/*		Mat A   = NULL;
		Vec b   = NULL,
		    x   = NULL;
		KSP ksp = NULL;

// Advection and Poisson (ToBeDeleted)
		solver_implicit_linear_system(&A,&b,&x,&ksp,0,PrintEnabled);
		solver_implicit_update_What(x);

		KSPDestroy(&ksp);
		finalize_ksp(&A,&b,&x,2);*/
	}
}



// "class" functions

struct S_solver_info constructor_solver_info (const bool display, const bool output, const bool adapt,
                                              const char imex_type, const unsigned int method)
{
	struct S_solver_info solver_info;

	solver_info.display = display;
	solver_info.output  = output;
	solver_info.adapt   = adapt;

	if (imex_type != 'E' && imex_type != 'I')
		EXIT_UNSUPPORTED;
	solver_info.imex_type = imex_type;

	if (method != METHOD_DG && method != METHOD_HDG)
		EXIT_UNSUPPORTED;
	solver_info.method = method;

	switch (DB.PDE_index) {
	case PDE_ADVECTION:
		solver_info.positivity = false;
		solver_info.symmetric  = false;
		solver_info.linear     = true;
		break;
	case PDE_POISSON:
		solver_info.positivity = false;
		solver_info.symmetric  = true;
		solver_info.linear     = true;
		break;
	case PDE_EULER:
	case PDE_NAVIERSTOKES:
		solver_info.positivity = true;
		solver_info.symmetric  = false;
		solver_info.linear     = false;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	return solver_info;
}
