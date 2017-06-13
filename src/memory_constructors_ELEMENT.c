// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_constructors_ELEMENT.h"

#include <stdio.h>

#include "S_DB.h"
#include "S_ELEMENT.h"

#include "memory_constructors_matrix.h"

/*
 *	Purpose:
 *		Allocate memory for components stored in the ELEMENT structs.
 *
 *	Comments:
 *		Move New_ELEMENT (from memory_constructors.c) to this file. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void constructors_ops            (struct S_ELEMENT *const ELEMENT);
static void constructors_ops_solver     (struct S_ELEMENT *const ELEMENT);

void constructors_ELEMENT (struct S_ELEMENT *const ELEMENT)
{
	constructors_ops(ELEMENT);
}

static void constructors_ops (struct S_ELEMENT *const ELEMENT)
{
	constructors_ops_solver(ELEMENT);
}


// Solver operators
static void constructors_ops_solver_DG  (struct S_ELEMENT *const ELEMENT);
static void constructors_ops_solver_HDG (struct S_ELEMENT *const ELEMENT);

static void constructors_ops_solver (struct S_ELEMENT *const ELEMENT)
{
	constructors_ops_solver_DG(ELEMENT);
	constructors_ops_solver_HDG(ELEMENT);
}

static void constructors_ops_solver_DG (struct S_ELEMENT *const ELEMENT)
{
	return;
	if (0) printf("%p\n",ELEMENT);
}

static void constructors_ops_solver_HDG (struct S_ELEMENT *const ELEMENT)
{
	unsigned int const NP = DB.NP;

	struct S_OPS_SOLVER_HDG *const HDG = &ELEMENT->ops.solver.HDG;

	HDG->ChiTRS_vIs = constructor2_mat(NP,NP); // free
	HDG->ChiTRS_vIc = constructor2_mat(NP,NP); // free
	HDG->Is_FF      = constructor2_mat(NP,NP); // free
	HDG->Ic_FF      = constructor2_mat(NP,NP); // free
}
