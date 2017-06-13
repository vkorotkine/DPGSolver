// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_destructors_ELEMENT.h"

#include <stdio.h>

#include "S_DB.h"
#include "S_ELEMENT.h"

#include "array_free.h"

/*
 *	Purpose:
 *		Deallocate memory for components stored in the ELEMENT structs.
 *
 *	Comments:
 *		Move memory_destructor_E (from memory_destructors.c) to this file. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void destructors_ops            (struct S_ELEMENT *const ELEMENT);
static void destructors_ops_solver     (struct S_ELEMENT *const ELEMENT);

void destructors_ELEMENT (struct S_ELEMENT *const ELEMENT)
{
	destructors_ops(ELEMENT);
}

static void destructors_ops (struct S_ELEMENT *const ELEMENT)
{
	destructors_ops_solver(ELEMENT);
}


// Solver operators
static void destructors_ops_solver_DG  (struct S_ELEMENT *const ELEMENT);
static void destructors_ops_solver_HDG (struct S_ELEMENT *const ELEMENT);

static void destructors_ops_solver (struct S_ELEMENT *const ELEMENT)
{
	destructors_ops_solver_DG(ELEMENT);
	destructors_ops_solver_HDG(ELEMENT);
}

static void destructors_ops_solver_DG (struct S_ELEMENT *const ELEMENT)
{
	return;
	if (0) printf("%p\n",ELEMENT);
}

static void destructors_ops_solver_HDG (struct S_ELEMENT *const ELEMENT)
{
	unsigned int const NP = DB.NP;

	struct S_OPS_SOLVER_HDG const* HDG = &ELEMENT->ops.solver.HDG;

	matrix_free2(NP,HDG->ChiTRS_vIs);
	matrix_free2(NP,HDG->ChiTRS_vIc);
	matrix_free2(NP,HDG->Is_FF);
	matrix_free2(NP,HDG->Ic_FF);
}
