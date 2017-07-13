// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_constructors_ELEMENT.h"

#include <stdio.h>

#include "Parameters.h"
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
	unsigned int const NP = DB.NP,
	                   d  = DB.d;

	struct S_OPS_SOLVER_DG *const DG = &ELEMENT->ops.solver.DG;

	// VOLUME
	DG->ChiS_vIs = constructor_matrix4_pointer(NP,NP,NVREFSFMAX); // free
	DG->ChiS_vIc = constructor_matrix4_pointer(NP,NP,NVREFSFMAX); // free
	DG->Ds_Weak_VV = constructor_matrix5_pointer(NP,NP,1,d); // free
	DG->Dc_Weak_VV = constructor_matrix5_pointer(NP,NP,1,d); // free
	DG->I_vGs_vIs = constructor_matrix4_pointer(NP,NP,NVREFSFMAX); // free
	DG->I_vGc_vIc = constructor_matrix4_pointer(NP,NP,NVREFSFMAX); // free

	// FACE
	DG->ChiS_fIs = constructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX); // free
	DG->ChiS_fIc = constructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX); // free
	DG->Is_Weak_FV = constructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX); // free
	DG->Ic_Weak_FV = constructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX); // free
}

static void constructors_ops_solver_HDG (struct S_ELEMENT *const ELEMENT)
{
	unsigned int const NP = DB.NP;

	struct S_OPS_SOLVER_HDG *const HDG = &ELEMENT->ops.solver.HDG;

	HDG->ChiTRS_vIs = constructor_matrix2_pointer(NP); // free
	HDG->ChiTRS_vIc = constructor_matrix2_pointer(NP); // free
	HDG->Is_FF      = constructor_matrix2_pointer(NP); // free
	HDG->Ic_FF      = constructor_matrix2_pointer(NP); // free
}
