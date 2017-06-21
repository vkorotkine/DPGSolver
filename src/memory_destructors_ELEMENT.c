// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_destructors_ELEMENT.h"

#include <stdio.h>

#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"

#include "matrix_structs.h"

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
	unsigned int const NP = DB.NP,
	                   d  = DB.d;

	struct S_OPS_SOLVER_DG *const DG = &ELEMENT->ops.solver.DG;

	// VOLUME
	destructor_matrix4_pointer(NP,NP,NVREFSFMAX,DG->ChiS_vIs);
	destructor_matrix4_pointer(NP,NP,NVREFSFMAX,DG->ChiS_vIc);
	destructor_matrix5_pointer(NP,NP,1,d,DG->Ds_Weak_VV);
	destructor_matrix5_pointer(NP,NP,1,d,DG->Dc_Weak_VV);
	destructor_matrix4_pointer(NP,NP,NVREFSFMAX,DG->I_vGs_vIs);
	destructor_matrix4_pointer(NP,NP,NVREFSFMAX,DG->I_vGc_vIc);

	// FACE
	destructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX,DG->ChiS_fIs);
	destructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX,DG->ChiS_fIc);
	destructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX,DG->Is_Weak_FV);
	destructor_matrix4_pointer(NP,NP,NFREFMAX*NFMAX,DG->Ic_Weak_FV);
}

static void destructors_ops_solver_HDG (struct S_ELEMENT *const ELEMENT)
{
	unsigned int const NP = DB.NP;

	struct S_OPS_SOLVER_HDG const* HDG = &ELEMENT->ops.solver.HDG;

	destructor_matrix2_pointer(NP,HDG->ChiTRS_vIs);
	destructor_matrix2_pointer(NP,HDG->ChiTRS_vIc);
	destructor_matrix2_pointer(NP,HDG->Is_FF);
	destructor_matrix2_pointer(NP,HDG->Ic_FF);
}
