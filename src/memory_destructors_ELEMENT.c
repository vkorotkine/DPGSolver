// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_destructors_ELEMENT.h"

#include <stdio.h>

#include "Parameters.h"
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
	unsigned int const NP = DB.NP,
	                   d  = DB.d;

	struct S_OPS_SOLVER_DG *const DG = &ELEMENT->ops.solver.DG;

	// VOLUME
	matrix_free4(NP,NP,NVREFSFMAX,DG->ChiS_vIs);
	matrix_free4(NP,NP,NVREFSFMAX,DG->ChiS_vIc);
	matrix_free5(NP,NP,1,d,DG->Ds_Weak_VV);
	matrix_free5(NP,NP,1,d,DG->Dc_Weak_VV);
	matrix_free4(NP,NP,NVREFSFMAX,DG->I_vGs_vIs);
	matrix_free4(NP,NP,NVREFSFMAX,DG->I_vGc_vIc);

	// FACE
	matrix_free4(NP,NP,NFREFMAX*NFMAX,DG->ChiS_fIs);
	matrix_free4(NP,NP,NFREFMAX*NFMAX,DG->ChiS_fIc);
	matrix_free4(NP,NP,NFREFMAX*NFMAX,DG->Is_Weak_FV);
	matrix_free4(NP,NP,NFREFMAX*NFMAX,DG->Ic_Weak_FV);
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
