// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_symmetric_functions.h"

#include <stdlib.h>
#include <complex.h>

#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Provide solver related functions for symmetric systems.
 *
 *	Comments:
 *		Currently this function is only used for the Poisson equation.
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnI;
	double const *w_vI;
};

static struct S_OPERATORS init_ops (struct S_VOLUME const *const VOLUME)
{
	struct S_OPERATORS OPS;

	struct S_ELEMENT const *const ELEMENT = get_ELEMENT_type(VOLUME->type);

	unsigned int P = VOLUME->P;

	if (!VOLUME->curved) {
		OPS.NvnI = ELEMENT->NvnIs[P];
		OPS.w_vI = ELEMENT->w_vIs[P];
	} else {
		OPS.NvnI = ELEMENT->NvnIc[P];
		OPS.w_vI = ELEMENT->w_vIc[P];
	}

	return OPS;
}

void correct_collocated_for_symmetry (void)
{
	/*
	 *	Purpose:
	 *		Premultiply RHS/LHS entries by diagonal weights for collocated schemes to recover teh symmetry of the global
	 *		system matrix.
	 *
	 *	Comments:
	 *		The correct must be applied after calling finalize_RHS so that source terms are treated correctly if
	 *		present.
	 *
	 *		See comments in setup_operators for why this correction is necessary.
	 */

	if (!DB.Collocated)
		return;

	unsigned int const Neq  = DB.Neq,
	                   Nvar = DB.Nvar;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		struct S_OPERATORS OPS = init_ops(VOLUME);

		unsigned int const NvnS = OPS.NvnI;

		mm_diag_d_inplace(NvnS,Neq,OPS.w_vI,VOLUME->RHS,1.0,'L','C');
		for (size_t eq  = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			size_t const IndLHS = (eq*Nvar+var)*NvnS*NvnS;
			mm_diag_d_inplace(NvnS,NvnS,OPS.w_vI,&VOLUME->LHS[IndLHS],1.0,'L','R');
		}}
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		struct S_OPERATORS OPSL = init_ops(FACE->VL);

		unsigned int const NvnSL = OPSL.NvnI;

		mm_diag_d_inplace(NvnSL,Neq,OPSL.w_vI,FACE->RHSL,1.0,'L','C');
		for (size_t eq  = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			size_t const IndLHS = (eq*Nvar+var)*NvnSL*NvnSL;
			mm_diag_d_inplace(NvnSL,NvnSL,OPSL.w_vI,&FACE->LHSLL[IndLHS],1.0,'L','R');
		}}

		if (!FACE->Boundary) {
			struct S_OPERATORS OPSR = init_ops(FACE->VR);

			unsigned int const NvnSR = OPSL.NvnI;

			mm_diag_d_inplace(NvnSR,Neq,OPSR.w_vI,FACE->RHSR,1.0,'L','C');
			for (size_t eq  = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				size_t IndLHS = (eq*Nvar+var)*NvnSL*NvnSR;
				mm_diag_d_inplace(NvnSL,NvnSR,OPSL.w_vI,&FACE->LHSRL[IndLHS],1.0,'L','R');
				mm_diag_d_inplace(NvnSR,NvnSL,OPSR.w_vI,&FACE->LHSLR[IndLHS],1.0,'L','R');
				       IndLHS = (eq*Nvar+var)*NvnSR*NvnSR;
				mm_diag_d_inplace(NvnSR,NvnSR,OPSR.w_vI,&FACE->LHSRR[IndLHS],1.0,'L','R');
			}}
		}
	}
}

void correct_collocated_for_symmetry_c (void)
{
	/*
	 *	Purpose:
	 *		Premultiply RHS_c entries by diagonal weights for collocated schemes to recover the symmetry of the global
	 *		system matrix.
	 *
	 *	Comments:
	 *		See comments in setup_operators for why this correction is necessary.
	 */

	if (!DB.Collocated)
		return;

	unsigned int const Neq  = DB.Neq;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		struct S_OPERATORS OPS = init_ops(VOLUME);

		unsigned int const NvnS = OPS.NvnI;

		mm_diag_dc_inplace(NvnS,Neq,OPS.w_vI,VOLUME->RHS_c,1.0,'L','C');
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		struct S_OPERATORS OPSL = init_ops(FACE->VL);

		unsigned int const NvnSL = OPSL.NvnI;

		mm_diag_dc_inplace(NvnSL,Neq,OPSL.w_vI,FACE->RHSL_c,1.0,'L','C');

		if (!FACE->Boundary) {
			struct S_OPERATORS OPSR = init_ops(FACE->VR);

			unsigned int const NvnSR = OPSL.NvnI;

			mm_diag_dc_inplace(NvnSR,Neq,OPSR.w_vI,FACE->RHSR_c,1.0,'L','C');
		}
	}
}
