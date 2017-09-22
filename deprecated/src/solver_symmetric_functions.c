// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_symmetric_functions.h"
#include "solver.h"
#include "S_VOLUME.h"

#include <stdlib.h>
#include <complex.h>

#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_FACE.h"

#include "explicit_info_c.h"
#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Provide solver related functions for symmetric systems.
 *
 *	Comments:
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
	 *		Premultiply RHS/LHS entries by diagonal weights for collocated schemes to recover the symmetry of the global
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

			unsigned int const NvnSR = OPSR.NvnI;

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

void correct_collocated_for_symmetry_local (struct S_LHS_info*const LHS_info)
{
	/*
	 *	Purpose:
	 *		Premultiply local LHS entry by diagonal weights for collocated schemes before adding to the petsc Mat.
	 *
	 *	Comments:
	 *		When collocation is enabled, symmetry of the global system matrix is lost for the diffusion operator due to
	 *		the premultiplication of operators by inverse cubature weights. To recover the symmetry, the RHS/LHS terms
	 *		are corrected before assembly of the global system matrix.
	 *
	 *		Some redundant operations are performed when using this function as the same diagonal scaling is applied to
	 *		several terms occupying the same position in the global system matrix. However, this is a relatively cheap
	 *		operation and it allows for the duplicated memory to no longer be needed as the data is stored directly in
	 *		the petsc Mat.
	 */

	if (!DB.Collocated)
		return;

	const unsigned int Neq  = DB.Neq,
	                   Nvar = DB.Nvar;

	const struct S_OPERATORS OPS = init_ops(LHS_info->VOLUME[0]);

	const unsigned int NvnSL = LHS_info->Nn[0],
	                   NvnSR = LHS_info->Nn[1];

	for (size_t eq  = 0; eq < Neq; eq++) {
	for (size_t var = 0; var < Nvar; var++) {
		size_t const IndLHS = (eq*Nvar+var)*NvnSL*NvnSR;
		mm_diag_d_inplace(NvnSL,NvnSR,OPS.w_vI,&LHS_info->LHS[IndLHS],1.0,'L','R');
	}}
}

void correct_collocated_for_symmetry_c (struct S_VOLUME *const VOLUME_perturbed, bool const correct_all,
                                        bool const correct_V, bool const correct_F)
{
	/*
	 *	Purpose:
	 *		Premultiply RHS_c entries by diagonal weights for collocated schemes to recover the symmetry of the global
	 *		system matrix.
	 *
	 *	Comments:
	 *		See comments in setup_operators for why this correction is necessary.
	 *
	 *	Notation:
	 *		correct_all : multiply (all) RHS terms by diagonal weights?
	 *		correct_V   : multiply (V)OLUME RHS terms?
	 *		correct_F   : multiply (F)ACE   RHS terms?
	 */

	if (!DB.Collocated)
		return;

	if (correct_all && !(correct_V || correct_F))
		EXIT_UNSUPPORTED;

	if (!(correct_V || correct_F))
		EXIT_UNSUPPORTED;

	unsigned int const Neq  = DB.Neq;

	struct S_LOCAL_MESH_ELEMENTS local_ELEMENTs;

	if (correct_V) {
		if (!correct_all)
			local_ELEMENTs = compute_local_ELEMENT_list(VOLUME_perturbed,'V');

		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (!correct_all && !is_VOLUME_in_local_list(VOLUME,&local_ELEMENTs))
				continue;

			struct S_OPERATORS OPS = init_ops(VOLUME);

			unsigned int const NvnS = OPS.NvnI;

			mm_diag_dc_inplace(NvnS,Neq,OPS.w_vI,VOLUME->RHS_c,1.0,'L','C');
		}
	}

	if (correct_F) {
		if (!correct_all)
			local_ELEMENTs = compute_local_ELEMENT_list(VOLUME_perturbed,'F');

		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			if (!correct_all && !is_FACE_in_local_list(FACE,&local_ELEMENTs))
				continue;

			struct S_OPERATORS OPSL = init_ops(FACE->VL);

			unsigned int const NvnSL = OPSL.NvnI;
			mm_diag_dc_inplace(NvnSL,Neq,OPSL.w_vI,FACE->RHSL_c,1.0,'L','C');

			if (!FACE->Boundary) {
				struct S_OPERATORS OPSR = init_ops(FACE->VR);

				unsigned int const NvnSR = OPSR.NvnI;
				mm_diag_dc_inplace(NvnSR,Neq,OPSR.w_vI,FACE->RHSR_c,1.0,'L','C');
			}
		}
	}
}
