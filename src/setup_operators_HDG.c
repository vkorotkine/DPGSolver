// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "setup_operators_HDG.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "matrix_structs.h"

#include "element_functions.h"
#include "select_functions.h"
#include "adaptation.h"
#include "matrix_functions.h"
#include "array_free.h"
#include "setup_operators_support.h"
#include "memory_constructors_matrix.h"

#include "cubature.h" // Needed only for the definitions of the structs
#include "bases.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Set up additional operators required for the HDG method.
 *
 *	Comments:
 *		The operators computed here are those required assuming that the traces are in L2 (non-conforming). As it is
 *		currently assumed that the FACE cubature nodes of the parent ELEMENTs are the same as the VOLUME cubature nodes
 *		of the current element, the operators computed below are identical to the corresponding VOLUME operators if the
 *		VOLUME and FACE cubature orders are equal.
 *
 *		For consistency of FACE trace terms with VOLUME terms which are interpolated to the FACE, the cubature order for
 *		the FACE operators is specified as that of PIf of the parent ELEMENT.
 *
 *		The I(s/c)_FF operators currently include the inverse of the cubature weights for Collocated = TRUE but this
 *		does not have an effect on the scheme.
 *
 *	Notation:
 *		Chi(Ref)(1)(2)_(3)(4)(5)[6] : (Ref)erence (optional) Basis functions (Chi) of type (1) which are (2) evaluated
 *		                              at (3) nodes of (4) type which are (5) of order [6]
 *		                              (1/4): TR(ace)S(olution), I(ntegration)
 *		                              (2/5): (s)traight, (c)urved
 *		                              (3)  : (v)olume,
 *
 *	References:
 */

static void setup_operators_HDG_std (unsigned int const EType)
{
	char const *const *const *const NodeTypeS   = (char const *const *const *const) DB.NodeTypeS,
	           *const *const *const NodeTypeIfs = (char const *const *const *const) DB.NodeTypeIfs,
	           *const *const *const NodeTypeIfc = (char const *const *const *const) DB.NodeTypeIfc;

	unsigned int const              Eclass = get_Eclass(EType),
	                   *const       PTRS   = DB.PTRS,
	                   *const *const PIfs  = (unsigned int const *const *const) DB.PIfs,
	                   *const *const PIfc  = (unsigned int const *const *const) DB.PIfc;

	struct S_ELEMENT *const ELEMENT = get_ELEMENT_type(EType);

	unsigned int const dE = ELEMENT->d;

	// Returned operators
	struct S_OPS_SOLVER_HDG *const HDG = &ELEMENT->ops.solver.HDG;

	struct S_MATRIX **ChiTRS_vIs = HDG->ChiTRS_vIs,
	                **ChiTRS_vIc = HDG->ChiTRS_vIc,
	                **Is_FF      = HDG->Is_FF,
	                **Ic_FF      = HDG->Ic_FF;

	cubature_tdef cubature;
	basis_tdef    basis;

	select_functions_cubature(&cubature,EType);
	select_functions_basis(&basis,EType);

	size_t PSMin, PSMax;
	get_PS_range(&PSMin,&PSMax);
	for (size_t P = PSMin; P <= PSMax; P++) {
		// Build transformation matrix
		struct S_CUBATURE *cub_vTRS = cub_constructor(false,false,NodeTypeS[P][Eclass],dE,PTRS[P],cubature); // free

		struct S_MATRIX *ChiRefTRS_vTRS    = basis_mat(PTRS[P],cub_vTRS,basis), // free
		                *ChiRefInvTRS_vTRS = inverse_mat(ChiRefTRS_vTRS);       // free

		struct S_MATRIX *ChiTRS_vTRS = NULL;
		if (strstr(DB.BasisType,"Modal"))
			ChiTRS_vTRS = ChiRefTRS_vTRS;
		else if (strstr(DB.BasisType,"Nodal"))
			ChiTRS_vTRS = identity_mat(cub_vTRS->Nn); // free
		else
			EXIT_UNSUPPORTED;

		cub_destructor(cub_vTRS);

		struct S_MATRIX *TTRS = mm_mat_alloc('R','N','N',ChiRefInvTRS_vTRS,ChiTRS_vTRS); // free

		matrix_free(ChiRefInvTRS_vTRS);
		matrix_free(ChiTRS_vTRS);

		// Compute returned operators (straight)
		struct S_CUBATURE *cub_vIs = cub_constructor(true,false,NodeTypeIfs[P][Eclass],dE,PIfs[P][Eclass],cubature); // free

		struct S_MATRIX *w_vIs = constructor_matrix1_move_d_1(cub_vIs->Nn,cub_vIs->w);

		struct S_MATRIX *ChiRefTRS_vIs = basis_mat(PTRS[P],cub_vIs,basis); // free

		ChiTRS_vIs[P] = mm_mat_alloc('R','N','N',ChiRefTRS_vIs,TTRS);       // keep
		Is_FF[P]      = mm_diag_mat_alloc('R','T','R',ChiTRS_vIs[P],w_vIs); // keep

		cub_destructor(cub_vIs);
		matrix_free(ChiRefTRS_vIs);


		// Compute returned operators (curved)
		struct S_CUBATURE *cub_vIc = cub_constructor(true,false,NodeTypeIfc[P][Eclass],dE,PIfc[P][Eclass],cubature); // free

		struct S_MATRIX *w_vIc = constructor_matrix1_move_d_1(cub_vIc->Nn,cub_vIc->w);

		struct S_MATRIX *ChiRefTRS_vIc = basis_mat(PTRS[P],cub_vIc,basis); // free

		ChiTRS_vIc[P] = mm_mat_alloc('R','N','N',ChiRefTRS_vIc,TTRS);       // keep
		Ic_FF[P]      = mm_diag_mat_alloc('R','T','R',ChiTRS_vIc[P],w_vIc); // keep

		cub_destructor(cub_vIc);
		matrix_free(ChiRefTRS_vIc);

		matrix_free(TTRS);
	}
}

static void setup_operators_HDG_TP(unsigned int const EType)
{
	if (EType != QUAD)
		EXIT_UNSUPPORTED;

	EXIT_UNSUPPORTED;
}

static void move_operators_to_mat_std (unsigned int const EType)
{
	/*
	 *	Purpose:
	 *		Convert required operators to mat format.
	 *
	 *	Comments:
	 *		This can be removed once setup_operators is modified to directly set up arrays in this format. (ToBeDeleted)
	 */

	struct S_ELEMENT       *const ELEMENT = get_ELEMENT_type(EType);
	struct S_OPS_SOLVER_DG *const DG      = &ELEMENT->ops.solver.DG;

	struct S_OP_RANGE op_range;
	op_range.d_range  = rd_op_t;
	op_range.ELEMENT  = ELEMENT;

	unsigned int const *const NvnS  = ELEMENT->NvnS;

	// VOLUME operators
	unsigned int const *const NvnIs = ELEMENT->NvnIs,
	                   *const NvnIc = ELEMENT->NvnIc,
	                   *const NvnGs = ELEMENT->NvnGs,
	                   *const NvnGc = ELEMENT->NvnGc;
	op_range.type_op  = 'V';

	op_range.PS_range = rP_op_t;
	op_range.Pb_range = P_op_t;
	op_range.vh_range = zero_op_t;
	move_pointers_matrix5('R','D',NvnS,NvnIs,ELEMENT->Ds_Weak_VV,DG->Ds_Weak_VV,&op_range);
	move_pointers_matrix5('R','D',NvnS,NvnIc,ELEMENT->Dc_Weak_VV,DG->Dc_Weak_VV,&op_range);

//	op_range.Pb_range = rP_op_t; // Required for assembly of sum factorized operators
//	op_range.vh_range = rhrefSF_op_t;
	move_pointers_matrix4('R','D',NvnIs,NvnS,NULL,NULL,ELEMENT->ChiS_vIs,DG->ChiS_vIs,&op_range);
	move_pointers_matrix4('R','D',NvnIc,NvnS,NULL,NULL,ELEMENT->ChiS_vIc,DG->ChiS_vIc,&op_range);

	op_range.PS_range = one_op_t;
	op_range.Pb_range = rP_op_t;
	move_pointers_matrix4('R','D',NvnIs,NvnGs,NULL,NULL,ELEMENT->I_vGs_vIs,DG->I_vGs_vIs,&op_range);

	op_range.PS_range = rP_op_t;
	move_pointers_matrix4('R','D',NvnIc,NvnGc,NULL,NULL,ELEMENT->I_vGc_vIc,DG->I_vGc_vIc,&op_range);

	// FACE operators
	unsigned int const *const *const NfnIs = (unsigned int const *const *const) ELEMENT->NfnIs,
	                   *const *const NfnIc = (unsigned int const *const *const) ELEMENT->NfnIc;
	op_range.type_op  = 'F';
	op_range.EType    = EType;

	op_range.Pb_range = rP_op_t;
	op_range.fh_range = rfh_op_t;
	op_range.e_to_e   = RfCv_op_r;
	move_pointers_matrix4('R','D',NULL,NvnS,NfnIs,NULL,ELEMENT->ChiS_fIs,DG->ChiS_fIs,&op_range);
	move_pointers_matrix4('R','D',NULL,NvnS,NfnIc,NULL,ELEMENT->ChiS_fIc,DG->ChiS_fIc,&op_range);

	op_range.e_to_e   = RvCf_op_r;
	move_pointers_matrix4('R','D',NvnS,NULL,NULL,NfnIs,ELEMENT->Is_Weak_FV,DG->Is_Weak_FV,&op_range);
	move_pointers_matrix4('R','D',NvnS,NULL,NULL,NfnIc,ELEMENT->Ic_Weak_FV,DG->Ic_Weak_FV,&op_range);
}

static void move_operators_to_mat_TP (unsigned int const EType)
{
	EXIT_UNSUPPORTED; // Add support
	printf("%d\n",EType);
}

void setup_operators_HDG(void)
{
	if (DB.Method != METHOD_HDG)
		return;

	unsigned int const d = DB.d;

	unsigned int EType;

	// LINE
	EType = LINE;
	setup_operators_HDG_std(EType);
	move_operators_to_mat_std(EType);

	if (d == 3) {
		if (is_ELEMENT_present(TET) || is_ELEMENT_present(WEDGE) || is_ELEMENT_present(PYR)) {
			EType = TRI;
			setup_operators_HDG_std(EType);
		} else if (is_ELEMENT_present(HEX) || is_ELEMENT_present(WEDGE) || is_ELEMENT_present(PYR)) {
			EType = QUAD;
			setup_operators_HDG_TP(EType);
		}
	}

	// TRI
	EType = TRI;
	if (is_ELEMENT_present(EType)) {
		move_operators_to_mat_std(EType);
	}

	// QUAD
	EType = QUAD;
	if (is_ELEMENT_present(EType)) {
		move_operators_to_mat_TP(EType);
	}
}
