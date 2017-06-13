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
	                   *const *const PIfs   = (unsigned int const *const *const) DB.PIfs,
	                   *const *const PIfc   = (unsigned int const *const *const) DB.PIfc;

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

	unsigned int PSMin, PSMax;
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

		struct S_VECTOR *w_vIs = vec_constructor_move(cub_vIs->Nn,cub_vIs->w);

		struct S_MATRIX *ChiRefTRS_vIs = basis_mat(PTRS[P],cub_vIs,basis); // free

		ChiTRS_vIs[P] = mm_mat_alloc('R','N','N',ChiRefTRS_vIs,TTRS);       // keep
		Is_FF[P]      = mm_diag_mat_alloc('R','T','R',ChiTRS_vIs[P],w_vIs); // keep

		cub_destructor(cub_vIs);
		matrix_free(ChiRefTRS_vIs);


		// Compute returned operators (curved)
		struct S_CUBATURE *cub_vIc = cub_constructor(true,false,NodeTypeIfc[P][Eclass],dE,PIfc[P][Eclass],cubature); // free

		struct S_VECTOR *w_vIc = vec_constructor_move(cub_vIc->Nn,cub_vIc->w);

		struct S_MATRIX *ChiRefTRS_vIc = basis_mat(PTRS[P],cub_vIc,basis); // free

// Likely currently overwriting allocated memory in constructors_matrix (Check this)
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
	 *		This can be removed once setup_operators is modified to directly set up arrays in this format.
	 */

printf("%d\n",EType);
EXIT_BASIC;
//	struct S_ELEMENT *const ELEMENT = get_ELEMENT_type(EType);

//		OPS->D_Weak =
//			mat_constructor2_move('R','D',ELEMENT->NvnS[P],ELEMENT->NvnIs[P],DB.d, ELEMENT->Ds_Weak_VV[P][P][0],NULL);
}

static void move_operators_to_mat_TP (unsigned int const EType)
{
	if (EType != QUAD)
		EXIT_UNSUPPORTED;

	EXIT_UNSUPPORTED;
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
			move_operators_to_mat_std(EType);
		} else if (is_ELEMENT_present(HEX) || is_ELEMENT_present(WEDGE) || is_ELEMENT_present(PYR)) {
			EType = QUAD;
			setup_operators_HDG_TP(EType);
			move_operators_to_mat_TP(EType);
		}
	}
}
