// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_functions_HDG.h"
#include "matrix_structs.h"
#include "fluxes_structs.h"
#include "S_VOLUME.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
//#include "S_OpCSR.h" // ToBeModified (Add when using sparse matrices)

#include "element_functions.h"
#include "matrix_functions.h"
#include "memory_constructors_matrix.h"

#include "fluxes_inviscid.h"
#include "jacobian_fluxes_inviscid.h"

#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Provide HDG solver related functions.
 *
 *	Comments:
 *		Many functions contain a (mem)ory_(op)eration parameter which indicates whether memory is to be allocated or
 *		freed.
 *
 *	Notation:
 *
 *	References:
 */


// Memory moving functions (ToBeDeleted after matrix structs are adopted)

void convert_to_multiarray_V (struct S_VOLUME *const VOLUME, char const mem_op)
{
	/*
	 *	Purpose:
	 *		Used to temporarily store array data in multiarray structs.
	 */

	if (mem_op == 'A') {
		struct S_ELEMENT const *const ELEMENT = get_ELEMENT_type(VOLUME->type);

		unsigned int const P      = VOLUME->P,
		                   curved = VOLUME->curved;

		unsigned int const NvnS = ELEMENT->NvnS[P],
		                   NvnI = { !curved ? ELEMENT->NvnIs[P] : ELEMENT->NvnIc[P] },
		                   NvnG = VOLUME->NvnG;

		unsigned int const Nvar = DB.Nvar,
		                   Neq  = DB.Neq,
		                   d    = DB.d;

		VOLUME->What_MA     = constructor_multiarray1_move_d_2('C',NvnS,Nvar,VOLUME->What);
		VOLUME->C_vI_MA     = constructor_multiarray1_move_d_3('C',NvnI,d,d,VOLUME->C_vI);
		VOLUME->detJV_vI_MA = constructor_multiarray1_move_d_2('C',NvnI,1,VOLUME->detJV_vI);
		VOLUME->XYZ_MA      = constructor_multiarray1_move_d_2('C',NvnG,d,VOLUME->XYZ);
		VOLUME->RHS_MA      = constructor_multiarray1_move_d_2('C',NvnS,Nvar,VOLUME->RHS);
		VOLUME->LHS_MA      = constructor_multiarray1_move_d_4('R',NvnS,NvnS,Nvar,Neq,VOLUME->LHS);
	} else if (mem_op == 'F') {
		destructor_multiarray1_default(VOLUME->What_MA);
		destructor_multiarray1_default(VOLUME->C_vI_MA);
		destructor_multiarray1_default(VOLUME->detJV_vI_MA);
		destructor_multiarray1_default(VOLUME->XYZ_MA);
		destructor_multiarray1_default(VOLUME->RHS_MA);
		destructor_multiarray1_default(VOLUME->LHS_MA);
	} else {
		EXIT_UNSUPPORTED;
	}
}


// VOLUME solver functions

struct S_OPERATORS_V *init_mat_ops_VOLUME (struct S_VOLUME const *const VOLUME)
{
	struct S_OPERATORS_V *OPS = malloc(sizeof *OPS); // returned

	unsigned int const P = VOLUME->P;

	struct S_ELEMENT       const *const ELEMENT = get_ELEMENT_type(VOLUME->type);
	struct S_OPS_SOLVER_DG const *const DG      = &ELEMENT->ops.solver.DG;

	if (!VOLUME->curved) {
		OPS->ChiS_vI = (struct S_MATRIX const *const)        DG->ChiS_vIs[P][P][0];
		OPS->D_Weak  = (struct S_MATRIX const *const *const) DG->Ds_Weak_VV[P][P][0];

		OPS->I_vG_vI = (struct S_MATRIX const *const)        DG->I_vGs_vIs[1][P][0];
	} else {
		EXIT_UNSUPPORTED; // Add support
	}

	return OPS;
}

void coef_to_values_vI_MA (struct S_VDATA *const VDATA, char const coef_type, char const mem_op)
{
	/*
	 *	Purpose:
	 *		Interpolate VOLUME coefficients (of type 'coef_type') to VOLUME cubature nodes.
	 *
	 *	Comments:
	 *
	 *	Notation:
	 *		To avoid confusion with (C)ofactor terms, (v)olume cubature nodes are denoted with the subscript (v)olume
	 *		(I)ntegration.
	 */

	if (DB.Collocated) { // use move constructor
		if (mem_op == 'A')
			VDATA->W_vI = VDATA->VOLUME->What_MA;
		else if (mem_op == 'F')
			return;
		else
			EXIT_UNSUPPORTED;
	} else {
		if (!(coef_type == 'W'))
			EXIT_UNSUPPORTED;

		if (mem_op == 'A') {
			struct S_OPERATORS_V const *const OPS    = VDATA->OPS;
			struct S_VOLUME      const *const VOLUME = VDATA->VOLUME;

			struct S_MATRIX const What_M = constructor_matrix1_move_multiarray2_2(VOLUME->What_MA);
			struct S_MATRIX       *const W_vI_M = mm_mat_alloc('C','N','N',OPS->ChiS_vI,&What_M); // free

			VDATA->W_vI = constructor_multiarray1_move_matrix_2(W_vI_M); // free (below)
			destructor_matrix1_default_const(W_vI_M);
		} else if (mem_op == 'F') {
			multiarray_free(VDATA->W_vI);
		} else {
			EXIT_UNSUPPORTED;
		}
	}
}

void compute_flux_inviscid_MA (struct S_VDATA *const VDATA, struct S_FLUX_MA *const FLUXDATA, char const imex_type,
                               const char mem_op)
{
	/*
	 *	Comments:
	 *		In the case of linear advection, the flux depends on the physical coordinates in general; this is not
	 *		currently the case for other PDEs.
	 *		The fluxes and flux jacobians are stored as single matrices where the ordering is as specified below. To
	 *		support these "higher-dimensional" matrices, an NColsSub index is also included in the matrix struct.
	 *			flux:          (C)olumn-major: (N)umber of (n)odes, (d)imensions, (N)umber of (eq)uations
	 *			flux jacobian: (C)olumn-major: (N)umber of (n)odes, (d)imensions, (N)umber of (eq)uations/(var)iables
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	if (mem_op == 'A') {
		if (DB.PDE_index == PDE_ADVECTION) {
			struct S_MATRIX const V_XYZ_M = constructor_matrix1_move_multiarray2_2(VDATA->VOLUME->XYZ_MA);
			struct S_MATRIX       *const XYZ_M = mm_mat_alloc('C','N','N',VDATA->OPS->I_vG_vI,&V_XYZ_M); // free

			FLUXDATA->XYZ = constructor_multiarray1_move_matrix_2(XYZ_M); // free (below)
			destructor_matrix1_default_const(XYZ_M);
		}

		FLUXDATA->d = DB.d;
		FLUXDATA->W = VDATA->W_vI;

		size_t const d    = DB.d,
		             Nvar = DB.Nvar,
		             Neq  = DB.Neq,
		             Nn   = FLUXDATA->W->extents[0];

		FLUXDATA->F = constructor_multiarray1_empty4('C',Nn,d,Neq,1); // free

		if (imex_type == 'E') {
			flux_inviscid_MA(FLUXDATA);
		} else if (imex_type == 'I') {
			FLUXDATA->dFdW = constructor_multiarray1_empty4('C',Nn,d,Neq,Nvar); // free
			jacobian_flux_inviscid_MA(FLUXDATA);
		}
	} else if (mem_op == 'F') {
		if (DB.PDE_index == PDE_ADVECTION)
			matrix_free((void *) FLUXDATA->XYZ);

		multiarray_free(FLUXDATA->F);
		if (imex_type == 'I')
			multiarray_free(FLUXDATA->dFdW);
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void dot_row (struct S_MATRIX const *const A, struct S_MATRIX const *const B, struct S_MATRIX *const C)
{
	/*
	 *	Purpose:
	 *		Computes the dot product of each row of the input arrays (A, B) and adds the results to the corresponding
	 *		row of the output array (C).
	 */

	size_t const m = A->extents[0],
	             n = A->extents[1];

	double const *const A_data = A->data,
	             *const B_data = B->data;
	double       *const C_data = C->data;

	if ((A->layout != 'C') || (B->layout != 'C') || (C->layout != 'C'))
		EXIT_UNSUPPORTED;

	for (size_t j = 0; j < n; j++) {
		size_t const IndA = j*m;
		for (size_t i = 0; i < m; i++) {
			C_data[i] += A_data[IndA+i]*B_data[IndA+i];
		}
	}
}

void compute_flux_ref_MA (struct S_MULTI_ARRAY const *const C, struct S_MULTI_ARRAY const *const Ap,
                          struct S_MULTI_ARRAY **Ar, const char mem_op)
{
	/*
	 *	Purpose:
	 *		Convert input physical (A)rray to (ref)erence space by multiplying by appropriate (C)ofactor terms.
	 *
	 *	Comments:
	 *		See comments in setup_geom_factors.c for storage convection of C.
	 *
	 *		The flux in reference space is ordered such that terms are grouped by dimension so that differentiation
	 *		operators in the solver can be applied blockwise. Note that this is not the same ordering as that used for
	 *		the flux in physical space.
	 *
	 *		Certain multiplications can be avoided when computing dFrdW from dFdW based on the sparsity of the flux
	 *		Jacobian, usage of this knowledge is currently not incorporated.
	 *
	 *	Notation:
	 *		Ap  : (A)rray (p)hysical
	 *		Ar  : (A)rray (r)eference
	 *		C   : (C)ofactor terms
	 *
	 *		Nvar == Ap->extents[3] is only equal to:
	 *			1       for the flux
	 *			DB.Nvar for the flux jacobian
	 *
	 *	References:
	 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_
	 *		                 Schemes (Eq. B.3)
	 */

	if (mem_op == 'A') {
		size_t const Nn   = Ap->extents[0],
		             Ndim = Ap->extents[1],
		             Neq  = Ap->extents[2],
		             Nvar = Ap->extents[3];

		if ((Nn != C->extents[0]) || (Ndim != DB.d) || (Neq != DB.Neq))
			EXIT_UNSUPPORTED;

		struct S_MULTI_ARRAY *Ar1 = constructor_multiarray1_empty4('C',Nn,Neq,Nvar,Ndim); // free
		set_to_zero_multiarray(Ar1);

		for (size_t dim = 0; dim < Ndim; dim++) {
			struct S_MATRIX const C_sub = constructor_matrix1_move_multiarray3_2(C,dim);
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				struct S_MATRIX const Ap_sub = constructor_matrix1_move_multiarray4_2(Ap,eq,var);
				struct S_MATRIX       Ar_sub = constructor_matrix1_move_multiarray4_1(Ar1,eq,var,dim);
				dot_row(&Ap_sub,&C_sub,&Ar_sub);
			}}
		}

		*Ar = Ar1;
	} else if (mem_op == 'F') {
		multiarray_free(*Ar);
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Inviscid_Weak_MA (struct S_MULTI_ARRAY const *const Ar_vI, struct S_MULTI_ARRAY *const RLHS,
                                       char const imex_type, struct S_VDATA const *const VDATA)
{
	/*
	 *	Purpose:
	 *		Compute the inviscid contribution to the the RHS/LHS terms for the weak formulation by applying the D_Weak
	 *		operator to the input array.
	 *
	 *	Comments:
	 *		The implementation of this function is optimized for explicit runs (i.e. where only the RHS needs to be
	 *		computed). In this case, the ordering of Fr_vI in memory is such that the RHS can be computed using d BLAS3
	 *		calls.
	 *
	 *		For implicit runs, the total runtime is currently dominated by the linear system solve, and this function
	 *		was thus not optimized. If fewer BLAS3 calls are to be used or if it is desired to compute the LHS terms
	 *		using the sum factorized operators, it is required to multiply dFrdW_vI into the ChiS_vI term (as opposed to
	 *		the derivative terms). However, this results in d times the total number of flops for the computation as the
	 *		new ChiS_vI terms are then different for each dimension. Note that computing the LHS in this manner changes
	 *		the storage format of LHS terms from Neq*Nvar blocks of size NvnS*NvnS to a single block of size
	 *		NvnS*(NvnS*Neq*Nvar) which is differently ordered in memory due to it being stored in row-major ordering.
	 *		This then requires an alternate version of the finalize_LHS function for the VOLUME term (Note that this
	 *		would not be a problem if the LHS was stored in column-major ordering).
	 *
	 *		In the case of a collocated scheme, ChiS_vI == I results in the possibility of directly computing LHS terms
	 *		without any matrix-matrix products.
	 *
	 *		Certain multiplications can be avoided when computing LHS terms based on the sparsity of the flux Jacobian,
	 *		usage of this knowledge is currently not incorporated.
	 */

	struct S_OPERATORS_V const *const OPS = VDATA->OPS;

	if (imex_type == 'E') {
		unsigned int const d = DB.d;

		struct S_MATRIX RLHS_M = constructor_matrix1_move_multiarray2_2(RLHS);

//		struct S_MATRIX Ar_vI_M;

		for (size_t dim = 0; dim < d; dim++) {
			struct S_MATRIX const Ar_vI_M = constructor_matrix1_move_multiarray4_2(Ar_vI,0,dim);
			mm_mat('C','N','N',1.0,1.0,OPS->D_Weak[dim],&Ar_vI_M,&RLHS_M);
		}
	} else if (imex_type == 'I') {
		size_t const Neq  = DB.Neq,
		             Nvar = DB.Nvar;
		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
		}}
	} else {
		EXIT_UNSUPPORTED;
	}
}
