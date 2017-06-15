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
 *
 *	Notation:
 *
 *	References:
 */


// Memory moving functions (ToBeDeleted)

void convert_to_mat_V (struct S_VOLUME *const VOLUME, char const mem_op)
{
	/*
	 *	Purpose:
	 *		Used to temporarily store array data in matrix structs.
	 */

	if (mem_op == 'A') {
		struct S_ELEMENT const *const ELEMENT = get_ELEMENT_type(VOLUME->type);

		unsigned int const P      = VOLUME->P,
		                   curved = VOLUME->curved;

		unsigned int const NvnS = ELEMENT->NvnS[P],
		                   NvnI = { !curved ? ELEMENT->NvnIs[P] : ELEMENT->NvnIc[P] },
		                   NvnG = VOLUME->NvnG;

		unsigned int const Nvar = DB.Nvar,
		                   d    = DB.d;

		constructor_move2_mat('C','D',NvnS,Nvar,VOLUME->What,&VOLUME->What_M);
		constructor_move2_mat('C','D',NvnI,d*d,VOLUME->C_vI,&VOLUME->C_vI_M);
		constructor_move2_mat('C','D',NvnI,1,VOLUME->detJV_vI,&VOLUME->detJV_vI_M);
		constructor_move2_mat('C','D',NvnG,d,VOLUME->XYZ,&VOLUME->XYZ_M);
	} else if (mem_op == 'F') {
		free_NULL(VOLUME->What_M);
		free_NULL(VOLUME->C_vI_M);
		free_NULL(VOLUME->detJV_vI_M);
		free_NULL(VOLUME->XYZ_M);
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

void coef_to_values_vI_M (struct S_VDATA *const VDATA, char const coef_type, char const mem_op)
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
			VDATA->W_vI = VDATA->VOLUME->What_M;
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

			VDATA->W_vI = mm_mat_alloc('C','N','N',OPS->ChiS_vI,VOLUME->What_M); // free
		} else if (mem_op == 'F') {
			matrix_free(VDATA->W_vI);
		} else {
			EXIT_UNSUPPORTED;
		}
	}
}

void compute_flux_inviscid_M (struct S_VDATA *const VDATA, struct S_FLUX_M *const FLUXDATA, char const imex_type,
                              const char mem_op)
{
	/*
	 *	Comments:
	 *		In the case of linear advection, the flux depends on the physical coordinates in general; this is not
	 *		currently the case for other PDEs.
	 *		The fluxes and flux jacobians are stored as single matrices where the ordering is as follows:
	 *			flux:          (C)olumn-major: (N)umber of (n)odes, (d)imensions, (N)umber of (eq)uations
	 *			flux jacobian: (C)olumn-major: (N)umber of (n)odes, (d)imensions, (N)umber of (eq)uations/(var)iables
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	if (mem_op == 'A') {
		if (DB.PDE_index == PDE_ADVECTION)
			FLUXDATA->XYZ = mm_mat_alloc('C','N','N',VDATA->OPS->I_vG_vI,VDATA->VOLUME->XYZ_M); // free

		FLUXDATA->d = DB.d;
		FLUXDATA->W = VDATA->W_vI;

		size_t const d    = DB.d,
		             Nvar = DB.Nvar,
		             Neq  = DB.Neq,
		             Nn   = FLUXDATA->W->NRows;

		FLUXDATA->F = constructor1_mat_D('C',Nn,d*Neq,d); // free

		if (imex_type == 'E') {
			flux_inviscid_M(FLUXDATA);
		} else if (imex_type == 'I') {
			FLUXDATA->dFdW = constructor1_mat_D('C',Nn,d*Neq*Nvar,d); // free
			jacobian_flux_inviscid_M(FLUXDATA);
		}
	} else if (mem_op == 'F') {
		if (DB.PDE_index == PDE_ADVECTION)
			matrix_free((struct S_MATRIX *) FLUXDATA->XYZ);

		matrix_free(FLUXDATA->F);
		if (imex_type == 'I')
			matrix_free(FLUXDATA->dFdW);
	} else {
		EXIT_UNSUPPORTED;
	}
}
/*
void convert_between_rp_M (struct S_MATRIX const *const C, struct S_MATRIX *const Ap, struct S_MATRIX *const Ar
                           char const *const conv_type)
{
	//
	 *	Purpose:
	 *		Convert input (A)rray between (p)hysical and (r)eference space by multiplying by appropriate (C)ofactor
	 *		terms.
	 *
	 *	Comments:
	 *		See comments in setup_geom_factors.c for storage convection of C.
	 *		Supported conversion types:
	 *			1) "FluxToRef": 1d arrays (values), column major ordering, input/output dims: [Nn*d*Nrc], [Nn*Nrc*d]
	 *				Nrc is equal either to Neq or Neq*Nvar.
	 *				The flux in reference space is ordered such that terms are grouped by dimension so that
	 *				differentiation operators in the solver can be applied blockwise. Note that this is not the same
	 *				ordering as that used for the flux in physical space.
	 *
	 *		Certain multiplications can be avoided when computing dFrdW from dFdW based on the sparsity of the flux
	 *		Jacobian, usage of this knowledge is currently not incorporated.
	 *
	 *	Notation:
	 *		Ap  : (A)rray (p)hysical
	 *		Ar  : (A)rray (r)eference
	 *		C   : (C)ofactor terms
	 *
	 *		conv_type : (conv)ersion (type)
	 *
	 *	References:
	 //

	unsigned int const d = DB.d;

	if (strstr(conv_type,"FluxToRef")) {
		memset(Ar,0.0,Nn*Nrc*d * sizeof *Ar);
		for (size_t col = 0; col < Nrc; col++) {
		for (size_t dim1 = 0; dim1 < d; dim1++) {
			size_t const IndAr = (Nrc*dim1+col)*Nn;
			for (size_t dim2 = 0; dim2 < d; dim2++) {
				size_t const IndAp = (col*d+dim2)*Nn,
				             IndC  = (dim1*d+dim2)*Nn;
				for (size_t n = 0; n < Nn; n++)
					Ar[IndAr+n] += Ap[IndAp+n]*C[IndC+n];
			}
		}}
	} else {
		EXIT_UNSUPPORTED;
	}
}
*/
