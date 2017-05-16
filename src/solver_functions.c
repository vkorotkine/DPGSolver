// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_OpCSR.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "sum_factorization.h"

#include "exact_solutions.h"
#include "variable_functions.h"
#include "boundary_conditions.h"
#include "fluxes_structs.h"
#include "fluxes_inviscid.h"
#include "fluxes_viscous.h"
#include "jacobian_boundary_conditions.h"
#include "jacobian_fluxes_inviscid.h"
#include "jacobian_fluxes_viscous.h"

#include "array_swap.h"
#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Provide solver related functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

// **************************************************************************************************** //
// General functions
// **************************************************************************************************** //

void manage_solver_memory(struct S_DATA *const DATA, char const mem_op, char const mem_type)
{
	/*
	 *	Purpose:
	 *		Allocate or free memory associated with solver functions.
	 *
	 *	Comments:
	 *		The Poisson and standard Viscous flux memory treatments are identical excepting that dnFNumdW is not needed
	 *		in the Poisson case.
	 *
	 *	Notation:
	 *		mem_op    : (mem)ory (op)eration.         Options: 'A'llocate, 'F'ree
	 *		feature   :                               Options: 'V'OLUME, 'F'ACE.
	 *		mem_type  : (mem)ory (type).
	 *		            Options: 'W' (Solution), 'Q' (Gradients), 'S' (Solution flux), 'I' (Inviscid flux),
	 *		                     'V' (Viscous flux), 'P' (Poisson Viscous flux)
	 *		imex_type : (im)plicit (ex)plicit (type). Options: 'I'mplicit, 'E'xplicit
	 */

	char const feature   = DATA->feature,
	           imex_type = DATA->imex_type;

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_VDATA         *const VDATA     = DATA->VDATA;
	struct S_FDATA         *const FDATAL    = DATA->FDATAL,
	                       *const FDATAR    = DATA->FDATAR;
	struct S_FLUX          *const FLUXDATA  = DATA->FLUXDATA;
	struct S_NumericalFlux *const NFluxData = DATA->NFluxData;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	if (mem_op == 'A') {
		if (feature == 'V') {
			unsigned int const NvnI = VDATA->OPS[0]->NvnI;

			if (mem_type == 'W') {
				VDATA->W_vI = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI)); // keep
			} else if (mem_type == 'Q') {
				VDATA->Q_vI = malloc(d * sizeof *(VDATA->Q_vI)); // keep
				for (size_t dim = 0; dim < d; dim++)
					VDATA->Q_vI[dim] = malloc(NvnI*Nvar * sizeof *(VDATA->Q_vI[dim])); // keep
			} else if (mem_type == 'I' || mem_type == 'V') {
				FLUXDATA->F  = malloc(NvnI*d*Neq * sizeof *(FLUXDATA->F));  // keep
				FLUXDATA->Fr = malloc(NvnI*Neq*d * sizeof *(FLUXDATA->Fr)); // keep

				if (imex_type == 'I') {
			       FLUXDATA->dFdW  = malloc(NvnI*d*Nvar*Neq * sizeof *(FLUXDATA->dFdW));  // keep
			       FLUXDATA->dFrdW = malloc(NvnI*Nvar*Neq*d * sizeof *(FLUXDATA->dFrdW)); // keep
				   if (mem_type == 'V') {
					   FLUXDATA->dFdQ  = malloc(d * sizeof *(FLUXDATA->dFdQ));  // keep
					   FLUXDATA->dFrdQ = malloc(d * sizeof *(FLUXDATA->dFrdQ)); // keep
					   for (size_t dim = 0; dim < d; dim++) {
						   FLUXDATA->dFdQ[dim]  = malloc(NvnI*d*Nvar*Neq * sizeof *(FLUXDATA->dFdQ[dim]));  // keep
						   FLUXDATA->dFrdQ[dim] = malloc(NvnI*Nvar*Neq*d * sizeof *(FLUXDATA->dFrdQ[dim])); // keep
					   }
				   }
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (feature == 'F') {
			struct S_FACE              *const        FACE = (struct S_FACE *const) FDATAL->FACE;
			struct S_OPERATORS_F const *const *const OPSF = FDATAL->OPS;
			unsigned int const IndFType = FDATAL->IndFType,
			                   NfnI     = OPSF[IndFType]->NfnI;

			if (mem_type == 'W') {
				FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // keep
				FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // keep
			} else if (mem_type == 'Q') {
				FDATAL->Qp_fIL = malloc(d * sizeof *(FDATAL->Qp_fIL)), // keep
				FDATAR->Qp_fIL = malloc(d * sizeof *(FDATAR->Qp_fIL)); // keep
				for (size_t dim = 0; dim < d; dim++) {
					FDATAL->Qp_fIL[dim] = malloc(NfnI*Nvar * sizeof *(FDATAL->Qp_fIL[dim])); // keep
					FDATAR->Qp_fIL[dim] = malloc(NfnI*Nvar * sizeof *(FDATAR->Qp_fIL[dim])); // keep
				}
			} else if (mem_type == 'S') {
				NFluxData->nSolNum_fI = malloc(d * sizeof *(NFluxData->nSolNum_fI)); // keep
				for (size_t dim = 0; dim < d; dim++)
					NFluxData->nSolNum_fI[dim] = malloc(NfnI*Neq * sizeof *(NFluxData->nSolNum_fI[dim])); // keep

				if (imex_type == 'I') {
					NFluxData->dnSolNumdWL_fI = malloc(d * sizeof *(NFluxData->dnSolNumdWL_fI)); // keep
					NFluxData->dnSolNumdWR_fI = malloc(d * sizeof *(NFluxData->dnSolNumdWR_fI)); // keep
					for (size_t dim = 0; dim < d; dim++) {
						if (FACE->Boundary) {
							NFluxData->dnSolNumdWL_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnSolNumdWL_fI[dim])); // keep
							NFluxData->dnSolNumdWR_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnSolNumdWR_fI[dim])); // keep
						} else {
							NFluxData->dnSolNumdWL_fI[dim] = malloc(NfnI * sizeof *(NFluxData->dnSolNumdWL_fI[dim])); // keep
							NFluxData->dnSolNumdWR_fI[dim] = malloc(NfnI * sizeof *(NFluxData->dnSolNumdWR_fI[dim])); // keep
						}
					}
				}
			} else if (mem_type == 'I') {
				NFluxData->nFluxNum_fI = malloc(NfnI*Neq * sizeof *(NFluxData->nFluxNum_fI)); // keep

				if (imex_type == 'I') {
					NFluxData->dnFluxNumdWL_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxNumdWL_fI)); // keep
					NFluxData->dnFluxNumdWR_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxNumdWR_fI)); // keep
				}
			} else if (mem_type == 'V') {
				NFluxData->nFluxViscNum_fI = malloc(NfnI*Neq * sizeof *(NFluxData->nFluxViscNum_fI)); // keep

				if (imex_type == 'I') {
					NFluxData->dnFluxViscNumdWL_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdWL_fI)); // keep
					NFluxData->dnFluxViscNumdWR_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdWR_fI)); // keep
					NFluxData->dnFluxViscNumdQL_fI = malloc(d * sizeof *(NFluxData->dnFluxViscNumdQL_fI)); // keep
					NFluxData->dnFluxViscNumdQR_fI = malloc(d * sizeof *(NFluxData->dnFluxViscNumdQR_fI)); // keep
					for (size_t dim = 0; dim < d; dim++) {
						NFluxData->dnFluxViscNumdQL_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdQL_fI[dim])); // keep
						NFluxData->dnFluxViscNumdQR_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdQR_fI[dim])); // keep
					}
				}
			} else if (mem_type == 'P') {
				NFluxData->nFluxViscNum_fI = malloc(NfnI*Neq * sizeof *(NFluxData->nFluxViscNum_fI)); // keep

				if (imex_type == 'I') {
					NFluxData->dnFluxViscNumdQL_fI = malloc(d * sizeof *(NFluxData->dnFluxViscNumdQL_fI)); // keep
					NFluxData->dnFluxViscNumdQR_fI = malloc(d * sizeof *(NFluxData->dnFluxViscNumdQR_fI)); // keep
					for (size_t dim = 0; dim < d; dim++) {
						NFluxData->dnFluxViscNumdQL_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdQL_fI[dim])); // keep
						NFluxData->dnFluxViscNumdQR_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdQR_fI[dim])); // keep
					}
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (mem_op == 'F') {
		if (feature == 'V') {
			if (mem_type == 'W') {
				free(VDATA->W_vI);
			} else if (mem_type == 'Q') {
				array_free2_d(d,VDATA->Q_vI);
			} else if (mem_type == 'I' || mem_type == 'V') {
				free(FLUXDATA->F);
				free(FLUXDATA->Fr);
				if (imex_type == 'I') {
			       free(FLUXDATA->dFdW);
			       free(FLUXDATA->dFrdW);
				   if (mem_type == 'V') {
					   array_free2_d(d,FLUXDATA->dFdQ);
					   array_free2_d(d,FLUXDATA->dFrdQ);
				   }
				}
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (feature == 'F') {
			if (mem_type == 'W') {
				free(FDATAL->W_fIL);
				free(FDATAR->W_fIL);
			} else if (mem_type == 'Q') {
				array_free2_d(d,FDATAL->Qp_fIL);
				array_free2_d(d,FDATAR->Qp_fIL);
			} else if (mem_type == 'S') {
				array_free2_d(d,NFluxData->nSolNum_fI);

				if (imex_type == 'I') {
					array_free2_d(d,NFluxData->dnSolNumdWL_fI);
					array_free2_d(d,NFluxData->dnSolNumdWR_fI);
				}
			} else if (mem_type == 'I') {
				free(NFluxData->nFluxNum_fI);

				if (imex_type == 'I') {
					free(NFluxData->dnFluxNumdWL_fI);
					free(NFluxData->dnFluxNumdWR_fI);
				}
			} else if (mem_type == 'V') {
				free(NFluxData->nFluxViscNum_fI);

				if (imex_type == 'I') {
					free(NFluxData->dnFluxViscNumdWL_fI);
					free(NFluxData->dnFluxViscNumdWR_fI);
					array_free2_d(d,NFluxData->dnFluxViscNumdQL_fI);
					array_free2_d(d,NFluxData->dnFluxViscNumdQR_fI);
				}
			} else if (mem_type == 'P') {
				free(NFluxData->nFluxViscNum_fI);

				if (imex_type == 'I') {
					array_free2_d(d,NFluxData->dnFluxViscNumdQL_fI);
					array_free2_d(d,NFluxData->dnFluxViscNumdQR_fI);
				}
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

// **************************************************************************************************** //
// VOLUME functions
// **************************************************************************************************** //

double *compute_Dxyz(struct S_Dxyz *DxyzInfo, unsigned int d)
{
	/*
	 *	Purpose:
	 *		Compute physical derivative operator matrices using the chain rule.
	 *
	 *	Comments:
	 *		Note the ordering of C specified in the comments of setup_geom_factors.
	 *
	 *		Including the cofactor terms in the differentiation operators as opposed to the fluxes results in a slightly
	 *		more clear linearization but a significant cost increase in assembly as the number of entries which must be
	 *		multiplied scales with the number of basis functions. If the cofactor terms are included in the flux, they
	 *		will always multiply only d flux terms. Furthermore, inclusion in the flux terms would allow for
	 *		vectorization of this function. A possible drawback is the required enforcement of the geometric
	 *		conservation if the strong form of the scheme is used. INVESTIGATE. (ToBeModified)
	 *
	 *	References:
	 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_
	 *		                 Schemes (eq. B.2)
	 */

	unsigned int const Nbf  = DxyzInfo->Nbf,
	                   Nn   = DxyzInfo->Nn,
	                   dim1 = DxyzInfo->dim;

	double const *const *const D = DxyzInfo->D,
	             *const        C = DxyzInfo->C;

	double *const Dxyz = calloc(Nbf*Nn, sizeof *Dxyz); // keep
	for (size_t dim2 = 0; dim2 < d; dim2++) {
		size_t const IndC = (dim1+dim2*d)*Nn;
		mm_diag_d(Nbf,Nn,&C[IndC],D[dim2],Dxyz,1.0,1.0,'R','R');
	}

	return Dxyz;
}

void init_ops_VOLUME(struct S_OPERATORS_V *const OPS, struct S_VOLUME const *const VOLUME, unsigned int const IndClass)
{
	// Initialize DB Parameters
	const unsigned int *const *const *const SF_BE = (const unsigned int *const *const *const) DB.SF_BE;

	// Standard datatypes
	unsigned int P, type, curved, Eclass;
	struct S_ELEMENT *ELEMENT, *ELEMENT_SF;

	// silence
	ELEMENT_SF = NULL;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;
	Eclass = VOLUME->Eclass;

	ELEMENT = get_ELEMENT_type(type);
	if ((Eclass == C_TP && SF_BE[P][0][0]) || (Eclass == C_WEDGE && SF_BE[P][1][0]))
		ELEMENT_SF = ELEMENT->ELEMENTclass[IndClass];
	else
		ELEMENT_SF = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[P];
	OPS->NvnS_SF = ELEMENT_SF->NvnS[P];
	if (!curved) {
		OPS->NvnI    = ELEMENT->NvnIs[P];
		OPS->NvnI_SF = ELEMENT_SF->NvnIs[P];

		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P][0];
		OPS->D_Weak  = (double const *const *const) ELEMENT->Ds_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT->Is_Weak_VV[P][P][0];

		OPS->ChiS_vI_SF = ELEMENT_SF->ChiS_vIs[P][P][0];
		OPS->D_Weak_SF  = (double const *const *const) ELEMENT_SF->Ds_Weak_VV[P][P][0];
		OPS->I_Weak_SF  = ELEMENT_SF->Is_Weak_VV[P][P][0];

		OPS->D_Weak_sp = (struct S_OpCSR const *const *const) ELEMENT->Ds_Weak_VV_sp[P][P][0];
	} else {
		OPS->NvnI    = ELEMENT->NvnIc[P];
		OPS->NvnI_SF = ELEMENT_SF->NvnIc[P];

		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];
		OPS->D_Weak  = (double const *const *const) ELEMENT->Dc_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT->Ic_Weak_VV[P][P][0];

		OPS->ChiS_vI_SF = ELEMENT_SF->ChiS_vIc[P][P][0];
		OPS->D_Weak_SF  = (double const *const *const) ELEMENT_SF->Dc_Weak_VV[P][P][0];
		OPS->I_Weak_SF  = ELEMENT_SF->Ic_Weak_VV[P][P][0];

		OPS->D_Weak_sp = (struct S_OpCSR const *const *const) ELEMENT->Dc_Weak_VV_sp[P][P][0];
	}
}

void init_VDATA(struct S_VDATA *const VDATA, struct S_VOLUME const *const VOLUME)
{
	VDATA->VOLUME = VOLUME;

	VDATA->P      = VOLUME->P;
	VDATA->Eclass = VOLUME->Eclass;
	init_ops_VOLUME((struct S_OPERATORS_V *const) VDATA->OPS[0],VOLUME,0);
	if (VOLUME->type == WEDGE)
		init_ops_VOLUME((struct S_OPERATORS_V *const) VDATA->OPS[1],VOLUME,1);
}

void coef_to_values_vI(struct S_VDATA const *const VDATA, char const coef_type)
{
	/*
	 *	Purpose:
	 *		Interpolate VOLUME coefficients (of type 'coef_type') to VOLUME cubature nodes.
	 *
	 *	Comments:
	 *		Various options are available for the computation:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
	 *
	 *	Notation:
	 *		To avoid confusion with (C)ofactor terms, (v)olume cubature nodes are denoted with the subscript (v)olume
	 *		(I)ntegration.
	 */

	if (DB.Collocated) {
		EXIT_UNSUPPORTED; // Set A_vI = VOLUME->Ahat in calling function.
	} else {
		if (!(coef_type == 'W' || coef_type == 'Q'))
			EXIT_UNSUPPORTED;

		unsigned int const d      = DB.d,
		                   Nvar   = DB.Nvar,
		                   P      = VDATA->P,
		                   Eclass = VDATA->Eclass,
	                       *const *const *const SF_BE = (const unsigned int *const *const *const) DB.SF_BE;

		struct S_OPERATORS_V const *const *const OPS    = (struct S_OPERATORS_V const *const *const) VDATA->OPS;
		struct S_VOLUME      const *const        VOLUME = VDATA->VOLUME;

		unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX];
		double const *OP[DMAX];

		if (Eclass == C_TP && SF_BE[P][0][0]) {
			get_sf_parameters(OPS[0]->NvnS_SF,OPS[0]->NvnI_SF,OPS[0]->ChiS_vI_SF,0,0,NULL,NIn,NOut,OP,d,3,Eclass);

			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 0;

			if (coef_type == 'W') {
				sf_apply_d(VOLUME->What,VDATA->W_vI,NIn,NOut,Nvar,OP,Diag,d);
			} else if (coef_type == 'Q') {
				for (size_t dim = 0; dim < d; dim++)
					sf_apply_d(VOLUME->Qhat[dim],VDATA->Q_vI[dim],NIn,NOut,Nvar,OP,Diag,d);
			}
		} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
			get_sf_parameters(OPS[0]->NvnS_SF,OPS[0]->NvnI_SF,OPS[0]->ChiS_vI_SF,
			                  OPS[1]->NvnS_SF,OPS[1]->NvnI_SF,OPS[1]->ChiS_vI_SF,NIn,NOut,OP,d,3,Eclass);

			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 0;
			Diag[1] = 2; // TRI operator (dim = 0) takes care of dim = 1 as well.

			if (coef_type == 'W') {
				sf_apply_d(VOLUME->What,VDATA->W_vI,NIn,NOut,Nvar,OP,Diag,d);
			} else if (coef_type == 'Q') {
				for (size_t dim = 0; dim < d; dim++)
					sf_apply_d(VOLUME->Qhat[dim],VDATA->Q_vI[dim],NIn,NOut,Nvar,OP,Diag,d);
			}
		} else {
			if (coef_type == 'W') {
				mm_CTN_d(OPS[0]->NvnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_vI,VOLUME->What,VDATA->W_vI);
			} else if (coef_type == 'Q') {
				for (size_t dim = 0; dim < d; dim++)
					mm_CTN_d(OPS[0]->NvnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_vI,VOLUME->Qhat[dim],VDATA->Q_vI[dim]);
			}
		}
	}
}

void convert_between_rp(unsigned int const Nn, unsigned int const Nrc, double const *const C, double *const Ap,
                        double *const Ar, char const *const conv_type)
{
	/*
	 *	Purpose:
	 *		Convert input (A)rray between (p)hysical and (r)eference space by multiplying by appropriate (C)ofactor
	 *		terms.
	 *
	 *	Comments:
	 *		See comments in setup_geom_factors.c for storage convection of C.
	 *		Supported conversion types:
	 *			1) "FluxToRef": 1d arrays, column major ordering, input/output dims: [Nn*d*Nrc], [Nn*Nrc*d]
	 *				The flux in reference space is ordered such that terms are grouped by dimension so that
	 *				differentiation operators in the solver can be applied blockwise. Note that this is not the same
	 *				ordering as that used for the flux in physical space.
	 *
	 *		Certain multiplications can be avoided when computing dFrdW from dFdW based on the sparsity of the flux
	 *		Jacobian, usage of this knowledge is currently not incorporated.
	 *
	 *	Notation:
	 *		Nn  : (N)umber of (n)odes
	 *		Nrc : (N)umber of (r)ows/(c)olumns
	 *		Ap  : (A)rray (p)hysical
	 *		Ar  : (A)rray (r)eference
	 *		C   : (C)ofactor terms
	 *
	 *		conv_type : (conv)ersion (type)
	 *
	 *	References:
	 */

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

void finalize_VOLUME_Inviscid_Weak(unsigned int const Nrc, double const *const Ar_vI, double *const RLHS,
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
	 *		Various options are available for performing the computation:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
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
	 *
	 *		If it is only desired to gain an understanding of what this function is doing is, it is suffient to look
	 *		only at the standard approach (i.e. the last condition in the if/else chain).
	 */

	struct S_OPERATORS_V const *const *const OPS = (struct S_OPERATORS_V const *const *const) VDATA->OPS;

	unsigned int const d          = DB.d,
	                   Collocated = DB.Collocated,
	                   P          = VDATA->P,
	                   Eclass     = VDATA->Eclass,
	                   NvnS       = OPS[0]->NvnS,
	                   NvnI       = OPS[0]->NvnI,
	                   *const *const *const SF_BE = (const unsigned int *const *const *const) DB.SF_BE;

	if (imex_type == 'E') {
		if (Eclass == C_TP && SF_BE[P][0][0]) {
			unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX];
			double const *OP[DMAX];

			double *const DAr = malloc(NvnS*Nrc * sizeof *DAr); // free
			for (size_t dim1 = 0; dim1 < d; dim1++) {
				get_sf_parameters(OPS[0]->NvnI_SF,OPS[0]->NvnS_SF,OPS[0]->I_Weak_SF,
				                  OPS[0]->NvnI_SF,OPS[0]->NvnS_SF,OPS[0]->D_Weak_SF[0],NIn,NOut,OP,d,dim1,Eclass);

				if (Collocated) {
					for (size_t dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 2;
					Diag[dim1] = 0;
				} else {
					for (size_t dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;
				}

				sf_apply_d(&Ar_vI[NvnI*Nrc*dim1],DAr,NIn,NOut,Nrc,OP,Diag,d);
				for (size_t i = 0, iMax = NvnS*Nrc; i < iMax; i++)
					RLHS[i] += DAr[i];
			}
			free(DAr);
		} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
			unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX];
			double const *OP[DMAX], *OP0, *OP1;

			double *const DAr = malloc(NvnS*Nrc * sizeof *DAr); // free
			for (size_t dim1 = 0; dim1 < d; dim1++) {
				if (dim1 < 2) OP0 = OPS[0]->D_Weak_SF[dim1], OP1 = OPS[1]->I_Weak_SF;
				else          OP0 = OPS[0]->I_Weak_SF,       OP1 = OPS[1]->D_Weak_SF[0];
				get_sf_parameters(OPS[0]->NvnI_SF,OPS[0]->NvnS_SF,OP0,
				                  OPS[1]->NvnI_SF,OPS[1]->NvnS_SF,OP1,NIn,NOut,OP,d,3,Eclass);

				if (Collocated) {
					for (size_t dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 2;
					if (dim1 < 2)
						Diag[0] = 0;
					else
						Diag[dim1] = 0;
				} else {
					for (size_t dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;
					Diag[1] = 2;
				}

				sf_apply_d(&Ar_vI[NvnI*Nrc*dim1],DAr,NIn,NOut,Nrc,OP,Diag,d);
				for (size_t i = 0, iMax = NvnS*Nrc; i < iMax; i++)
					RLHS[i] += DAr[i];
			}
			free(DAr);
		} else if (Collocated && (Eclass == C_TP || Eclass == C_WEDGE) && DB.AllowSparseVOL) {
			double *const DAr = malloc(NvnS*Nrc * sizeof *DAr); // free
			for (size_t dim = 0; dim < d; dim++) {
				mm_CTN_CSR_d(NvnS,Nrc,NvnI,OPS[0]->D_Weak_sp[dim],&Ar_vI[NvnI*Nrc*dim],DAr);
				for (size_t i = 0, iMax = NvnS*Nrc; i < iMax; i++)
					RLHS[i] += DAr[i];
			}
			free(DAr);
		} else {
			for (size_t dim = 0; dim < d; dim++)
				mm_d(CBCM,CBT,CBNT,NvnS,Nrc,NvnI,1.0,1.0,OPS[0]->D_Weak[dim],&Ar_vI[NvnI*Nrc*dim],RLHS);
		}
	} else if (imex_type == 'I') {
		unsigned int const Neq  = DB.Neq,
		                   Nvar = DB.Nvar;

		double const *Ar_vI_ptr[d];
		for (size_t dim = 0; dim < d; dim++)
			Ar_vI_ptr[dim] = &Ar_vI[Nvar*Neq*NvnI*dim];

		if (Collocated) {
			if (Eclass == C_TP || Eclass == C_WEDGE) {
				struct S_OpCSR const *const *const D = OPS[0]->D_Weak_sp;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					size_t const IndAr  = (eq*Nvar+var)*NvnI,
					             IndLHS = (eq*Nvar+var)*NvnS*NvnS;

					for (size_t dim = 0; dim < d; dim++) {
						unsigned int const *rowIndex = D[dim]->rowIndex,
						                   *columns  = D[dim]->columns;
						double       const *values   = D[dim]->values;

						for (size_t i = 0; i < NvnS; i++) {
							size_t const IndD = i*NvnI;
							for (size_t Indj = *rowIndex, IndjMax = *(++rowIndex); Indj < IndjMax; Indj++) {
								RLHS[IndLHS+IndD+(*columns)] += (*values++)*Ar_vI_ptr[dim][IndAr+(*columns)];
								columns++;
							}
						}
					}
				}}
			} else {
				double const *const *const D = OPS[0]->D_Weak;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					size_t const IndAr  = (eq*Nvar+var)*NvnI,
					             IndLHS = (eq*Nvar+var)*NvnS*NvnS;

					for (size_t dim = 0; dim < d; dim++) {
						for (size_t i = 0; i < NvnS; i++) {
							size_t const IndD = i*NvnI;
							for (size_t j = 0; j < NvnI; j++)
								RLHS[IndLHS+IndD+j] += D[dim][IndD+j]*Ar_vI_ptr[dim][IndAr+j];
						}
					}
				}}
			}
		} else {
			double const *const *const D = OPS[0]->D_Weak;

			double *const DAr_vI = malloc(NvnS*NvnI * sizeof *DAr_vI); // free
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				size_t const IndAr = (eq*Nvar+var)*NvnI;

				memset(DAr_vI,0.0,NvnS*NvnI * sizeof *DAr_vI);
				for (size_t dim = 0; dim < d; dim++)
					mm_diag_d(NvnS,NvnI,&Ar_vI_ptr[dim][IndAr],D[dim],DAr_vI,1.0,1.0,'R','R');

				size_t const IndLHS = (eq*Nvar+var)*NvnS*NvnS;
				mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,1.0,DAr_vI,OPS[0]->ChiS_vI,&RLHS[IndLHS]);
			}}
			free(DAr_vI);

		/*
			// d BLAS3 calls with increased total flops
			// Ensure that finalize_LHS for VOLUME terms is modified if this is used.
			double *ChiS_vI = OPS[0]->ChiS_vI;
			double *ArChiS_vI = calloc(NvnS*NvnI*Neq*Nvar , sizeof *ArChiS_vI); // free
			memset(RLHS,0.0,NvnS*Nrc * sizeof *RLHS);
			for (size_t dim = 0; dim < d; dim++) {
				double *ArChiS_vI_ptr = ArChiS_vI;
				for (size_t i = 0; i < NvnI; i++) {
					for (size_t eq = 0; eq < Neq; eq++) {
					for (size_t var = 0; var < Nvar; var++) {
						double *ChiS_vI_ptr = &ChiS_vI[i*NvnS];
						double const *Ar_vI_ptr2 = &Ar_vI_ptr[dim][(eq*Nvar+var)*NvnI+i];
						for (size_t j = 0; j < NvnS; j++)
							*ArChiS_vI_ptr++ = (*ChiS_vI_ptr++)*(*Ar_vI_ptr2);
					}}
				}

				mm_d(CBRM,CBNT,CBNT,NvnS,NvnS*Neq*Nvar,NvnI,1.0,1.0,D[dim],ArChiS_vI,RLHS);
			}

			free(ArChiS_vI);
		*/
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Viscous_Weak(unsigned int const Nrc, double *const Ar_vI, double *const RLHS, char const imex_type,
                                  struct S_VDATA const *const VDATA)
{
	/*
	 *	Purpose:
	 *		Compute the viscous contribution to the the RHS/LHS terms for the weak formulation by applying the D_Weak
	 *		operator to the input array.
	 *
	 *	Comments:
	 *		The required operation is identical to that of finalize_VOLUME_Inviscid_Weak as the viscous flux has already
	 *		been negated.
	 */

	finalize_VOLUME_Inviscid_Weak(Nrc,Ar_vI,RLHS,imex_type,VDATA);
}

void initialize_VOLUME_LHSQ_Weak(unsigned int const Nrc, double const *const *const dFrdQ_vI, double *const *const LHSQ,
                                 struct S_VDATA const *const VDATA)
{
	/*
	 *	Purpose:
	 *		Initialize operator for VOLUME terms which have dependence on the current VOLUME as well as adjacent
	 *		VOLUMEs.
	 *
	 *	Comments:
	 *		The operator computed here is: D_Weak * dFrdQ_vI * dQ/dQhat.
	 *
	 *		Potentially ToBeModified:
	 *		The VOLUME contribution (QhatV_What) is added to VOLUME->LHS in finalize_VOLUME_LHSQV_Weak.c.
	 *		The FACE   contributions:
	 *			1) Qhat_What(LL/RR) are added VL/VR->LHS
	 *			2) Qhat_What(RL/LR) are added FACE->LHS(RL/LR)
	 *		in finalize_VOLUME_LHSQF_Weak.c.
	 */

	if (!DB.Viscous)
		return;

	unsigned int const d = DB.d;

	for (size_t dim = 0; dim < d; dim++)
		finalize_VOLUME_Inviscid_Weak(Nrc,dFrdQ_vI[dim],LHSQ[dim],'I',VDATA);
}

void finalize_VOLUME_LHSQV_Weak(struct S_VOLUME *const VOLUME)
{
	if (!DB.Viscous)
		return;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq,
	                   NvnS = VOLUME->NvnS;

	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t varQ = 0; varQ < Nvar; varQ++) {
	for (size_t var = 0; var < Nvar; var++) {
		// dQhat/dWhat is block diagonal for this term
		if (var != varQ)
			continue;

		size_t const Indev = (eq*Nvar+var)*NvnS*NvnS;
		for (size_t dim = 0; dim < d; dim++)
			mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnS,1.0,1.0,&VOLUME->LHSQ[dim][Indev],VOLUME->QhatV_What[dim],&VOLUME->LHS[Indev]);
	}}}
}

// **************************************************************************************************** //
// FACE functions
// **************************************************************************************************** //

void init_ops_FACE(struct S_OPERATORS_F *const OPS, struct S_VOLUME const *const VOLUME, struct S_FACE const *const FACE,
                   unsigned int const IndFType)
{
	/*
	 *	Comments:
	 *		For WEDGE ELEMENTs, the nOrd arrays for QUAD FACEs are stored with the TRI OPs, while those for TRI FACEs
	 *		are stored with the LINE OPs. While this is not logical, it precludes the need for an additional OP struct.
	 */

	// Initialize DB Parameters
	unsigned int const *const *const *const SF_BE = (const unsigned int *const *const *const) DB.SF_BE;

	// Standard datatypes
	unsigned int PV, PF, Vtype, Eclass, FtypeInt, IndOrdLR, IndOrdRL;

	struct S_ELEMENT *ELEMENT, *ELEMENT_SF, *ELEMENT_FACE;

	// silence
	ELEMENT_SF = NULL;

	PV       = VOLUME->P;
	PF       = FACE->P;
	Vtype    = VOLUME->type;
	Eclass   = VOLUME->Eclass;

	FtypeInt = FACE->typeInt;
	IndOrdLR = FACE->IndOrdLR;
	IndOrdRL = FACE->IndOrdRL;

	ELEMENT      = get_ELEMENT_type(Vtype);
	ELEMENT_FACE = get_ELEMENT_FACE(Vtype,IndFType);
	if ((Eclass == C_TP && SF_BE[PF][0][1]) || (Eclass == C_WEDGE && SF_BE[PF][1][1]))
		ELEMENT_SF = ELEMENT->ELEMENTclass[IndFType];
	else
		ELEMENT_SF = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[PV];
	OPS->NvnS_SF = ELEMENT_SF->NvnS[PV];
	if (FtypeInt == 's') {
		// Straight FACE Integration
		OPS->NfnI    = ELEMENT->NfnIs[PF][IndFType];
		OPS->NfnI_SF = ELEMENT_SF->NfnIs[PF][0];
		OPS->NvnI_SF = ELEMENT_SF->NvnIs[PF];

		OPS->w_fI    = ELEMENT->w_fIs[PF][IndFType];

		OPS->ChiS_fI   = (double const *const *const) ELEMENT->ChiS_fIs[PV][PF];
		OPS->I_Weak_FF = (double const *const *const) ELEMENT->Is_Weak_FF[PV][PF];

		OPS->ChiS_fI_SF   = (double const *const *const) ELEMENT_SF->ChiS_fIs[PV][PF];
		OPS->ChiS_vI_SF   = (double const *const *const) ELEMENT_SF->ChiS_vIs[PV][PF];
		OPS->I_Weak_FF_SF = (double const *const *const) ELEMENT_SF->Is_Weak_FF[PV][PF];
		OPS->I_Weak_VV_SF = (double const *const *const) ELEMENT_SF->Is_Weak_VV[PV][PF];

		OPS->ChiS_fI_sp   = (struct S_OpCSR const *const *const) ELEMENT->ChiS_fIs_sp[PV][PF];
		OPS->I_Weak_FF_sp = (struct S_OpCSR const *const *const) ELEMENT->Is_Weak_FF_sp[PV][PF];

		OPS->nOrdLR = ELEMENT_FACE->nOrd_fIs[PF][IndOrdLR];
		OPS->nOrdRL = ELEMENT_FACE->nOrd_fIs[PF][IndOrdRL];
	} else {
		// Curved FACE Integration
		OPS->NfnI    = ELEMENT->NfnIc[PF][IndFType];
		OPS->NfnI_SF = ELEMENT_SF->NfnIc[PF][0];
		OPS->NvnI_SF = ELEMENT_SF->NvnIc[PF];

		OPS->w_fI    = ELEMENT->w_fIc[PF][IndFType];

		OPS->ChiS_fI   = (double const *const *const) ELEMENT->ChiS_fIc[PV][PF];
		OPS->I_Weak_FF = (double const *const *const) ELEMENT->Ic_Weak_FF[PV][PF];

		OPS->ChiS_fI_SF   = (double const *const *const) ELEMENT_SF->ChiS_fIc[PV][PF];
		OPS->ChiS_vI_SF   = (double const *const *const) ELEMENT_SF->ChiS_vIc[PV][PF];
		OPS->I_Weak_FF_SF = (double const *const *const) ELEMENT_SF->Ic_Weak_FF[PV][PF];
		OPS->I_Weak_VV_SF = (double const *const *const) ELEMENT_SF->Ic_Weak_VV[PV][PF];

		OPS->ChiS_fI_sp   = (struct S_OpCSR const *const *const) ELEMENT->ChiS_fIc_sp[PV][PF];
		OPS->I_Weak_FF_sp = (struct S_OpCSR const *const *const) ELEMENT->Ic_Weak_FF_sp[PV][PF];

		OPS->nOrdLR = ELEMENT_FACE->nOrd_fIc[PF][IndOrdLR];
		OPS->nOrdRL = ELEMENT_FACE->nOrd_fIc[PF][IndOrdRL];
	}
}

void init_FDATA(struct S_FDATA *const FDATA, struct S_FACE const *const FACE, char const side)
{
	// Initialize DB Parameters
	unsigned int const Collocated = DB.Collocated;

	FDATA->FACE = FACE;
	FDATA->side = side;

	FDATA->P = FACE->P;

	if (side == 'L') {
		FDATA->VOLUME  = FACE->VL;
		FDATA->Vf      = FACE->VfL;
		FDATA->QhatF   = FACE->QhatL;
		FDATA->QhatF_c = FACE->QhatL_c;
	} else if (side == 'R') {
		FDATA->VOLUME  = FACE->VR;
		FDATA->Vf      = FACE->VfR;
		FDATA->QhatF   = FACE->QhatR;
		FDATA->QhatF_c = FACE->QhatR_c;
	} else {
		EXIT_UNSUPPORTED;
	}

	FDATA->f    = (FDATA->Vf)/NFREFMAX;
	FDATA->SpOp = Collocated && ((FDATA->Vf) % NFREFMAX == 0 && FDATA->VOLUME->P == FDATA->P);

	FDATA->Eclass = FDATA->VOLUME->Eclass;
	FDATA->IndFType = get_IndFType(FDATA->Eclass,FDATA->f);

	init_ops_FACE((struct S_OPERATORS_F *const) FDATA->OPS[0],FDATA->VOLUME,FACE,0);
	if (FDATA->VOLUME->type == WEDGE || FDATA->VOLUME->type == PYR)
		// Needed for sum factorized operators and alternate FACE operators (TRIs/QUADs)
		init_ops_FACE((struct S_OPERATORS_F *const) FDATA->OPS[1],FDATA->VOLUME,FACE,1);
}

void coef_to_values_fI(struct S_FDATA *const FDATA, char const coef_type, char const imex_type)
{
	/*
	 *	Purpose:
	 *		Interpolate VOLUME coefficients (of type 'coef_type') to FACE cubature nodes. In the case of the solution
	 *		gradients, also add the appropriate partial correction based on the selected viscous numerical flux.
	 *
	 *	Comments:
	 *		Various options are available for the computation:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
	 *
	 *		imex_type is only used to compute Jacobian terms for coef_type = 'Q'.
	 *
	 *		The CDG2 and BR2 fluxes are determined according to Brdar(2012, eq. (4.3)) following the comments of section
	 *		4.1 (i.e. setting eta = 0 and taking chi according to Theorem 2 part b. The area switch (eq. (4.5)) is used
	 *		such that only one of the two VOLUME contributions must be accounted for in the CDG2 flux. Note, as the
	 *		viscous fluxes are linear in the gradients, that this formulation corresponds exactly to that typically
	 *		presented for the BR2 flux (such as in eq. (10) in Bassi(2010)) with the stabilization parameter selected
	 *		according to the guidelines above. Note also that this is the analogue of the original form of the BR2 flux
	 *		(eq. (21) in Bassi(2000)) when the scaling is added.
	 *
	 *		It is currently unclear to me where the cost savings arise when using CDG2 flux as compared to the BR2 flux
	 *		as all terms must be computed for the full contribution to Qhat used in the VOLUME term. Savings were stated
	 *		as being as high as 10% in Brdar(2012). (ToBeModified)
	 *
	 *		Fidkowski(2016, p. 81) takes chi >= 2*Nf for BR2. There is also discussion of boundary condition
	 *		implementation. ToBeModified.
	 *
	 *		It is currently uncertain whether the boundary gradients should also be corrected. Currently, they are
	 *		corrected, for consistency with the internal formulation. INVESTIGATE. (ToBeModified)
	 *
	 *		For several of the viscous boundary conditions, including NoSlip_Dirichlet_T and NoSlip_Adiabatic (as they
	 *		are currently implemented), there is no dependence of boundary values on variables other than the the
	 *		variable under consideration (i.e. dWB/dWL is block diagonal). This means that only the block diagonal
	 *		entries of Qhat_What, Q_What are non-zero. This is currently not exploited. (ToBeModified)
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 *
	 *		The allowed options for coef_type are 'W' (conserved variables), 'Q' (partially corrected gradients (Qp))
	 *
	 *		To avoid confusion with (C)ofactor terms, (f)ace cubature nodes are denoted with the subscript (f)ace
	 *		(I)ntegration.
	 *
	 *	References:
	 *		Brdar(2012)-Compact_and_Stable_Discontinuous_Galerkin_Methods_for_Convection-Diffusion_Problems
	 *		Fidkowski(2016)-A_Hybridized_Discontinuous_Galerkin_Method_on_Mapped_Deforming_Domains
	 *		Bassi(2000)-A_High_Order_Discontinuous_Galerking_Method_for_Compressible_Turbulent_Flows
	 *		Bassi(2010)-Very_High-Order_Accurate_Discontinuous_Galerkin_Computation_of_Transonic_Turbulent_Flows_on_
	 *		            Aeronautical_Configurations
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	if (!(coef_type == 'W' || coef_type == 'Q'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS    = FDATA->OPS;
	struct S_VOLUME      const *const        VOLUME = FDATA->VOLUME;
	struct S_FACE              *const        FACE   = (struct S_FACE *const) FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   Vf       = FDATA->Vf,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   NvnS     = OPS[0]->NvnS;

	double **Qphat = NULL;
	if (coef_type == 'Q') {
		Qphat = malloc(d * sizeof *Qphat);

		// Add VOLUME contribution
		for (size_t dim = 0; dim < d; dim++) {
			Qphat[dim] = malloc(NvnS*Nvar * sizeof *Qphat[dim]); // free
			for (size_t i = 0; i < NvnS*Nvar; i++)
				Qphat[dim][i] = VOLUME->QhatV[dim][i];
		}

		if (DB.ViscousFluxType == FLUX_CDG2) {
			// Evaluate from which side scaling should be computed based on area switch (Brdar(2012), eq. (4.5))
			if (FDATA->side == 'L') {
				if (FACE->Boundary) {
					FACE->CDG2_side = 'L';
				} else {
					struct S_VOLUME const *const VL = FACE->VL,
					                      *const VR = FACE->VR;

					unsigned int const NvnSL = VL->NvnS,
					                   NvnSR = VR->NvnS;

					double const *const detJVL_vIL = VL->detJV_vI,
					             *const detJVR_vIR = VR->detJV_vI;

// Potentially scale by volume of reference element (necessary for mixed meshes?) (ToBeDeleted)
					double VolL = 0.0;
					for (size_t n = 0; n < NvnSL; n++)
						VolL = max(VolL,detJVL_vIL[n]);

					double VolR = 0.0;
					for (size_t n = 0; n < NvnSR; n++)
						VolR = max(VolR,detJVR_vIR[n]);

					if (VolL <= VolR)
						FACE->CDG2_side = 'L';
					else
						FACE->CDG2_side = 'R';
				}
			}
		}

		unsigned int const Boundary = FACE->Boundary;

		double chi = 0.0;
		if (DB.ViscousFluxType == FLUX_CDG2) {
			if (FDATA->side == FACE->CDG2_side) {
				// Multiplied by 2.0 as the viscous flux terms from each side are subsequently averaged.
				if (Boundary)
					chi =     (DB.NfMax);
				else
					chi = 2.0*(DB.NfMax);
			}
		} else if (DB.ViscousFluxType == FLUX_BR2) {
			if (Boundary)
				chi =     (DB.NfMax);
			else
				chi =     (DB.NfMax);
		} else {
			EXIT_UNSUPPORTED;
		}

		// Add partial FACE contribution
		if (chi != 0.0) {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0; i < NvnS*Nvar; i++)
					Qphat[dim][i] += chi*FDATA->QhatF[dim][i];
			}
		}

		// Compute linearized contributions of partially corrected Q (Qp_What**)
		if (imex_type == 'I') {
			struct S_VOLUME const *const VL = FACE->VL,
			                      *const VR = FACE->VR;

			unsigned int const NvnSL = VL->NvnS,
			                   NvnSR = VR->NvnS;

			// Qp_What** arrays are freed in finalize_implicit_FACE_Q_Weak.
			if (Boundary) {
				double *const Qphat_What = malloc(NvnSL*NvnSL * sizeof *Qphat_What); // free

				double **Qp_What = malloc(d * sizeof *Qp_What); // keep
				for (size_t dim = 0; dim < d; dim++) {
					Qp_What[dim] = malloc(NfnI*NvnSL*Nvar*Nvar * sizeof *Qp_What[dim]); // keep

					for (size_t varQ = 0; varQ < Nvar; varQ++) {
					for (size_t var = 0; var < Nvar; var++) {
						size_t const Indvv = varQ*Nvar+var;
						if (varQ == var) {
							for (size_t i = 0; i < NvnSL*NvnSL; i++)
								Qphat_What[i] = VL->QhatV_What[dim][i];
						} else {
							memset(Qphat_What,0.0,NvnSL*NvnSL * sizeof *Qphat_What);
						}

						if (chi != 0.0) {
							for (size_t i = 0; i < NvnSL*NvnSL; i++)
								Qphat_What[i] += chi*(FACE->QhatL_WhatL[dim][Indvv*NvnSL*NvnSL+i]);
						}
						mm_d(CBRM,CBNT,CBNT,NfnI,NvnSL,NvnSL,1.0,0.0,OPS[0]->ChiS_fI[Vf],Qphat_What,&Qp_What[dim][Indvv*NfnI*NvnSL]);
					}}
				}
				free(Qphat_What);
				FDATA->Qp_WhatL = Qp_What;
			} else {
				if (FDATA->side == 'L') {
					double *const QphatL_WhatL = malloc(NvnSL*NvnSL * sizeof *QphatL_WhatL); // free

					double **Qp_WhatL = malloc(d * sizeof *Qp_WhatL), // keep
					       **Qp_WhatR = malloc(d * sizeof *Qp_WhatR); // keep
					for (size_t dim = 0; dim < d; dim++) {
						Qp_WhatL[dim] = malloc(NfnI*NvnSL * sizeof *Qp_WhatL[dim]); // keep
						for (size_t i = 0; i < NvnSL*NvnSL; i++)
							QphatL_WhatL[i] = VL->QhatV_What[dim][i];
						if (chi != 0.0) {
							for (size_t i = 0; i < NvnSL*NvnSL; i++)
								QphatL_WhatL[i] += chi*(FACE->QhatL_WhatL[dim][i]);
						}
						mm_d(CBRM,CBNT,CBNT,NfnI,NvnSL,NvnSL,1.0,0.0,OPS[0]->ChiS_fI[Vf],QphatL_WhatL,Qp_WhatL[dim]);

						Qp_WhatR[dim] = calloc(NfnI*NvnSR , sizeof *Qp_WhatR[dim]); // keep
						if (chi != 0.0)
							mm_d(CBRM,CBNT,CBNT,NfnI,NvnSR,NvnSL,chi,0.0,OPS[0]->ChiS_fI[Vf],FACE->QhatL_WhatR[dim],Qp_WhatR[dim]);
					}
					free(QphatL_WhatL);

					FDATA->Qp_WhatL = Qp_WhatL;
					FDATA->Qp_WhatR = Qp_WhatR;
				} else if (FDATA->side == 'R') {
					double *const QphatR_WhatR = malloc(NvnSR*NvnSR * sizeof *QphatR_WhatR); // free

					double **Qp_WhatL = malloc(d * sizeof *Qp_WhatL), // keep
					       **Qp_WhatR = malloc(d * sizeof *Qp_WhatR); // keep
					for (size_t dim = 0; dim < d; dim++) {
						Qp_WhatL[dim] = calloc(NfnI*NvnSL , sizeof *Qp_WhatL[dim]); // keep
						if (chi != 0.0)
							mm_d(CBRM,CBNT,CBNT,NfnI,NvnSL,NvnSR,chi,0.0,OPS[0]->ChiS_fI[Vf],FACE->QhatR_WhatL[dim],Qp_WhatL[dim]);

						Qp_WhatR[dim] = malloc(NfnI*NvnSR * sizeof *Qp_WhatR[dim]); // keep
						for (size_t i = 0; i < NvnSR*NvnSR; i++)
							QphatR_WhatR[i] = VR->QhatV_What[dim][i];
						if (chi != 0.0) {
							for (size_t i = 0; i < NvnSR*NvnSR; i++)
								QphatR_WhatR[i] += chi*(FACE->QhatR_WhatR[dim][i]);
						}
						mm_d(CBRM,CBNT,CBNT,NfnI,NvnSR,NvnSR,1.0,0.0,OPS[0]->ChiS_fI[Vf],QphatR_WhatR,Qp_WhatR[dim]);
					}
					free(QphatR_WhatR);

					FDATA->Qp_WhatL = Qp_WhatL;
					FDATA->Qp_WhatR = Qp_WhatR;
				}
			}
		}
	}

	unsigned int const P      = FDATA->P,
	                   Eclass = FDATA->Eclass,
	                   f      = FDATA->f,
	                   SpOp   = FDATA->SpOp,
	                   *const VFPartUnity         = DB.VFPartUnity,
	                   *const *const *const SF_BE = (const unsigned int *const *const *const) DB.SF_BE;

	unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX], NOut0, NOut1;
	double const *OP[DMAX], *const *OP0, *const *OP1;

	if (Eclass == C_TP && SF_BE[P][0][1]) {
		get_sf_parametersF(OPS[0]->NvnS_SF,OPS[0]->NvnI_SF,OPS[0]->ChiS_vI_SF,
		                   OPS[0]->NvnS_SF,OPS[0]->NfnI_SF,OPS[0]->ChiS_fI_SF,NIn,NOut,OP,d,Vf,C_TP);

		if (SpOp) {
			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 2;
			Diag[f/2] = 0;
		} else {
			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 0;
		}

		if (coef_type == 'W') {
			sf_apply_d(VOLUME->What,FDATA->W_fIL,NIn,NOut,Nvar,OP,Diag,d);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				sf_apply_d(Qphat[dim],FDATA->Qp_fIL[dim],NIn,NOut,Nvar,OP,Diag,d);
		}
	} else if (Eclass == C_WEDGE && SF_BE[P][1][1]) {
		if (f < 3) { OP0   = OPS[0]->ChiS_fI_SF, OP1   = OPS[1]->ChiS_vI_SF;
		             NOut0 = OPS[0]->NfnI_SF,    NOut1 = OPS[1]->NvnI_SF;
		} else {     OP0   = OPS[0]->ChiS_vI_SF, OP1   = OPS[1]->ChiS_fI_SF;
		             NOut0 = OPS[0]->NvnI_SF,    NOut1 = OPS[1]->NfnI_SF; }
		get_sf_parametersF(OPS[0]->NvnS_SF,NOut0,OP0,OPS[1]->NvnS_SF,NOut1,OP1,NIn,NOut,OP,d,Vf,C_WEDGE);

		if (SpOp) {
			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 2;
			if (f < 3)
				Diag[0] = 0;
			else
				Diag[2] = 0;
		} else {
			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 0;
			Diag[1] = 2;
		}

		if (coef_type == 'W') {
			sf_apply_d(VOLUME->What,FDATA->W_fIL,NIn,NOut,Nvar,OP,Diag,d);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				sf_apply_d(Qphat[dim],FDATA->Qp_fIL[dim],NIn,NOut,Nvar,OP,Diag,d);
		}
	} else if (((SpOp && (Eclass == C_TP || Eclass == C_WEDGE)) || (VFPartUnity[Eclass])) && DB.AllowSparseFACE) {
		if (coef_type == 'W') {
			mm_CTN_CSR_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI_sp[Vf],VOLUME->What,FDATA->W_fIL);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				mm_CTN_CSR_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI_sp[Vf],Qphat[dim],FDATA->Qp_fIL[dim]);
		}
	} else {
		if (coef_type == 'W') {
			mm_CTN_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI[Vf],VOLUME->What,FDATA->W_fIL);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				mm_CTN_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI[Vf],Qphat[dim],FDATA->Qp_fIL[dim]);
		}
	}

	if (coef_type == 'Q')
		array_free2_d(d,Qphat);
}

static void evaluate_boundary(struct S_BC *const BCdata, bool const ComputeGradient)
{
	unsigned int const BC = BCdata->BC;

	if (ComputeGradient) {
		if (!(BC % BC_STEP_SC == BC_NOSLIP_T         ||
		      BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC ||
		      BC % BC_STEP_SC == BC_DIRICHLET        ||
		      BC % BC_STEP_SC == BC_NEUMANN)) {
			EXIT_UNSUPPORTED;
		}
	}

	compute_boundary_values(BCdata);
}

static void evaluate_jacobian_boundary(struct S_BC *const BCdata, bool const ComputeGradient)
{
	if (ComputeGradient)
		EXIT_UNSUPPORTED;

	compute_jacobian_boundary_values(BCdata);
}

void compute_WR_fIL(struct S_FDATA const *const FDATA, double const *const WL_fIL, double *const WR_fIL)
{
	/*
	 *	Purpose:
	 *		Compute W of the (R)ight VOLUME at the FACE cubature nodes corresponding to the (L)eft VOLUME.
	 *
	 *	Comments:
	 *		For internal and periodic BCs this is done by interpolating What from the right VOLUME to the FACE cubature
	 *		nodes and then rearranging the result such that the ordering corresponds to that seen from the left VOLUME
	 *		(recall that nOrdLR gives the (n)ode (Ord)ering from (L)eft to (R)ight).
	 *		For other boundary conditions, the solution is computed directly with the correct ordering using the
	 *		appropriate boundary condition functions.
	 */

	struct S_OPERATORS_F const *const *const OPS  = FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const Nvar          = DB.Nvar,
	                   IndFType      = FDATA->IndFType,
	                   BC            = FACE->BC,
	                   NfnI          = OPS[IndFType]->NfnI,
	                   *const nOrdRL = OPS[IndFType]->nOrdRL;

	double const *const XYZ_fIL = FACE->XYZ_fI,
	             *const n_fIL   = FACE->n_fI;

	if (!FACE->Boundary) { // Internal/Periodic FACE
		coef_to_values_fI((struct S_FDATA *const) FDATA,'W','E');
		array_rearrange_d(NfnI,Nvar,nOrdRL,'C',WR_fIL);
	} else { // Boundary FACE
		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;
		BCdata->BC  = BC;

		BCdata->XYZ = XYZ_fIL;
		BCdata->nL  = n_fIL;
		BCdata->WL  = WL_fIL;
		BCdata->WB  = WR_fIL;
		BCdata->QL  = NULL;
		BCdata->QB  = NULL;

		evaluate_boundary(BCdata,0);
		free(BCdata);
	}
}

void compute_WR_QpR_fIL(struct S_FDATA const *const FDATA, double const *const WL_fIL, double *const WR_fIL,
                           double const *const *const QpL_fIL, double *const *const QpR_fIL, char const imex_type)
{
	/*
	 *	Purpose:
	 *		Compute W and Qp of the (R)ight VOLUME at the FACE cubature nodes corresponding to the (L)eft VOLUME.
	 */

	struct S_OPERATORS_F const *const *const OPS  = FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   BC       = FACE->BC,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   *const nOrdRL = OPS[IndFType]->nOrdRL;

	if (!FACE->Boundary) { // Internal/Periodic FACE
		compute_WR_fIL(FDATA,WL_fIL,WR_fIL);

		coef_to_values_fI((struct S_FDATA *const) FDATA,'Q',imex_type);
		for (size_t dim = 0; dim < d; dim++)
			array_rearrange_d(NfnI,Nvar,nOrdRL,'C',QpR_fIL[dim]);
	} else { // Boundary FACE
		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;
		BCdata->BC  = BC;

		BCdata->XYZ = FACE->XYZ_fI;
		BCdata->nL  = FACE->n_fI;
		BCdata->WL  = WL_fIL;
		BCdata->WB  = WR_fIL;
		BCdata->QL  = QpL_fIL;
		BCdata->QB  = QpR_fIL;

		evaluate_boundary(BCdata,1);
		free(BCdata);
	}
}

void compute_numerical_flux(struct S_FDATA const *const FDATA, char const imex_type)
{
	/*
	 *	Purpose:
	 *		Compute the numerical flux (and its Jacobians wrt to WL and WR if applicable) evaluated at the FACE cubature
	 *		nodes as seen from the left VOLUME.
	 *
	 *	Comments:
	 *		If on a boundary FACE, the boundary Jacobian contributions are included in dnFluxNumdWL_fIL.
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = d+2,
	                   Neq      = d+2,
	                   IndFType = FDATA->IndFType,
	                   Boundary = FACE->Boundary,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL            = FACE->n_fI,
	             *const WL_fIL           = FDATA->NFluxData->WL_fIL,
	             *const WR_fIL           = FDATA->NFluxData->WR_fIL;
	double       *const nFluxNum_fIL     = FDATA->NFluxData->nFluxNum_fI,
	             *const dnFluxNumdWL_fIL = FDATA->NFluxData->dnFluxNumdWL_fI,
	             *const dnFluxNumdWR_fIL = FDATA->NFluxData->dnFluxNumdWR_fI;

	// Consider combining numerical flux and jacobian computations (ToBeDeleted)
	// Consder passing struct to numerical flux functions (as for standard fluxes) (ToBeDeleted)
	switch (DB.InviscidFluxType) {
	case FLUX_LF:
		flux_LF(NfnI,1,WL_fIL,WR_fIL,nFluxNum_fIL,n_fIL,d,Neq);
		if (imex_type == 'I') {
			jacobian_flux_LF(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWL_fIL,n_fIL,d,Neq,'L');
			jacobian_flux_LF(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWR_fIL,n_fIL,d,Neq,'R');
		}
		break;
	case FLUX_ROE:
		flux_Roe(NfnI,1,WL_fIL,WR_fIL,nFluxNum_fIL,n_fIL,d,Neq);
		if (imex_type == 'I') {
			jacobian_flux_Roe(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWL_fIL,n_fIL,d,Neq,'L');
			jacobian_flux_Roe(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWR_fIL,n_fIL,d,Neq,'R');
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	// Include the BC information in dnFluxNumWL_fIL if on a boundary
	if (imex_type == 'I' && Boundary) {
		unsigned int const BC = FACE->BC;

		double const *const XYZ_fIL    = FACE->XYZ_fI;
		double       *const dWRdWL_fIL = malloc(NfnI*Nvar*Nvar * sizeof *dWRdWL_fIL); // free

		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;
		BCdata->BC  = BC;

		BCdata->XYZ = XYZ_fIL;
		BCdata->nL  = n_fIL;
		BCdata->WL  = WL_fIL;
		BCdata->QL  = NULL;

		BCdata->dWBdWL = dWRdWL_fIL;

		evaluate_jacobian_boundary(BCdata,0);
		free(BCdata);

		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			size_t const InddnFdWL = (eq*Nvar+var)*NfnI;

			for (size_t i = 0; i < Nvar; i++) {
				size_t const InddnFdWR = (eq*Neq+i)*NfnI,
				             InddWRdWL = (var*Nvar+i)*NfnI;
				for (size_t n = 0; n < NfnI; n++)
					dnFluxNumdWL_fIL[InddnFdWL+n] += dnFluxNumdWR_fIL[InddnFdWR+n]*dWRdWL_fIL[InddWRdWL+n];
			}
		}}

		free(dWRdWL_fIL);
	}
}

void compute_numerical_solution(struct S_FDATA const *const FDATA, char const imex_type)
{
	/*
	 *	Purpose:
	 *		Compute the numerical solution (and its Jacobians wrt to WL and WR if applicable) evaluated at the FACE
	 *		cubature nodes as seen from the left VOLUME (used for the FACE contribution to the weak gradient).
	 *
	 *	Comments:
	 *		It is currently hard-coded that a central numerical solution is used. As the numerical solution is linear in
	 *		the solution variables, its Jacobian is constant for all variables/equations. However, the potential
	 *		dependence of the boundary conditions on all variables results in different linearized terms being required
	 *		for each variable of each equation. This is why the dnSolNumdW(L/R) terms are potentially redundantly
	 *		computed.
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   Neq      = DB.Neq,
	                   IndFType = FDATA->IndFType,
	                   Boundary = FACE->Boundary,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL  = FACE->n_fI,
	             *const WL_fIL = FDATA->NFluxData->WL_fIL,
	             *const WR_fIL = FDATA->NFluxData->WR_fIL;

	double       *const *const nSolNum_fIL     = FDATA->NFluxData->nSolNum_fI,
	             *const *const dnSolNumdWL_fIL = FDATA->NFluxData->dnSolNumdWL_fI,
	             *const *const dnSolNumdWR_fIL = FDATA->NFluxData->dnSolNumdWR_fI;

	for (size_t dim = 0; dim < d; dim++) {
	for (size_t var = 0; var < Nvar; var++) {
	for (size_t n = 0; n < NfnI; n++) {
		nSolNum_fIL[dim][var*NfnI+n] = n_fIL[n*d+dim]*0.5*(WL_fIL[var*NfnI+n]+WR_fIL[var*NfnI+n]);
	}}}

	if (imex_type == 'I') {
		unsigned int eqMax, varMax;
		if (Boundary) {
			eqMax  = Neq;
			varMax = Nvar;
		} else {
			eqMax  = 1;
			varMax = 1;
		}

		for (size_t dim = 0; dim < d; dim++) {
		for (size_t eq = 0; eq < eqMax; eq++) {
		for (size_t var = 0; var < varMax; var++) {
			size_t const Indeqvar = (eq*Nvar+var)*NfnI;

			if (eq != var) {
				memset(&dnSolNumdWL_fIL[dim][Indeqvar],0.0,NfnI * sizeof dnSolNumdWL_fIL[dim][Indeqvar]);
				memset(&dnSolNumdWR_fIL[dim][Indeqvar],0.0,NfnI * sizeof dnSolNumdWR_fIL[dim][Indeqvar]);
				continue;
			}

			for (size_t n = 0; n < NfnI; n++) {
				dnSolNumdWL_fIL[dim][Indeqvar+n] = n_fIL[n*d+dim]*0.5;
				dnSolNumdWR_fIL[dim][Indeqvar+n] = n_fIL[n*d+dim]*0.5;
			}
		}}}
	}

	// Include the BC information in dnSolNumWL_fIL if on a boundary (Only include contribution from W)
	if (imex_type == 'I' && Boundary) {
		unsigned int const BC = FACE->BC;

		double const *const XYZ_fIL    = FACE->XYZ_fI;
		double       *const dWRdWL_fIL = malloc(NfnI*Nvar*Nvar * sizeof *dWRdWL_fIL); // free

		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;
		BCdata->BC  = BC;

		BCdata->XYZ = XYZ_fIL;
		BCdata->nL  = n_fIL;
		BCdata->WL  = WL_fIL;
		BCdata->WB  = (double *const) WR_fIL;
		BCdata->QL  = NULL;
		BCdata->QB  = NULL;

		BCdata->dWBdWL = dWRdWL_fIL;
		evaluate_jacobian_boundary(BCdata,0);
		free(BCdata);

		for (size_t dim = 0; dim < d; dim++) {
		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			size_t const InddnSdWL = (eq*Nvar+var)*NfnI;

			for (size_t i = 0; i < Nvar; i++) {
				size_t const InddnSdWR = (eq*Neq+i)*NfnI,
				             InddWRdWL = (var*Nvar+i)*NfnI;
				for (size_t n = 0; n < NfnI; n++)
					dnSolNumdWL_fIL[dim][InddnSdWL+n] += dnSolNumdWR_fIL[dim][InddnSdWR+n]*dWRdWL_fIL[InddWRdWL+n];
			}
		}}}

		free(dWRdWL_fIL);
	}
}

static void correct_numerical_solution_strong(struct S_FDATA const *const FDATA, char const imex_type, char const side)
{
	/*
	 *	Purpose:
	 *		Add internal contribution to nSolNum_fIL when the strong form is used for the first equation of the mixed
	 *		form.
	 *
	 *	Comments:
	 *		Scale by Jacobian term here as well as this was added to nSolNum_fIL for the Left VOLUME contribution.
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL     = FACE->n_fI,
	             *const detJF_fIL = FACE->detJF_fI,
	             *const WL_fIL    = FDATA->NFluxData->WL_fIL,
	             *const WR_fIL    = FDATA->NFluxData->WR_fIL;

	double       *const *const nSolNum_fIL     = FDATA->NFluxData->nSolNum_fI,
	             *const *const dnSolNumdWL_fIL = FDATA->NFluxData->dnSolNumdWL_fI,
	             *const *const dnSolNumdWR_fIL = FDATA->NFluxData->dnSolNumdWR_fI;

	if (side == 'L') {
		if (imex_type == 'E') {
			for (size_t dim = 0; dim < d; dim++) {
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NfnI; n++) {
				nSolNum_fIL[dim][var*NfnI+n] -= n_fIL[n*d+dim]*detJF_fIL[n]*WL_fIL[var*NfnI+n];
			}}}
		} else if (imex_type == 'I') {
			unsigned int eqMax, varMax;
			if (FACE->Boundary) {
				eqMax  = DB.Neq;
				varMax = DB.Nvar;
			} else {
				eqMax  = 1;
				varMax = 1;
			}

			// Only need to correct dnSolNumdWL_fIL
			for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < eqMax; eq++) {
			for (size_t var = 0; var < varMax; var++) {
				if (eq != var)
					continue;

				size_t const Indeqvar = (eq*Nvar+var)*NfnI;
				for (size_t n = 0; n < NfnI; n++) {
					dnSolNumdWL_fIL[dim][Indeqvar+n] -= n_fIL[n*d+dim]*detJF_fIL[n];
				}
			}}}
		}
	} else if (side == 'R') {
		if (FACE->Boundary)
			EXIT_UNSUPPORTED;

		if (imex_type == 'E') {
			for (size_t dim = 0; dim < d; dim++) {
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NfnI; n++) {
				nSolNum_fIL[dim][var*NfnI+n] += n_fIL[n*d+dim]*detJF_fIL[n]*(WL_fIL[var*NfnI+n]-WR_fIL[var*NfnI+n]);
			}}}
		} else if (imex_type == 'I') {
			unsigned int const eqMax  = 1,
			                   varMax = 1;

			for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < eqMax; eq++) {
			for (size_t var = 0; var < varMax; var++) {
				size_t const Indeqvar = (eq*Nvar+var)*NfnI;
				for (size_t n = 0; n < NfnI; n++) {
					dnSolNumdWL_fIL[dim][Indeqvar+n] += n_fIL[n*d+dim]*detJF_fIL[n];
					dnSolNumdWR_fIL[dim][Indeqvar+n] -= n_fIL[n*d+dim]*detJF_fIL[n];
				}
			}}}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void dot_with_normal(unsigned int const Nn, unsigned int const NCol, double const *const n_fIL,
                            double const *const ANum_fIL, double *const nANum_fIL)
{
	unsigned int const d = DB.d;

	memset(nANum_fIL,0.0,Nn*NCol * sizeof *nANum_fIL);
	for (size_t col = 0; col < NCol; col++) {
		double const *const ANum_ptr = &ANum_fIL[Nn*d*col];
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t n = 0; n < Nn; n++) {
			nANum_fIL[n+Nn*col] += n_fIL[n*d+dim]*ANum_ptr[Nn*dim+n];
		}}
	}
}

void compute_numerical_flux_viscous(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                    char const imex_type)
{
	/*
	 *	Purpose:
	 *		Compute the numerical viscous flux (and its Jacobians wrt to WL and WR if applicable) evaluated at the FACE
	 *		cubature nodes as seen from the left VOLUME.
	 *
	 *	Comments:
	 *		If on a boundary FACE, the boundary Jacobian contributions are included in dnFluxViscNumdWL_fIL for the
	 *		solution only. The boundary Jacobian contributions with respect to Q are treated here directly.
	 *
	 *		The CDG2 and BR2 fluxes are determined according to Brdar(2012, eq. (4.3)) following the comments of section
	 *		4.1 (i.e. setting eta = 0 and taking chi according to Theorem 2 part b. The area switch (eq. (4.5)) is used
	 *		such that only one of the two VOLUME contributions must be accounted for in the CDG2 flux. Note, as the
	 *		viscous fluxes are linear in the gradients, that this formulation corresponds exactly to that typically
	 *		presented for the BR2 flux (such as in eq. (10) in Bassi(2010)) with the stabilization parameter selected
	 *		according to the guidelines above. Note also that this is the analogue of the original form of the BR2 flux
	 *		(eq. (21) in Bassi(2000)) when the scaling is added.
	 *
	 *		It is currently unclear to me where the cost savings arise when using CDG2 flux as compared to the BR2 flux
	 *		as all terms must be computed for the full contribution to Qhat used in the VOLUME term. Savings were stated
	 *		as being as high as 10% in Brdar(2012). (ToBeModified)
	 *
	 *		Fidkowski(2016, p. 81) takes chi >= 2*Nf for BR2. There is also discussion of boundary condition
	 *		implementation. ToBeModified.
	 *
	 *
	 *		In the treatment for the linearization with respect to the boundary conditions, it was assumed that:
	 *			1) dQR/dQL = I (Identity);
	 *			2) dQR/dWL = 0 (Zero).
	 *
	 *		When this is true the Q contribution to the full Jacobian can be simplified as follows:
	 *			dnFNumVisc/dWhatL = df/dQ*(dQ/dQL*dQL/dWhatL + dQ/dQR*(dQR/dWL*dWL/dWhatL + dQR/dQL*dQL/dWhatL))
	 *			                  = df/dQ*( 0.5  * QL_WhatL  +  0.5  *(  0.0  *  ChiS_vI  +    I   * QL_WhatL ))
	 *			                  = df/dQ*QL_WhatL
	 *		Note that only the dfdQ contribution (with the normal) is computed here.
	 *
	 *		For the Poisson Neumann boundary condition, dQR/dQL = -I resulting in dnFNumVisc/dWhatL = 0.
	 *
	 *		If this is not the case for this boundary condition, update the boundary treatment appropriately.
	 *
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 *
	 *	References:
	 *		Brdar(2012)-Compact_and_Stable_Discontinuous_Galerkin_Methods_for_Convection-Diffusion_Problems
	 *		Fidkowski(2016)-A_Hybridized_Discontinuous_Galerkin_Method_on_Mapped_Deforming_Domains
	 *		Bassi(2000)-A_High_Order_Discontinuous_Galerking_Method_for_Compressible_Turbulent_Flows
	 *		Bassi(2010)-Very_High-Order_Accurate_Discontinuous_Galerkin_Computation_of_Transonic_Turbulent_Flows_on_
	 *		            Aeronautical_Configurations
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPSL = (struct S_OPERATORS_F const *const *const) FDATAL->OPS;
	struct S_FACE        const *const        FACE = FDATAL->FACE;

	char const *PDE = DB.PDE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   Neq      = DB.Neq,
	                   IndFType = FDATAL->IndFType,
	                   Boundary = FACE->Boundary,
	                   NfnI     = OPSL[IndFType]->NfnI;

	double const *const n_fIL                       = FACE->n_fI,
	             *const WL_fIL                      = FDATAL->NFluxData->WL_fIL,
	             *const WR_fIL                      = FDATAL->NFluxData->WR_fIL,
	             *const *const QpL_fIL              = (double const *const *const) FDATAL->Qp_fIL,
	             *const *const QpR_fIL              = (double const *const *const) FDATAR->Qp_fIL;
	double       *const nFluxViscNum_fIL            = FDATAL->NFluxData->nFluxViscNum_fI,
	             *const dnFluxViscNumdWL_fIL        = FDATAL->NFluxData->dnFluxViscNumdWL_fI,
	             *const dnFluxViscNumdWR_fIL        = FDATAL->NFluxData->dnFluxViscNumdWR_fI,
	             *const *const dnFluxViscNumdQL_fIL = (double *const *const) FDATAL->NFluxData->dnFluxViscNumdQL_fI,
	             *const *const dnFluxViscNumdQR_fIL = (double *const *const) FDATAL->NFluxData->dnFluxViscNumdQR_fI;

	double *const FluxViscNum_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscNum_fIL); // free

	double *dFluxViscNumdWL_fIL = NULL, **dFluxViscNumdQL_fIL = NULL,
	       *dFluxViscNumdWR_fIL = NULL, **dFluxViscNumdQR_fIL = NULL;
	if (imex_type == 'I') {
		dFluxViscNumdWL_fIL  = malloc(NfnI*d*Neq*Nvar * sizeof *dFluxViscNumdWL_fIL);  // free
		dFluxViscNumdWR_fIL  = malloc(NfnI*d*Neq*Nvar * sizeof *dFluxViscNumdWR_fIL);  // free
		dFluxViscNumdQL_fIL  = malloc(d               * sizeof *dFluxViscNumdQL_fIL);  // free
		dFluxViscNumdQR_fIL  = malloc(d               * sizeof *dFluxViscNumdQR_fIL);  // free
		for (size_t dim = 0; dim < d; dim++) {
			dFluxViscNumdQL_fIL[dim]  = malloc(NfnI*d*Neq*Nvar * sizeof *dFluxViscNumdQL_fIL[dim]);  // free
			dFluxViscNumdQR_fIL[dim]  = malloc(NfnI*d*Neq*Nvar * sizeof *dFluxViscNumdQR_fIL[dim]);  // free
		}
	}

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = NfnI;
	FLUXDATA->Nel = 1;

	if (Boundary) {
		double *const W_fIL = malloc(NfnI*Nvar * sizeof *W_fIL); // free
		for (size_t i = 0; i < NfnI*Nvar; i++)
			W_fIL[i] = 0.5*(WL_fIL[i]+WR_fIL[i]);

		double **Qp_fIL = malloc(d * sizeof *Qp_fIL); // free
		for (size_t dim = 0; dim < d; dim++) {
			Qp_fIL[dim] = malloc(NfnI*Nvar * sizeof *Qp_fIL[dim]); // free
			for (size_t i = 0; i < NfnI*Nvar; i++)
				Qp_fIL[dim][i] = 0.5*(QpL_fIL[dim][i]+QpR_fIL[dim][i]);
		}

		if (strstr(PDE,"NavierStokes")) {
			// Note: This term can be computed in jacobian_flux_viscous for imex_type == 'I' (ToBeModified)
			FLUXDATA->W = W_fIL;
			FLUXDATA->Q = (double const *const *const) Qp_fIL;
			FLUXDATA->F = FluxViscNum_fIL;
			flux_viscous(FLUXDATA);
		} else if (strstr(PDE,"Poisson")) {
			if (imex_type != 'I')
				EXIT_UNSUPPORTED;

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0, iMax = NfnI*Nvar; i < iMax; i++)
					FluxViscNum_fIL[iMax*dim+i] = Qp_fIL[dim][i];
			}
		} else {
			EXIT_UNSUPPORTED;
		}

		if (imex_type == 'I') {
			if (strstr(PDE,"NavierStokes")) {
				FLUXDATA->F    = NULL;
				FLUXDATA->dFdW = dFluxViscNumdWL_fIL;
				FLUXDATA->dFdQ = dFluxViscNumdQL_fIL;
				jacobian_flux_viscous(FLUXDATA);

				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdWL_fIL,dnFluxViscNumdWL_fIL);
				for (size_t i = 0; i < NfnI*Neq*Nvar; i++) {
					dnFluxViscNumdWL_fIL[i] *= 0.5;
					dnFluxViscNumdWR_fIL[i]  = dnFluxViscNumdWL_fIL[i];
				}

				/* Commented as this is identical to doing nothing to dnFluxViscNumQL_fIL. See comments.
				for (size_t dim = 0; dim < d; dim++) {
					for (size_t i = 0; i < NfnI*Neq*Nvar; i++) {
						// dQ/dQL
						dnFluxViscNumdQL_fIL[dim][i] *= 0.5;

						// dQ/dQR
						dnFluxViscNumdQR_fIL[dim][i]  = dnFluxViscNumdQL_fIL[dim][i];

						// Add effect of boundary Jacobian to dnFluxViscNumQL_fIL
						dnFluxViscNumdQL_fIL[dim][i] += 1.0*dnFluxViscNumdQR_fIL[dim][i];
					}
				} */
				for (size_t dim = 0; dim < d; dim++)
					dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdQL_fIL[dim],dnFluxViscNumdQL_fIL[dim]);
			} else if (strstr(PDE,"Poisson")) {
				unsigned int const BC = FACE->BC;
		      	if (BC % BC_STEP_SC == BC_DIRICHLET) {
					for (size_t dim = 0; dim < d; dim++) {
						memset(dFluxViscNumdQL_fIL[dim],0.0,NfnI*d*Neq*Nvar * sizeof *dFluxViscNumdQL_fIL[dim]);
						for (size_t i = 0, iMax = NfnI*Nvar; i < iMax; i++)
							dFluxViscNumdQL_fIL[dim][iMax*dim+i] = 1.0;
					}

					for (size_t dim = 0; dim < d; dim++)
						dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdQL_fIL[dim],dnFluxViscNumdQL_fIL[dim]);
				} else if (BC % BC_STEP_SC == BC_NEUMANN) {
					for (size_t dim = 0; dim < d; dim++)
						memset(dnFluxViscNumdQL_fIL[dim],0.0,NfnI*Neq*Nvar * sizeof *dnFluxViscNumdQL_fIL[dim]);
				} else {
					EXIT_UNSUPPORTED;
				}
			} else {
				EXIT_UNSUPPORTED;
			}

		}
		free(W_fIL);
		array_free2_d(d,Qp_fIL);
	} else {
		if (strstr(PDE,"NavierStokes")) {
			double *const FluxViscL_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscL_fIL); // free
			double *const FluxViscR_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscR_fIL); // free

			FLUXDATA->W = WL_fIL;
			FLUXDATA->Q = QpL_fIL;
			FLUXDATA->F = FluxViscL_fIL;
			// Note: This term can be computed in jacobian_flux_viscous for imex_type == 'I' (ToBeModified)
			flux_viscous(FLUXDATA);

			FLUXDATA->W = WR_fIL;
			FLUXDATA->Q = QpR_fIL;
			FLUXDATA->F = FluxViscR_fIL;
			// Note: This term can be computed in jacobian_flux_viscous for imex_type == 'I' (ToBeModified)
			flux_viscous(FLUXDATA);

			for (size_t i = 0; i < NfnI*d*Neq; i++)
				FluxViscNum_fIL[i] = 0.5*(FluxViscL_fIL[i]+FluxViscR_fIL[i]);

			free(FluxViscL_fIL);
			free(FluxViscR_fIL);
		} else if (strstr(PDE,"Poisson")) {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0, iMax = NfnI*Nvar; i < iMax; i++)
					FluxViscNum_fIL[iMax*dim+i] = 0.5*(QpL_fIL[dim][i] + QpR_fIL[dim][i]);
			}
		} else {
			EXIT_UNSUPPORTED;
		}

		if (imex_type == 'I') {
			if (strstr(PDE,"NavierStokes")) {
				FLUXDATA->F = NULL;

				FLUXDATA->W    = WL_fIL;
				FLUXDATA->Q    = QpL_fIL;
				FLUXDATA->dFdW = dFluxViscNumdWL_fIL;
				FLUXDATA->dFdQ = dFluxViscNumdQL_fIL;
				jacobian_flux_viscous(FLUXDATA);

				FLUXDATA->W    = WR_fIL;
				FLUXDATA->Q    = QpR_fIL;
				FLUXDATA->dFdW = dFluxViscNumdWR_fIL;
				FLUXDATA->dFdQ = dFluxViscNumdQR_fIL;
				jacobian_flux_viscous(FLUXDATA);

				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdWL_fIL,dnFluxViscNumdWL_fIL);
				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdWR_fIL,dnFluxViscNumdWR_fIL);

				for (size_t i = 0; i < NfnI*Neq*Nvar; i++) {
					dnFluxViscNumdWL_fIL[i] *= 0.5;
					dnFluxViscNumdWR_fIL[i] *= 0.5;
				}
			} else if (strstr(PDE,"Poisson")) {
				for (size_t dim = 0; dim < d; dim++) {
					memset(dFluxViscNumdQL_fIL[dim],0.0,NfnI*d*Neq*Nvar * sizeof *dFluxViscNumdQL_fIL[dim]);
					memset(dFluxViscNumdQR_fIL[dim],0.0,NfnI*d*Neq*Nvar * sizeof *dFluxViscNumdQR_fIL[dim]);
					for (size_t i = 0, iMax = NfnI*Nvar; i < iMax; i++) {
						dFluxViscNumdQL_fIL[dim][iMax*dim+i] = 1.0;
						dFluxViscNumdQR_fIL[dim][iMax*dim+i] = 1.0;
					}
				}
			} else {
				EXIT_UNSUPPORTED;
			}

			for (size_t dim = 0; dim < d; dim++) {
				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdQL_fIL[dim],dnFluxViscNumdQL_fIL[dim]);
				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdQR_fIL[dim],dnFluxViscNumdQR_fIL[dim]);

				for (size_t i = 0; i < NfnI*Neq*Nvar; i++) {
					dnFluxViscNumdQL_fIL[dim][i] *= 0.5;
					dnFluxViscNumdQR_fIL[dim][i] *= 0.5;
				}
			}
		}
	}
	free(FLUXDATA);

	dot_with_normal(NfnI,Neq,n_fIL,FluxViscNum_fIL,nFluxViscNum_fIL);
	free(FluxViscNum_fIL);

	if (imex_type == 'I') {
		free(dFluxViscNumdWL_fIL);
		free(dFluxViscNumdWR_fIL);
		array_free2_d(d,dFluxViscNumdQL_fIL);
		array_free2_d(d,dFluxViscNumdQR_fIL);
	}

	// Modify nFluxViscNum_fIL to account for boundary conditions (if necessary)
	unsigned int const BC = FACE->BC;

	if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
		// Set last component of nFluxViscNum to zero
		size_t const eq = Neq-1;
		for (size_t n = 0; n < NfnI; n++)
			nFluxViscNum_fIL[n+NfnI*eq] = 0.0;

		if (imex_type == 'I') {
			for (size_t n = 0; n < NfnI*Nvar; n++) {
				dnFluxViscNumdWL_fIL[n+NfnI*Nvar*eq] = 0.0;
				dnFluxViscNumdWR_fIL[n+NfnI*Nvar*eq] = 0.0;
			}

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t n = 0; n < NfnI*Nvar; n++) {
					dnFluxViscNumdQL_fIL[dim][n+NfnI*Nvar*eq] = 0.0;
					dnFluxViscNumdQR_fIL[dim][n+NfnI*Nvar*eq] = 0.0;
				}
			}
		}
	}

	if (imex_type == 'I' && Boundary) {
		if (!(BC % BC_STEP_SC == BC_NOSLIP_T         ||
		      BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC ||
		      BC % BC_STEP_SC == BC_DIRICHLET        ||
		      BC % BC_STEP_SC == BC_NEUMANN)) {
			printf("%d\n",BC);
			EXIT_UNSUPPORTED;
			// Update the boundary treatment if the assumptions on the treatment for the linearization with respect to
			// Q are not valid.
		}
	}

	// Include the BC information in dnFluxViscNumdWL_fIL if on a boundary.
	// Not currently required for the Poisson solver as dFNumdW = 0.
	bool update_dnFluxViscNumdWL_fIL = 0;
	if (strstr(PDE,"NavierStokes")) {
		update_dnFluxViscNumdWL_fIL = 1;
	} else if (strstr(PDE,"Poisson")) {
		update_dnFluxViscNumdWL_fIL = 0;
	} else {
		EXIT_UNSUPPORTED;
	}

	if (imex_type == 'I' && Boundary && update_dnFluxViscNumdWL_fIL) {
		// This is only done for dWBdW, the contribution from dQd(W/Q) is assumed to have been properly treated above.
		unsigned int const BC = FACE->BC;

		double const *const XYZ_fIL    = FACE->XYZ_fI;
		double       *const dWBdWL_fIL = malloc(NfnI*Nvar*Nvar * sizeof *dWBdWL_fIL); // free

		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;
		BCdata->BC  = BC;

		BCdata->XYZ = XYZ_fIL;
		BCdata->nL  = n_fIL;
		BCdata->WL  = WL_fIL;
		BCdata->QL  = NULL;

		BCdata->dWBdWL = dWBdWL_fIL;

		evaluate_jacobian_boundary(BCdata,0);
		free(BCdata);

		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			size_t const InddnFdWL = (eq*Nvar+var)*NfnI;

			for (size_t i = 0; i < Nvar; i++) {
				size_t const InddnFdWR = (eq*Neq+i)*NfnI,
				             InddWBdWL = (var*Nvar+i)*NfnI;
				for (size_t n = 0; n < NfnI; n++)
					dnFluxViscNumdWL_fIL[InddnFdWL+n] += dnFluxViscNumdWR_fIL[InddnFdWR+n]*dWBdWL_fIL[InddWBdWL+n];
			}
		}}
		free(dWBdWL_fIL);
	}
}

void add_Jacobian_scaling_FACE(struct S_FDATA const *const FDATA, char const imex_type, char const coef_type)
{
	/*
	 *	Purpose:
	 *		Add the Jacobian scaling introduced from transferring to the reference FACE for integration (with
	 *		orientation as seen from the left VOLUME).
	 *
	 *	Comments:
	 *		Potentially don't need to add scaling to d*d*R if on boundary. CHANGE BELOW (ToBeDeleted)
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 *		coef_type : (coef)icient (type) to be updated.
	 *		            Options:
	 *		            	'W' : Flux for inviscid contribution
	 *		            	'Q' : Flux for viscous  contribution to the 1st equation of the mixed form
	 *		            	'V' : Flux for viscous  contribution to the 2nd equation of the mixed form
	 *		            	'P' : Flux Jacobians for viscous contribution to 2nd equation of the mixed form from Q
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Neq      = DB.Neq,
	                   Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const detJF_fIL = FACE->detJF_fI;

	if (coef_type == 'W' || coef_type == 'V') {
		double *nFluxNum_fIL, *dnFluxNumdWL_fIL, *dnFluxNumdWR_fIL;
		if (coef_type == 'W') {
			nFluxNum_fIL     = (double *const) FDATA->NFluxData->nFluxNum_fI,
			dnFluxNumdWL_fIL = (double *const) FDATA->NFluxData->dnFluxNumdWL_fI,
			dnFluxNumdWR_fIL = (double *const) FDATA->NFluxData->dnFluxNumdWR_fI;
		} else if (coef_type == 'V') {
			nFluxNum_fIL     = (double *const) FDATA->NFluxData->nFluxViscNum_fI,
			dnFluxNumdWL_fIL = (double *const) FDATA->NFluxData->dnFluxViscNumdWL_fI,
			dnFluxNumdWR_fIL = (double *const) FDATA->NFluxData->dnFluxViscNumdWR_fI;
		}

		for (size_t eq = 0; eq < Neq; eq++) {
			size_t const IndnF = eq*NfnI;
			for (size_t n = 0; n < NfnI; n++)
				nFluxNum_fIL[IndnF+n] *= detJF_fIL[n];

			if (imex_type == 'I') {
				for (size_t var = 0; var < Nvar; var++) {
					size_t const InddnFdWIn = (eq*Nvar+var)*NfnI;
					for (size_t n = 0; n < NfnI; n++) {
						dnFluxNumdWL_fIL[InddnFdWIn+n] *= detJF_fIL[n];
						dnFluxNumdWR_fIL[InddnFdWIn+n] *= detJF_fIL[n];
					}
				}
			}
		}
	} else if (coef_type == 'Q') {
		double *const *const nSolNum_fIL     = FDATA->NFluxData->nSolNum_fI,
		       *const *const dnSolNumdWL_fIL = FDATA->NFluxData->dnSolNumdWL_fI,
		       *const *const dnSolNumdWR_fIL = FDATA->NFluxData->dnSolNumdWR_fI;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const IndnF = eq*NfnI;
				for (size_t n = 0; n < NfnI; n++)
					nSolNum_fIL[dim][IndnF+n] *= detJF_fIL[n];
			}
		}

		if (imex_type == 'I') {
			unsigned int eqMax, varMax;
			if (FACE->Boundary) {
				eqMax  = DB.Neq;
				varMax = DB.Nvar;
			} else {
				eqMax  = 1;
				varMax = 1;
			}

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t eq = 0; eq < eqMax; eq++) {
				for (size_t var = 0; var < varMax; var++) {
					size_t const InddnFdWIn = (eq*Nvar+var)*NfnI;
					for (size_t n = 0; n < NfnI; n++) {
						dnSolNumdWL_fIL[dim][InddnFdWIn+n] *= detJF_fIL[n];
						dnSolNumdWR_fIL[dim][InddnFdWIn+n] *= detJF_fIL[n];
					}
				}}
			}
		}
	} else if (coef_type == 'P') {
		double *const *const dnFluxViscNumdQL_fIL = FDATA->NFluxData->dnFluxViscNumdQL_fI,
		       *const *const dnFluxViscNumdQR_fIL = FDATA->NFluxData->dnFluxViscNumdQR_fI;

		if (imex_type == 'E') {
			EXIT_UNSUPPORTED;
		} else if (imex_type == 'I') {
			if (FACE->Boundary) {
				for (size_t dim = 0; dim < d; dim++) {
					for (size_t eq = 0; eq < Neq; eq++) {
					for (size_t varQ = 0; varQ < Nvar; varQ++) {
						size_t const InddnFdQ = (eq*Nvar+varQ)*NfnI;
						for (size_t n = 0; n < NfnI; n++) {
							dnFluxViscNumdQL_fIL[dim][InddnFdQ+n] *= detJF_fIL[n];
						}
					}}
				}
			} else {
				for (size_t dim = 0; dim < d; dim++) {
					for (size_t eq = 0; eq < Neq; eq++) {
					for (size_t varQ = 0; varQ < Nvar; varQ++) {
						size_t const InddnFdQ = (eq*Nvar+varQ)*NfnI;
						for (size_t n = 0; n < NfnI; n++) {
							dnFluxViscNumdQL_fIL[dim][InddnFdQ+n] *= detJF_fIL[n];
							dnFluxViscNumdQR_fIL[dim][InddnFdQ+n] *= detJF_fIL[n];
						}
					}}
				}
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void swap_FACE_orientation(struct S_FDATA const *const FDATA, char const imex_type, char const coef_type)
{
	/*
	 *	Purpose:
	 *		Change orientation of FACE operators to correspond to that of the right VOLUME.
	 *
	 *	Comments:
	 *		Note that the arrays are negated to account for the normal being negative when seen by the opposite VOLUME.
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 *		coef_type : (coef)icient (type) to be updated.
	 *		            Options:
	 *		            	'W' : Flux for inviscid contribution
	 *		            	'Q' : Flux for viscous  contribution to the 1st equation of the mixed form
	 *		            	'V' : Flux for viscous  contribution to the 2nd equation of the mixed form
	 *		            	'P' : Flux Jacobians for viscous contribution to 2nd equation of the mixed form from Q
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d             = DB.d,
	                   Neq           = DB.Neq,
	                   Nvar          = DB.Nvar,
	                   IndFType      = FDATA->IndFType,
	                   NfnI          = OPS[IndFType]->NfnI,
	                   *const nOrdLR = OPS[IndFType]->nOrdLR;

	if (coef_type == 'W' || coef_type == 'V') {
		double *nFluxNum_fI, *dnFluxNumdWL_fI, *dnFluxNumdWR_fI;
		if (coef_type == 'W') {
			nFluxNum_fI     = (double *const) FDATA->NFluxData->nFluxNum_fI,
			dnFluxNumdWL_fI = (double *const) FDATA->NFluxData->dnFluxNumdWL_fI,
			dnFluxNumdWR_fI = (double *const) FDATA->NFluxData->dnFluxNumdWR_fI;
		} else if (coef_type == 'V') {
			nFluxNum_fI     = (double *const) FDATA->NFluxData->nFluxViscNum_fI,
			dnFluxNumdWL_fI = (double *const) FDATA->NFluxData->dnFluxViscNumdWL_fI,
			dnFluxNumdWR_fI = (double *const) FDATA->NFluxData->dnFluxViscNumdWR_fI;
		}

		if (imex_type == 'E') {
			for (size_t i = 0, iMax = Neq*NfnI; i < iMax; i++)
				nFluxNum_fI[i] *= -1.0;

			array_rearrange_d(NfnI,Neq,nOrdLR,'C',nFluxNum_fI);
		} else if (imex_type == 'I') {
			for (size_t i = 0, iMax = Neq*Nvar*NfnI; i < iMax; i++) {
				dnFluxNumdWL_fI[i] *= -1.0;
				dnFluxNumdWR_fI[i] *= -1.0;
			}

			array_rearrange_d(NfnI,Neq*Nvar,nOrdLR,'C',dnFluxNumdWL_fI);
			array_rearrange_d(NfnI,Neq*Nvar,nOrdLR,'C',dnFluxNumdWR_fI);
		}
	} else if (coef_type == 'Q') {
		if (imex_type == 'E') {
			double *const *const nSolNum_fI = FDATA->NFluxData->nSolNum_fI;

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0, iMax = Neq*NfnI; i < iMax; i++)
					nSolNum_fI[dim][i] *= -1.0;

				array_rearrange_d(NfnI,Neq,nOrdLR,'C',nSolNum_fI[dim]);
			}
		} else if (imex_type == 'I') {
			double *const *const dnSolNumdWL_fI = FDATA->NFluxData->dnSolNumdWL_fI,
			       *const *const dnSolNumdWR_fI = FDATA->NFluxData->dnSolNumdWR_fI;

			unsigned int eqMax, varMax;
			if (FACE->Boundary) {
				eqMax  = d+2;
				varMax = d+2;
			} else {
				eqMax  = 1;
				varMax = 1;
			}

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0, iMax = eqMax*varMax*NfnI; i < iMax; i++) {
					dnSolNumdWL_fI[dim][i] *= -1.0;
					dnSolNumdWR_fI[dim][i] *= -1.0;
				}

				array_rearrange_d(NfnI,eqMax*varMax,nOrdLR,'C',dnSolNumdWL_fI[dim]);
				array_rearrange_d(NfnI,eqMax*varMax,nOrdLR,'C',dnSolNumdWR_fI[dim]);
			}
		}
	} else if (coef_type == 'P') {
		double *const *const dnFluxViscNumdQL_fI = FDATA->NFluxData->dnFluxViscNumdQL_fI,
		       *const *const dnFluxViscNumdQR_fI = FDATA->NFluxData->dnFluxViscNumdQR_fI;

		if (imex_type == 'E') {
			EXIT_UNSUPPORTED;
		} else if (imex_type == 'I') {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0, iMax = Neq*Nvar*NfnI; i < iMax; i++) {
					dnFluxViscNumdQL_fI[dim][i] *= -1.0;
					dnFluxViscNumdQR_fI[dim][i] *= -1.0;
				}

				array_rearrange_d(NfnI,Neq*Nvar,nOrdLR,'C',dnFluxViscNumdQL_fI[dim]);
				array_rearrange_d(NfnI,Neq*Nvar,nOrdLR,'C',dnFluxViscNumdQR_fI[dim]);
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void compute_LHS_FACE_Inviscid_Weak(unsigned int const NRows, unsigned int const NCols, unsigned int const Nn,
                                           double const *const I_FF, double const *const dnFluxNumdW_fI,
                                           double const *const ChiS_fI, double *const IdnFdW, double *const LHS)
{
	unsigned int const Neq  = DB.Neq,
	                   Nvar = DB.Nvar;

	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t var = 0; var < Nvar; var++) {
		size_t const Indeqvar = eq*Nvar+var;
		mm_diag_d(NRows,Nn,&dnFluxNumdW_fI[Indeqvar*Nn],I_FF,IdnFdW,1.0,0.0,'R','R');

		size_t const IndLHS = Indeqvar*NRows*NCols;
		mm_d(CBRM,CBNT,CBNT,NRows,NCols,Nn,1.0,1.0,IdnFdW,ChiS_fI,&LHS[IndLHS]);
	}}
}

void finalize_FACE_Inviscid_Weak(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                 double const *const nANumL_fI, double const *const nANumR_fI, char const side,
                                 char const imex_type, char const coef_type)
{
	/*
	 *	Purpose:
	 *		Finalize the RHS/LHS term contributions from FACEs.
	 *
	 *	Comments:
	 *		Various options are available for performing the computation of RHS terms:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
	 *
	 *		The inefficiency of allocating memory for the computed values and them simply copying to RHS/LHS could be
	 *		avoided if functionality is provided to allow for non-zero beta for sf_apply_?, mm_CTN_CSR_? and mm_CTN_?.
	 *		The added cost is likely negligible compared to the actual operations however. (ToBeModified)
	 *
	 *		This function is called to compute both the inviscid and viscous contributions to the RHS/LHS of the
	 *		equations which is why the general (A) is used instead of a specific name for a given normal flux
	 *		contribution (see Notation). For the viscous contribution, the numerical flux is simply negated. When this
	 *		function is called for explicit runs, nANumR_fI is not used.
	 *
	 *		If it is only desired to gain an understanding of what this contribution is, it is suffient to look only at
	 *		the standard approach (i.e. the last condition in the if/else chain).
	 *
	 *	Notation:
	 *		n[1]Num[2]_fI : (n)ormal [1] (== A) (Num)erical [2] evaluated at (f)ace (I)ntegration nodes
	 *			[1] : Flux, dFluxdW[2]
	 *			[2] : (L)eft, (R)ight
	 */

	if (imex_type == 'E') {
		struct S_FDATA const *FDATA = NULL;
		double               *RHS   = NULL;

		if (side == 'L') {
			FDATA = FDATAL;
			RHS   = FDATA->FACE->RHSL;
		} else if (side == 'R') {
			FDATA = FDATAR;
			RHS   = FDATA->FACE->RHSR;

			swap_FACE_orientation(FDATAR,imex_type,coef_type);
		} else {
			EXIT_UNSUPPORTED;
		}

		struct S_OPERATORS_F const *const *const OPS = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

		unsigned int const d        = DB.d,
		                   Neq      = DB.Neq,
		                   Nvar     = DB.Nvar,
		                   P        = FDATA->P,
		                   Eclass   = FDATA->Eclass,
		                   Vf       = FDATA->Vf,
		                   f        = FDATA->f,
		                   SpOp     = FDATA->SpOp,
		                   IndFType = FDATA->IndFType,
		                   NfnI     = OPS[IndFType]->NfnI,
		                   NvnS     = OPS[0]->NvnS,
		                   *const VFPartUnity         = DB.VFPartUnity,
		                   *const *const *const SF_BE = (unsigned int const *const *const *const) DB.SF_BE;

		unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX], NIn0, NIn1;
		double const *OP[DMAX], *const *OP0, *const *OP1;

		double *const InFNum = malloc(NvnS*Nvar * sizeof *InFNum); // free
		if (Eclass == C_TP && SF_BE[P][0][1]) {
			get_sf_parametersF(OPS[0]->NvnI_SF,OPS[0]->NvnS_SF,OPS[0]->I_Weak_VV_SF,
			                   OPS[0]->NfnI_SF,OPS[0]->NvnS_SF,OPS[0]->I_Weak_FF_SF,NIn,NOut,OP,d,Vf,C_TP);

			if (SpOp) {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 2;
				Diag[f/2] = 0;
			} else {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 0;
			}

			sf_apply_d(nANumL_fI,InFNum,NIn,NOut,Neq,OP,Diag,d);
		} else if (Eclass == C_WEDGE && SF_BE[P][1][1]) {
			if (f < 3) { OP0  = OPS[0]->I_Weak_FF_SF, OP1  = OPS[1]->I_Weak_VV_SF;
			             NIn0 = OPS[0]->NfnI_SF,      NIn1 = OPS[1]->NvnI_SF;
			} else {     OP0  = OPS[0]->I_Weak_VV_SF, OP1  = OPS[1]->I_Weak_FF_SF;
			             NIn0 = OPS[0]->NvnI_SF,      NIn1 = OPS[1]->NfnI_SF; }
			get_sf_parametersF(NIn0,OPS[0]->NvnS_SF,OP0,NIn1,OPS[1]->NvnS_SF,OP1,NIn,NOut,OP,d,Vf,C_WEDGE);

			if (SpOp) {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 2;
				if (f < 3)
					Diag[0] = 0;
				else
					Diag[2] = 0;
			} else {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 0;
				Diag[1] = 2;
			}

			sf_apply_d(nANumL_fI,InFNum,NIn,NOut,Neq,OP,Diag,d);
		} else if ((SpOp && (Eclass == C_TP || Eclass == C_WEDGE)) || (VFPartUnity[Eclass])) {
			mm_CTN_CSR_d(OPS[0]->NvnS,Neq,NfnI,OPS[0]->I_Weak_FF_sp[Vf],nANumL_fI,InFNum);
		} else  {
			mm_CTN_d(OPS[0]->NvnS,Neq,NfnI,OPS[0]->I_Weak_FF[Vf],nANumL_fI,InFNum);
		}

		for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
			RHS[i] += InFNum[i];
		free(InFNum);
	} else if (imex_type == 'I') {
		struct S_OPERATORS_F const *const *const OPSL  = (struct S_OPERATORS_F const *const *const) FDATAL->OPS,
		                           *const *const OPSR  = (struct S_OPERATORS_F const *const *const) FDATAR->OPS;
		struct S_FACE        const *const        FACE = FDATAL->FACE;

		unsigned int const VfL      = FDATAL->Vf,
		                   VfR      = FDATAR->Vf,
		                   IndFType = FDATAL->IndFType,
		                   NfnI     = OPSL[IndFType]->NfnI;

		double       *IdnFdW;
		double const *I_FF;

		if (side == 'L') {
			unsigned int const NvnSL = OPSL[0]->NvnS;

			// LHSLL (Effect of (L)eft VOLUME on (L)eft VOLUME)
			I_FF   = OPSL[0]->I_Weak_FF[VfL];
			IdnFdW = malloc(NvnSL*NfnI * sizeof *IdnFdW); // free

			double const *const ChiSL_fIL = OPSL[0]->ChiS_fI[VfL];
			compute_LHS_FACE_Inviscid_Weak(NvnSL,NvnSL,NfnI,I_FF,nANumL_fI,ChiSL_fIL,IdnFdW,FACE->LHSLL);

			free(IdnFdW);
		} else if (side == 'R') {
			double const * ChiS_fI;

			unsigned int const NvnSL         = OPSL[0]->NvnS,
			                   NvnSR         = OPSR[0]->NvnS,
			                   IndFType      = FDATAL->IndFType,
			                   *const nOrdLR = OPSL[IndFType]->nOrdLR,
			                   *const nOrdRL = OPSL[IndFType]->nOrdRL;

			I_FF   = OPSL[0]->I_Weak_FF[VfL];
			IdnFdW = malloc(NvnSL*NfnI * sizeof *IdnFdW); // free

			// LHSRL (Effect of (R)ight VOLUME on (L)eft VOLUME)
			double *const ChiSR_fIL = malloc(NvnSR*NfnI * sizeof *ChiSR_fIL); // free

			ChiS_fI = OPSR[0]->ChiS_fI[VfR];
			for (size_t i = 0; i < NfnI; i++) {
			for (size_t j = 0; j < NvnSR; j++) {
				ChiSR_fIL[i*NvnSR+j] = ChiS_fI[nOrdRL[i]*NvnSR+j];
			}}

			compute_LHS_FACE_Inviscid_Weak(NvnSL,NvnSR,NfnI,I_FF,nANumR_fI,ChiSR_fIL,IdnFdW,FACE->LHSRL);

			free(ChiSR_fIL);
			free(IdnFdW);

			// Swap orientation of numerical flux Jacobian terms
			swap_FACE_orientation(FDATAR,imex_type,coef_type);

			I_FF   = OPSR[0]->I_Weak_FF[VfR];
			IdnFdW = malloc(NvnSR*NfnI * sizeof *IdnFdW); // free

			// LHSLR (Effect of (L)eft VOLUME on (R)ight VOLUME)
			double *const ChiSL_fIR = malloc(NvnSL*NfnI * sizeof *ChiSL_fIR); // free

			ChiS_fI = OPSL[0]->ChiS_fI[VfL];
			for (size_t i = 0; i < NfnI; i++) {
			for (size_t j = 0; j < NvnSL; j++) {
				ChiSL_fIR[i*NvnSL+j] = ChiS_fI[nOrdLR[i]*NvnSL+j];
			}}

			compute_LHS_FACE_Inviscid_Weak(NvnSR,NvnSL,NfnI,I_FF,nANumL_fI,ChiSL_fIR,IdnFdW,FACE->LHSLR);

			free(ChiSL_fIR);

			// LHSRR (Effect of (R)ight VOLUME on (R)ight VOLUME)
			double const *const ChiSR_fIR = OPSR[0]->ChiS_fI[VfR];

			compute_LHS_FACE_Inviscid_Weak(NvnSR,NvnSR,NfnI,I_FF,nANumR_fI,ChiSR_fIR,IdnFdW,FACE->LHSRR);

			free(IdnFdW);
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void compute_LHS_QhatF_Weak(unsigned int const NRows, unsigned int const NCols, unsigned int const Nn,
                                   double const *const I_FF, double const *const *const dnSolNumdW_fI,
                                   double const *const ChiS_fI, double *const IdnFdW, double *const *const Qhat_What,
                                   bool const Boundary)
{
	unsigned int const d = DB.d;

	unsigned int eqMax, varMax;
	if (Boundary) {
		eqMax  = DB.Neq;
		varMax = DB.Nvar;
	} else {
		eqMax  = 1;
		varMax = 1;
	}

	for (size_t dim = 0; dim < d; dim++) {
	for (size_t eq = 0; eq < eqMax; eq++) {
	for (size_t var = 0; var < varMax; var++) {
		size_t const Indeqvar = eq*varMax+var;
		mm_diag_d(NRows,Nn,&dnSolNumdW_fI[dim][Indeqvar*Nn],I_FF,IdnFdW,1.0,0.0,'R','R');

		size_t const IndQ_W = Indeqvar*NRows*NCols;
		// Note that there is a minus sign included in the definition of I_Weak_FF.
		mm_d(CBRM,CBNT,CBNT,NRows,NCols,Nn,-1.0,0.0,IdnFdW,ChiS_fI,&Qhat_What[dim][IndQ_W]);
	}}}
}

void finalize_QhatF_Weak(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR, char const side,
                         char const imex_type)
{
	/*
	 *	Purpose:
	 *		Finalize the RHS/LHS term contributions to Qhat from FACEs.
	 *
	 *	Comments:
	 *		Except for the "-ve" sign, the need to compute d terms instead of 1, and the storage of the use of different
	 *		arrays to be operated on and stored in, this function is identical to finalize_FACE_Inviscid_Weak.
	 *
	 *		The local contribution is added to the numerical flux before evaluating QhatF.
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_FACE const *const FACE = FDATAL->FACE;

	if (imex_type == 'E') {
		correct_numerical_solution_strong(FDATAL,'E',side);

		struct S_FDATA const *FDATA = NULL;

		double *const *QhatF;
		if (side == 'L') {
			FDATA = FDATAL;
			QhatF = (double *const *const) FACE->QhatL;
		} else if (side == 'R') {
			FDATA = FDATAR;
			QhatF = (double *const *const) FACE->QhatR;
			swap_FACE_orientation(FDATAR,imex_type,'Q');
		} else {
			EXIT_UNSUPPORTED;
		}


		struct S_OPERATORS_F const *const *const OPS = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

		unsigned int const d        = DB.d,
		                   Neq      = DB.Neq,
		                   Nvar     = DB.Nvar,
		                   P        = FDATA->P,
		                   Eclass   = FDATA->Eclass,
		                   Vf       = FDATA->Vf,
		                   f        = FDATA->f,
		                   SpOp     = FDATA->SpOp,
		                   IndFType = FDATA->IndFType,
		                   NfnI     = OPS[IndFType]->NfnI,
		                   NvnS     = OPS[0]->NvnS,
		                   *const VFPartUnity         = DB.VFPartUnity,
		                   *const *const *const SF_BE = (unsigned int const *const *const *const) DB.SF_BE;

		double const *const *const nSolNum_fI = (double const *const *const) FDATA->NFluxData->nSolNum_fI;

		unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX], NIn0, NIn1;
		double const *OP[DMAX], *const *OP0, *const *OP1;

		double *const InSNum = malloc(NvnS*Nvar * sizeof *InSNum); // free
		if (Eclass == C_TP && SF_BE[P][0][1]) {
			get_sf_parametersF(OPS[0]->NvnI_SF,OPS[0]->NvnS_SF,OPS[0]->I_Weak_VV_SF,
			                   OPS[0]->NfnI_SF,OPS[0]->NvnS_SF,OPS[0]->I_Weak_FF_SF,NIn,NOut,OP,d,Vf,C_TP);

			if (SpOp) {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 2;
				Diag[f/2] = 0;
			} else {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 0;
			}

			for (size_t dim = 0; dim < d; dim++) {
				// Note that there is a minus sign included in the definition of I_Weak_FF.
				sf_apply_d(nSolNum_fI[dim],InSNum,NIn,NOut,Neq,OP,Diag,d);
				for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
					QhatF[dim][i] = -InSNum[i];
			}
		} else if (Eclass == C_WEDGE && SF_BE[P][1][1]) {
			if (f < 3) { OP0  = OPS[0]->I_Weak_FF_SF, OP1  = OPS[1]->I_Weak_VV_SF;
			             NIn0 = OPS[0]->NfnI_SF,      NIn1 = OPS[1]->NvnI_SF;
			} else {     OP0  = OPS[0]->I_Weak_VV_SF, OP1  = OPS[1]->I_Weak_FF_SF;
			             NIn0 = OPS[0]->NvnI_SF,      NIn1 = OPS[1]->NfnI_SF; }
			get_sf_parametersF(NIn0,OPS[0]->NvnS_SF,OP0,NIn1,OPS[1]->NvnS_SF,OP1,NIn,NOut,OP,d,Vf,C_WEDGE);

			if (SpOp) {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 2;
				if (f < 3)
					Diag[0] = 0;
				else
					Diag[2] = 0;
			} else {
				for (size_t dim = 0; dim < d; dim++)
					Diag[dim] = 0;
				Diag[1] = 2;
			}

			for (size_t dim = 0; dim < d; dim++) {
				// Note that there is a minus sign included in the definition of I_Weak_FF.
				sf_apply_d(nSolNum_fI[dim],InSNum,NIn,NOut,Neq,OP,Diag,d);
				for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
					QhatF[dim][i] = -InSNum[i];
			}
		} else if ((SpOp && (Eclass == C_TP || Eclass == C_WEDGE)) || (VFPartUnity[Eclass])) {
			for (size_t dim = 0; dim < d; dim++) {
				// Note that there is a minus sign included in the definition of I_Weak_FF.
				mm_CTN_CSR_d(OPS[0]->NvnS,Neq,NfnI,OPS[0]->I_Weak_FF_sp[Vf],nSolNum_fI[dim],InSNum);
				for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
					QhatF[dim][i] = -InSNum[i];
			}
		} else  {
			for (size_t dim = 0; dim < d; dim++) {
				// Note that there is a minus sign included in the definition of I_Weak_FF.
				mm_CTN_d(OPS[0]->NvnS,Neq,NfnI,OPS[0]->I_Weak_FF[Vf],nSolNum_fI[dim],InSNum);
				for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
					QhatF[dim][i] = -InSNum[i];
			}
		}
		free(InSNum);
	} else if (imex_type == 'I') {
		struct S_OPERATORS_F const *const *const OPSL = (struct S_OPERATORS_F const *const *const) FDATAL->OPS,
		                           *const *const OPSR = (struct S_OPERATORS_F const *const *const) FDATAR->OPS;

		unsigned int const VfL = FDATAL->Vf,
		                   VfR = FDATAR->Vf,
		                   Boundary = FACE->Boundary,
		                   IndFType = FDATAL->IndFType,
		                   NfnI     = OPSL[IndFType]->NfnI;

		double       *IdnFdW;
		double const *I_FF;

		if (side == 'L') {
			correct_numerical_solution_strong(FDATAL,'I',side);

			unsigned int const NvnSL = OPSL[0]->NvnS;

			double const *const *const dnSolNumdWL_fI = (double const *const *const) FDATAL->NFluxData->dnSolNumdWL_fI;

			// QhatL_WhatL (Effect of (L)eft VOLUME on (L)eft VOLUME)
			I_FF   = OPSL[0]->I_Weak_FF[VfL];
			IdnFdW = malloc(NvnSL*NfnI * sizeof *IdnFdW); // free

			double const *const ChiSL_fIL = OPSL[0]->ChiS_fI[VfL];
			compute_LHS_QhatF_Weak(NvnSL,NvnSL,NfnI,I_FF,dnSolNumdWL_fI,ChiSL_fIL,IdnFdW,FACE->QhatL_WhatL,Boundary);

			free(IdnFdW);
		} else if (side == 'R') {
			if (FACE->Boundary)
				EXIT_UNSUPPORTED;

			double const *      ChiS_fI;
			double const *const *const dnSolNumdWL_fI = (double const *const *const) FDATAL->NFluxData->dnSolNumdWL_fI;
			double const *const *const dnSolNumdWR_fI = (double const *const *const) FDATAL->NFluxData->dnSolNumdWR_fI;

			unsigned int const NvnSL         = OPSL[0]->NvnS,
			                   NvnSR         = OPSR[0]->NvnS,
			                   IndFType      = FDATAL->IndFType,
			                   *const nOrdLR = OPSL[IndFType]->nOrdLR,
			                   *const nOrdRL = OPSL[IndFType]->nOrdRL;

			I_FF   = OPSL[0]->I_Weak_FF[VfL];
			IdnFdW = malloc(NvnSL*NfnI * sizeof *IdnFdW); // free

			// QhatL_WhatR (Effect of (R)ight VOLUME on (L)eft VOLUME)
			double *const ChiSR_fIL = malloc(NvnSR*NfnI * sizeof *ChiSR_fIL); // free

			ChiS_fI = OPSR[0]->ChiS_fI[VfR];
			for (size_t i = 0; i < NfnI; i++) {
			for (size_t j = 0; j < NvnSR; j++) {
				ChiSR_fIL[i*NvnSR+j] = ChiS_fI[nOrdRL[i]*NvnSR+j];
			}}

			compute_LHS_QhatF_Weak(NvnSL,NvnSR,NfnI,I_FF,dnSolNumdWR_fI,ChiSR_fIL,IdnFdW,FACE->QhatL_WhatR,Boundary);

			free(ChiSR_fIL);
			free(IdnFdW);

			// Swap orientation of numerical flux Jacobian terms
			correct_numerical_solution_strong(FDATAR,'I',side);
			swap_FACE_orientation(FDATAR,'I','Q');

			I_FF   = OPSR[0]->I_Weak_FF[VfR];
			IdnFdW = malloc(NvnSR*NfnI * sizeof *IdnFdW); // free

			// QhatR_WhatL (Effect of (L)eft VOLUME on (R)ight VOLUME)
			double *const ChiSL_fIR = malloc(NvnSL*NfnI * sizeof *ChiSL_fIR); // free

			ChiS_fI = OPSL[0]->ChiS_fI[VfL];
			for (size_t i = 0; i < NfnI; i++) {
			for (size_t j = 0; j < NvnSL; j++) {
				ChiSL_fIR[i*NvnSL+j] = ChiS_fI[nOrdLR[i]*NvnSL+j];
			}}

			compute_LHS_QhatF_Weak(NvnSR,NvnSL,NfnI,I_FF,dnSolNumdWL_fI,ChiSL_fIR,IdnFdW,FACE->QhatR_WhatL,Boundary);

			free(ChiSL_fIR);

			// QhatR_WhatR (Effect of (R)ight VOLUME on (R)ight VOLUME)
			double const *const ChiSR_fIR = OPSR[0]->ChiS_fI[VfR];

			compute_LHS_QhatF_Weak(NvnSR,NvnSR,NfnI,I_FF,dnSolNumdWR_fI,ChiSR_fIR,IdnFdW,FACE->QhatR_WhatR,Boundary);

			free(IdnFdW);
		} else {
			EXIT_UNSUPPORTED;
		}
	}
}

void finalize_FACE_Viscous_Weak(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                double *const nANumL_fI, double *const nANumR_fI, char const side, char const imex_type,
                                char const coef_type)
{
	/*
	 *	Purpose:
	 *		Finalize the viscous RHS/LHS term contributions from FACEs.
	 *
	 *	Comments:
	 *		The required operation is identical to that of finalize_FACE_Inviscid_Weak as the viscous flux has already
	 *		been negated.
	 */

	finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,nANumL_fI,nANumR_fI,side,imex_type,coef_type);
}

static void compute_LHS_FACE_Q_Weak(unsigned int const NRows, unsigned int const NCols, unsigned int const Nn,
                                    double const *const I_FF, double const *const *const dnFluxViscNumdQ_fI,
                                    double const *const *const Q_What, double *const IdnFdQ, double *const LHS,
                                    bool const Boundary)
{
	unsigned int const d    = DB.d,
	                   Neq  = DB.Neq,
	                   Nvar = DB.Nvar;

	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t varQ = 0; varQ < Nvar; varQ++) {
		for (size_t dim = 0; dim < d; dim++) {
			size_t const InddFdQ = (eq*Nvar+varQ)*Nn;
			mm_diag_d(NRows,Nn,&dnFluxViscNumdQ_fI[dim][InddFdQ],I_FF,IdnFdQ,1.0,0.0,'R','R');

			for (size_t var = 0; var < Nvar; var++) {
				if (!Boundary && (var != varQ))
					continue; // dQ/dWhat is block diagonal when not on a boundary.

				size_t IndQ_What = 0;
				if (Boundary)
					IndQ_What = (varQ*Nvar+var)*Nn*NCols;
				else
					IndQ_What = 0;

				size_t const IndLHS = (eq*Nvar+var)*NRows*NCols;
				mm_d(CBRM,CBNT,CBNT,NRows,NCols,Nn,1.0,1.0,IdnFdQ,&Q_What[dim][IndQ_What],&LHS[IndLHS]);
			}
		}
	}}
}

void finalize_implicit_FACE_Q_Weak(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                   char const side)
{
	/*
	 *	Purpose:
	 *		Finalize the LHS contribution to FACE->LHS from linearization of the numerical flux with respect to Q.
	 *
	 *	Comments:
	 *		Again, very similar to other finalize_* functions where the linearization with respect to Q is used and the
	 *		right operator is given by Qp_What (as opposed to ChiS_fI).
	 */

	struct S_FACE const *const FACE = FDATAL->FACE;

	struct S_OPERATORS_F const *const *const OPSL = (struct S_OPERATORS_F const *const *const) FDATAL->OPS,
	                           *const *const OPSR = (struct S_OPERATORS_F const *const *const) FDATAR->OPS;
	struct S_NumericalFlux const *const NFluxData = (struct S_NumericalFlux const *const) FDATAL->NFluxData;

	unsigned int const d        = DB.d,
	                   VfL      = FDATAL->Vf,
	                   VfR      = FDATAR->Vf,
	                   Boundary = FACE->Boundary,
	                   IndFType = FDATAL->IndFType,
	                   NfnI     = OPSL[IndFType]->NfnI;

	double       *IdnFdQ;
	double const *I_FF;

	if (side == 'L') {
		unsigned int const NvnSL = OPSL[0]->NvnS;

		double const *const *const dnFluxViscNumdQL_fI = (double const *const *const) NFluxData->dnFluxViscNumdQL_fI;

		// QL_WhatL (Effect of QL(WL,...))
		I_FF   = OPSL[0]->I_Weak_FF[VfL];
		IdnFdQ = malloc(NvnSL*NfnI * sizeof *IdnFdQ); // free

		double const *const *const Qp_What = (double const *const *const) FDATAL->Qp_WhatL;
		compute_LHS_FACE_Q_Weak(NvnSL,NvnSL,NfnI,I_FF,dnFluxViscNumdQL_fI,Qp_What,IdnFdQ,FACE->LHSLL,Boundary);
		if (Boundary)
			array_free2_d(d,(double **) Qp_What);
		free(IdnFdQ);
	} else if (side == 'R') {
		if (Boundary)
			EXIT_UNSUPPORTED;

		unsigned int const *const nOrdRL = OPSL[IndFType]->nOrdRL,
		                   *const nOrdLR = OPSL[IndFType]->nOrdLR;

		double const *const *      Qp_What;
		double const *const *const dnFluxViscNumdQL_fI = (double const *const *const) NFluxData->dnFluxViscNumdQL_fI;
		double const *const *const dnFluxViscNumdQR_fI = (double const *const *const) NFluxData->dnFluxViscNumdQR_fI;

		unsigned int const NvnSL = OPSL[0]->NvnS,
		                   NvnSR = OPSR[0]->NvnS;

		I_FF   = OPSL[0]->I_Weak_FF[VfL];
		IdnFdQ = malloc(NvnSL*NfnI * sizeof *IdnFdQ); // free

		// Rearrange such that dQRdW(L/R) are seen at the nodes as ordered from the left VOLUME
		for (size_t dim = 0; dim < d; dim++) {
			array_rearrange_d(NfnI,NvnSL,nOrdRL,'R',FDATAR->Qp_WhatL[dim]);
			array_rearrange_d(NfnI,NvnSR,nOrdRL,'R',FDATAR->Qp_WhatR[dim]);
		}

		// LHSLL (Effect of QR(WL,...))
		Qp_What = (double const *const *const) FDATAR->Qp_WhatL;
		compute_LHS_FACE_Q_Weak(NvnSL,NvnSL,NfnI,I_FF,dnFluxViscNumdQR_fI,Qp_What,IdnFdQ,FACE->LHSLL,Boundary);

		// LHSRL (Effect of QL(...,WR))
		Qp_What = (double const *const *const) FDATAL->Qp_WhatR;
		compute_LHS_FACE_Q_Weak(NvnSL,NvnSR,NfnI,I_FF,dnFluxViscNumdQL_fI,Qp_What,IdnFdQ,FACE->LHSRL,Boundary);

		//       (Effect of QR(...,WR))
		Qp_What = (double const *const *const) FDATAR->Qp_WhatR;
		compute_LHS_FACE_Q_Weak(NvnSL,NvnSR,NfnI,I_FF,dnFluxViscNumdQR_fI,Qp_What,IdnFdQ,FACE->LHSRL,Boundary);
		free(IdnFdQ);

		// Rearrange such that dQ(R/L)dW(L/R) are seen at the nodes as ordered from the right VOLUME
		for (size_t dim = 0; dim < d; dim++) {
			array_rearrange_d(NfnI,NvnSL,nOrdLR,'R',FDATAL->Qp_WhatL[dim]);
			array_rearrange_d(NfnI,NvnSR,nOrdLR,'R',FDATAL->Qp_WhatR[dim]);
			array_rearrange_d(NfnI,NvnSL,nOrdLR,'R',FDATAR->Qp_WhatL[dim]);
			array_rearrange_d(NfnI,NvnSR,nOrdLR,'R',FDATAR->Qp_WhatR[dim]);
		}

		// Swap orientation of numerical flux Jacobian terms
		swap_FACE_orientation(FDATAR,'I','P');

		I_FF   = OPSR[0]->I_Weak_FF[VfR];
		IdnFdQ = malloc(NvnSR*NfnI * sizeof *IdnFdQ); // free

		// LHSLR (Effect of QL(WL,...))
		Qp_What = (double const *const *const) FDATAL->Qp_WhatL;
		compute_LHS_FACE_Q_Weak(NvnSR,NvnSL,NfnI,I_FF,dnFluxViscNumdQL_fI,Qp_What,IdnFdQ,FACE->LHSLR,Boundary);

		//       (Effect of QR(WL,...))
		Qp_What = (double const *const *const) FDATAR->Qp_WhatL;
		compute_LHS_FACE_Q_Weak(NvnSR,NvnSL,NfnI,I_FF,dnFluxViscNumdQR_fI,Qp_What,IdnFdQ,FACE->LHSLR,Boundary);

		// LHSRR (Effect of QL(...,WR))
		Qp_What = (double const *const *const) FDATAL->Qp_WhatR;
		compute_LHS_FACE_Q_Weak(NvnSR,NvnSR,NfnI,I_FF,dnFluxViscNumdQL_fI,Qp_What,IdnFdQ,FACE->LHSRR,Boundary);

		//       (Effect of QR(...,WR))
		Qp_What = (double const *const *const) FDATAR->Qp_WhatR;
		compute_LHS_FACE_Q_Weak(NvnSR,NvnSR,NfnI,I_FF,dnFluxViscNumdQR_fI,Qp_What,IdnFdQ,FACE->LHSRR,Boundary);
		free(IdnFdQ);

		array_free2_d(d,(double **) FDATAL->Qp_WhatL);
		array_free2_d(d,(double **) FDATAL->Qp_WhatR);
		array_free2_d(d,(double **) FDATAR->Qp_WhatL);
		array_free2_d(d,(double **) FDATAR->Qp_WhatR);
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_LHSQF_Weak(struct S_FACE *const FACE)
{
	if (!DB.Viscous)
		return;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	struct S_VOLUME *const VL = FACE->VL,
	                *const VR = FACE->VR;

	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t varQ = 0; varQ < Nvar; varQ++) {
	for (size_t var = 0; var < Nvar; var++) {
		unsigned int const NvnSL = VL->NvnS;
		if (FACE->Boundary) {
			size_t const IndLHSQ = (eq*Nvar+varQ)*NvnSL*NvnSL,
			             IndQ_W  = (varQ*Nvar+var)*NvnSL*NvnSL,
			             IndLHS  = (eq*Nvar+var)*NvnSL*NvnSL;
			for (size_t dim = 0; dim < d; dim++)
				mm_d(CBRM,CBNT,CBNT,NvnSL,NvnSL,NvnSL,1.0,1.0,&VL->LHSQ[dim][IndLHSQ],&FACE->QhatL_WhatL[dim][IndQ_W],&VL->LHS[IndLHS]);
		} else {
			// dQhat/dWhat is block diagonal for this term
			if (var != varQ)
				continue;

			size_t const Indeqvar = (eq*Nvar+var);
			size_t IndLHS, IndLHSQ;

			unsigned int const NvnSR = VR->NvnS;
			for (size_t dim = 0; dim < d; dim++) {
				IndLHS = Indeqvar*NvnSL*NvnSL;
				mm_d(CBRM,CBNT,CBNT,NvnSL,NvnSL,NvnSL,1.0,1.0,&VL->LHSQ[dim][IndLHS],FACE->QhatL_WhatL[dim],&VL->LHS[IndLHS]);

				IndLHS = Indeqvar*NvnSL*NvnSR;
				IndLHSQ = Indeqvar*NvnSL*NvnSL;
				mm_d(CBRM,CBNT,CBNT,NvnSL,NvnSR,NvnSL,1.0,1.0,&VL->LHSQ[dim][IndLHSQ],FACE->QhatL_WhatR[dim],&FACE->LHSRL[IndLHS]);

				IndLHSQ = Indeqvar*NvnSR*NvnSR;
				mm_d(CBRM,CBNT,CBNT,NvnSR,NvnSL,NvnSR,1.0,1.0,&VR->LHSQ[dim][IndLHSQ],FACE->QhatR_WhatL[dim],&FACE->LHSLR[IndLHS]);

				IndLHS = Indeqvar*NvnSR*NvnSR;
				mm_d(CBRM,CBNT,CBNT,NvnSR,NvnSR,NvnSR,1.0,1.0,&VR->LHSQ[dim][IndLHS],FACE->QhatR_WhatR[dim],&VR->LHS[IndLHS]);
			}
		}
	}}}
}
