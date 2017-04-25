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
#include "fluxes_inviscid.h"
#include "fluxes_viscous.h"
#include "jacobian_boundary_conditions.h"
#include "jacobian_fluxes_inviscid.h"
#include "jacobian_fluxes_viscous.h"

#include "array_swap.h"
#include "array_free.h"
#include "array_print.h"

#include "array_norm.h" // ToBeDeleted

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

		unsigned int const d        = DB.d,
		                   Nvar     = d+2,
		                   P        = VDATA->P,
		                   Eclass   = VDATA->Eclass,
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
	 *		Convert input (A)rray between (p)hysical to (r)eference space by multiplying by by appropriate (C)ofactor
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
	 *		Various options are available for performing the computation:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
	 *
	 *		Certain multiplications can be avoided when computing LHS terms based on the sparsity of the flux Jacobian,
	 *		usage of this knowledge is currently not incorporated.
	 *
	 *		If it is only desired to gain an understanding of what this contribution is, it is suffient to look only at
	 *		the standard approach (i.e. the last condition in the if/else chain).
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
		unsigned int const Neq  = d+2,
		                   Nvar = d+2;

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
	                   Nvar = d+2,
	                   Neq  = d+2,
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
	IndOrdLR = FACE->IndOrdInOut;
	IndOrdRL = FACE->IndOrdOutIn;

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
		FDATA->VOLUME  = FACE->VIn;
		FDATA->Vf      = FACE->VfIn;
		FDATA->QhatF   = FACE->QhatL;
		FDATA->QhatF_c = FACE->QhatL_c;
	} else if (side == 'R') {
		FDATA->VOLUME  = FACE->VOut;
		FDATA->Vf      = FACE->VfOut;
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
	 *		Increasing the dissipation parameter beyond NfMax appears to potentially give converged L2 errors.
	 *		INVESTIGATE. (ToBeModified)
	 *
	 *		For several of the viscous boundary conditions, including NoSlip_Dirichlet_T and NoSlip_Adiabatic (as they
	 *		are currently implemented), there is no dependence of boundary values on variables other than the the
	 *		variable under consideration (i.e. dWB/dWL is block diagonal). This means that only the block diagonal
	 *		entries of Qhat_What, Q_What are non-zero, resulting in potentially significant cost savings if exploited.
	 *		That non-boundary contributions are already block diagonal may reduce the importance of this observation
	 *		when many interior elements are present. ToBeModified.
	 *
	 *	Notation:
	 *		imex_type : (im)plicit (ex)plicit (type) indicates whether this function is being called for an implicit or
	 *		            explicit computation.
	 *
	 *		The allowed options for coef_type are 'W' (conserved variables), 'Q' ( == GradW here)
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

	unsigned int const d            = DB.d,
	                   Nvar         = d+2,
	                   P            = FDATA->P,
	                   Eclass       = FDATA->Eclass,
	                   Vf           = FDATA->Vf,
	                   f            = FDATA->f,
	                   SpOp         = FDATA->SpOp,
	                   IndFType     = FDATA->IndFType,
	                   NfnI         = OPS[IndFType]->NfnI,
	                   NvnS         = OPS[0]->NvnS,
	                   *const VFPartUnity         = DB.VFPartUnity,
	                   *const *const *const SF_BE = (const unsigned int *const *const *const) DB.SF_BE;

	double **Qhat = NULL;
	if (coef_type == 'Q') {
		Qhat = malloc(d * sizeof *Qhat);
		for (size_t dim = 0; dim < d; dim++) {
			Qhat[dim] = malloc(NvnS*Nvar * sizeof *Qhat[dim]); // free
			for (size_t i = 0; i < NvnS*Nvar; i++)
				Qhat[dim][i] = VOLUME->QhatV[dim][i];
		}

		if (DB.ViscousFluxType == FLUX_CDG2) {
			// Evaluate from which side scaling should be computed based on area switch (Brdar(2012), eq. (4.5))
			if (FDATA->side == 'L') {
				if (FACE->Boundary) {
					FACE->CDG2_side = 'L';
				} else {
					struct S_VOLUME const *const VL = FACE->VIn,
					                      *const VR = FACE->VOut;

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

		// Partially correct Qhat
		if (chi != 0.0) {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0; i < NvnS*Nvar; i++)
					Qhat[dim][i] += chi*FDATA->QhatF[dim][i];
			}
		}

		// Compute linearized contributions of partially corrected Q (Q_What**)
		if (imex_type == 'I') {
			struct S_VOLUME const *const VL = FACE->VIn,
			                      *const VR = FACE->VOut;

			unsigned int const NvnSL = VL->NvnS,
			                   NvnSR = VR->NvnS;

			// Q_What** arrays are freed in finalize_implicit_FACE_Q_Weak.
			if (Boundary) {
				double *const Qhat_What = malloc(NvnSL*NvnSL * sizeof *Qhat_What); // free

				double **Q_What = malloc(d * sizeof *Q_What); // keep
				for (size_t dim = 0; dim < d; dim++) {
					Q_What[dim] = malloc(NfnI*NvnSL*Nvar*Nvar * sizeof *Q_What[dim]); // keep

					for (size_t varQ = 0; varQ < Nvar; varQ++) {
					for (size_t var = 0; var < Nvar; var++) {
						size_t const Indvv = varQ*Nvar+var;
						if (varQ == var) {
							for (size_t i = 0; i < NvnSL*NvnSL; i++)
								Qhat_What[i] = VL->QhatV_What[dim][i];
						} else {
							memset(Qhat_What,0.0,NvnSL*NvnSL * sizeof *Qhat_What);
						}

						if (chi != 0.0) {
							for (size_t i = 0; i < NvnSL*NvnSL; i++)
								Qhat_What[i] += chi*(FACE->Qhat_WhatLL[dim][Indvv*NvnSL*NvnSL+i]);
						}
						mm_d(CBRM,CBNT,CBNT,NfnI,NvnSL,NvnSL,1.0,0.0,OPS[0]->ChiS_fI[Vf],Qhat_What,&Q_What[dim][Indvv*NfnI*NvnSL]);
					}}
				}
				free(Qhat_What);
				FDATA->Q_WhatLL = Q_What;
			} else {
				if (FDATA->side == 'L') {
					double *const Qhat_WhatLL = malloc(NvnSL*NvnSL * sizeof *Qhat_WhatLL); // free

					double **Q_WhatLL = malloc(d * sizeof *Q_WhatLL), // keep
					       **Q_WhatRL = malloc(d * sizeof *Q_WhatRL); // keep
					for (size_t dim = 0; dim < d; dim++) {
						Q_WhatLL[dim] = malloc(NfnI*NvnSL * sizeof *Q_WhatLL[dim]); // keep
						for (size_t i = 0; i < NvnSL*NvnSL; i++)
							Qhat_WhatLL[i] = VL->QhatV_What[dim][i];
						if (chi != 0.0) {
							for (size_t i = 0; i < NvnSL*NvnSL; i++)
								Qhat_WhatLL[i] += chi*(FACE->Qhat_WhatLL[dim][i]);
						}
						mm_d(CBRM,CBNT,CBNT,NfnI,NvnSL,NvnSL,1.0,0.0,OPS[0]->ChiS_fI[Vf],Qhat_WhatLL,Q_WhatLL[dim]);

						Q_WhatRL[dim] = calloc(NfnI*NvnSR , sizeof *Q_WhatRL[dim]); // keep
						if (chi != 0.0)
							mm_d(CBRM,CBNT,CBNT,NfnI,NvnSR,NvnSL,chi,0.0,OPS[0]->ChiS_fI[Vf],FACE->Qhat_WhatRL[dim],Q_WhatRL[dim]);
					}
					free(Qhat_WhatLL);

					FDATA->Q_WhatLL = Q_WhatLL;
					FDATA->Q_WhatRL = Q_WhatRL;
				} else if (FDATA->side == 'R') {
					double *const Qhat_WhatRR = malloc(NvnSR*NvnSR * sizeof *Qhat_WhatRR); // free

					double **Q_WhatLR = malloc(d * sizeof *Q_WhatLR), // keep
					       **Q_WhatRR = malloc(d * sizeof *Q_WhatRR); // keep
					for (size_t dim = 0; dim < d; dim++) {
						Q_WhatLR[dim] = calloc(NfnI*NvnSL , sizeof *Q_WhatLR[dim]); // keep
						if (chi != 0.0)
							mm_d(CBRM,CBNT,CBNT,NfnI,NvnSL,NvnSR,chi,0.0,OPS[0]->ChiS_fI[Vf],FACE->Qhat_WhatLR[dim],Q_WhatLR[dim]);

						Q_WhatRR[dim] = malloc(NfnI*NvnSR * sizeof *Q_WhatRR[dim]); // keep
						for (size_t i = 0; i < NvnSR*NvnSR; i++)
							Qhat_WhatRR[i] = VR->QhatV_What[dim][i];
						if (chi != 0.0) {
							for (size_t i = 0; i < NvnSR*NvnSR; i++)
								Qhat_WhatRR[i] += chi*(FACE->Qhat_WhatRR[dim][i]);
						}
						mm_d(CBRM,CBNT,CBNT,NfnI,NvnSR,NvnSR,1.0,0.0,OPS[0]->ChiS_fI[Vf],Qhat_WhatRR,Q_WhatRR[dim]);
					}
					free(Qhat_WhatRR);

					FDATA->Q_WhatLR = Q_WhatLR;
					FDATA->Q_WhatRR = Q_WhatRR;
				}
			}
		}
	}

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
				sf_apply_d(Qhat[dim],FDATA->GradW_fIL[dim],NIn,NOut,Nvar,OP,Diag,d);
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
				sf_apply_d(Qhat[dim],FDATA->GradW_fIL[dim],NIn,NOut,Nvar,OP,Diag,d);
		}
	} else if (((SpOp && (Eclass == C_TP || Eclass == C_WEDGE)) || (VFPartUnity[Eclass])) && DB.AllowSparseFACE) {
		if (coef_type == 'W') {
			mm_CTN_CSR_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI_sp[Vf],VOLUME->What,FDATA->W_fIL);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				mm_CTN_CSR_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI_sp[Vf],Qhat[dim],FDATA->GradW_fIL[dim]);
		}
	} else {
		if (coef_type == 'W') {
			mm_CTN_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI[Vf],VOLUME->What,FDATA->W_fIL);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				mm_CTN_d(NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI[Vf],Qhat[dim],FDATA->GradW_fIL[dim]);
		}
	}

	if (coef_type == 'Q')
		array_free2_d(d,Qhat);
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
	 *
	 *		The ExactSlipWall flag is present to enable computation of the exact solution, which may be used to
	 *		demonstrate that it is the SlipWall BC which is responsible for the suboptimal convergence using
	 *		isoparametric geometry representation for the Euler equations on curved domains.
	 */

	struct S_OPERATORS_F const *const *const OPS  = FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d             = DB.d,
	                   Nvar          = d+2,
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

		BCdata->XYZ    = XYZ_fIL;
		BCdata->nL     = n_fIL;
		BCdata->WL     = WL_fIL;
		BCdata->WB     = WR_fIL;
		BCdata->GradWL = NULL;
		BCdata->GradWB = NULL;

		if (BC % BC_STEP_SC == BC_RIEMANN) {
			boundary_Riemann(NfnI,1,XYZ_fIL,WL_fIL,NULL,WR_fIL,n_fIL,d);
		} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
			bool const ExactSlipWall = 0;
			if (ExactSlipWall) {
				double *const UR_fIL = malloc(NVAR3D*NfnI * sizeof *UR_fIL); // free

				compute_exact_solution(NfnI,XYZ_fIL,UR_fIL,0);
				convert_variables(UR_fIL,WR_fIL,DMAX,d,NfnI,1,'p','c');

				free(UR_fIL);
			} else {
				boundary_SlipWall(NfnI,1,WL_fIL,WR_fIL,n_fIL,d);
			}
		} else if (BC % BC_STEP_SC == BC_BACKPRESSURE) {
			boundary_BackPressure(NfnI,1,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_TOTAL_TP) {
			boundary_Total_TP(NfnI,1,XYZ_fIL,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_SUPERSONIC_IN) {
			boundary_SupersonicInflow(NfnI,1,XYZ_fIL,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_SUPERSONIC_OUT) {
			boundary_SupersonicOutflow(NfnI,1,XYZ_fIL,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_NOSLIP_T) {
			boundary_NoSlip_Dirichlet(BCdata);
		} else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
			boundary_NoSlip_Adiabatic(BCdata);
		} else {
			EXIT_UNSUPPORTED;
		}
		free(BCdata);
	}
}

void compute_WR_GradWR_fIL(struct S_FDATA const *const FDATA, double const *const WL_fIL, double *const WR_fIL,
                           double const *const *const GradWL_fIL, double *const *const GradWR_fIL, char const imex_type)
{
	/*
	 *	Purpose:
	 *		Compute W and GradW of the (R)ight VOLUME at the FACE cubature nodes corresponding to the (L)eft VOLUME.
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
			array_rearrange_d(NfnI,Nvar,nOrdRL,'C',GradWR_fIL[dim]);
	} else { // Boundary FACE
		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;

		BCdata->XYZ    = FACE->XYZ_fI;
		BCdata->nL     = FACE->n_fI;
		BCdata->WL     = WL_fIL;
		BCdata->WB     = WR_fIL;
		BCdata->GradWL = GradWL_fIL;
		BCdata->GradWB = GradWR_fIL;

		if (BC % BC_STEP_SC == BC_NOSLIP_T) {
			boundary_NoSlip_Dirichlet(BCdata);
		} else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
			boundary_NoSlip_Adiabatic(BCdata);
		} else {
			EXIT_UNSUPPORTED;
		}
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

		BCdata->XYZ    = XYZ_fIL;
		BCdata->nL     = n_fIL;
		BCdata->WL     = WL_fIL;
		BCdata->GradWL = NULL;

		BCdata->dWBdWL = dWRdWL_fIL;

		if (BC % BC_STEP_SC == BC_RIEMANN)
			jacobian_boundary_Riemann(NfnI,1,XYZ_fIL,WL_fIL,NULL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_SLIPWALL)
			jacobian_boundary_SlipWall(NfnI,1,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_BACKPRESSURE)
			jacobian_boundary_BackPressure(NfnI,1,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_TOTAL_TP)
			jacobian_boundary_Total_TP(NfnI,1,XYZ_fIL,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_SUPERSONIC_IN)
			jacobian_boundary_SupersonicInflow(NfnI,1,XYZ_fIL,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_SUPERSONIC_OUT)
			jacobian_boundary_SupersonicOutflow(NfnI,1,XYZ_fIL,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_NOSLIP_T)
			jacobian_boundary_NoSlip_Dirichlet(BCdata);
		else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC)
			jacobian_boundary_NoSlip_Adiabatic(BCdata);
		else
			EXIT_UNSUPPORTED;
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
	                   Nvar     = d+2,
	                   Neq      = d+2,
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
			eqMax  = d+2;
			varMax = d+2;
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

		BCdata->XYZ    = XYZ_fIL;
		BCdata->nL     = n_fIL;
		BCdata->WL     = WL_fIL;
		BCdata->WB     = (double *const) WR_fIL;
		BCdata->GradWL = NULL;
		BCdata->GradWB = NULL;

		BCdata->dWBdWL = dWRdWL_fIL;

		if (BC % BC_STEP_SC == BC_RIEMANN)
			jacobian_boundary_Riemann(NfnI,1,XYZ_fIL,WL_fIL,NULL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_SLIPWALL)
			jacobian_boundary_SlipWall(NfnI,1,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_BACKPRESSURE)
			jacobian_boundary_BackPressure(NfnI,1,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_TOTAL_TP)
			jacobian_boundary_Total_TP(NfnI,1,XYZ_fIL,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_SUPERSONIC_IN)
			jacobian_boundary_SupersonicInflow(NfnI,1,XYZ_fIL,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_SUPERSONIC_OUT)
			jacobian_boundary_SupersonicOutflow(NfnI,1,XYZ_fIL,WL_fIL,dWRdWL_fIL,n_fIL,d,Neq);
		else if (BC % BC_STEP_SC == BC_NOSLIP_T)
			jacobian_boundary_NoSlip_Dirichlet(BCdata);
		else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC)
			jacobian_boundary_NoSlip_Adiabatic(BCdata);
		else
			EXIT_UNSUPPORTED;
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

static void correct_numerical_solution_strong(struct S_FDATA const *const FDATA, char const imex_type, char const side,
                                              char const Form_MF1)
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

	if (Form_MF1 != 'S')
		EXIT_UNSUPPORTED;

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = d+2,
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
				eqMax  = d+2;
				varMax = d+2;
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
	 *		When this is true the full Jacobian can be simplified as follows:
	 *			dnFNumVisc/dWhatL = df/dQ*(dQ/dQL*dQL/dWhatL + dQ/dQR*(dQR/dWL*dWL/dWhatL + dQR/dQL*dQL/dWhatL))
	 *			                  = df/dQ*( 0.5  * Q_WhatLL  +  0.5  *(  0.0  *  ChiS_vI  +    I   * Q_WhatLL ))
	 *			                  = df/dQ*Q_WhatLL
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

	unsigned int const d        = DB.d,
	                   Nvar     = d+2,
	                   Neq      = d+2,
	                   IndFType = FDATAL->IndFType,
	                   Boundary = FACE->Boundary,
	                   NfnI     = OPSL[IndFType]->NfnI;

	double const *const n_fIL                = FACE->n_fI,
	             *const WL_fIL               = FDATAL->NFluxData->WL_fIL,
	             *const WR_fIL               = FDATAL->NFluxData->WR_fIL,
	             *const *const GradWL_fIL    = (double const *const *const) FDATAL->GradW_fIL,
	             *const *const GradWR_fIL    = (double const *const *const) FDATAR->GradW_fIL;
	double       *const nFluxViscNum_fIL     = FDATAL->NFluxData->nFluxViscNum_fI,
	             *const dnFluxViscNumdWL_fIL = FDATAL->NFluxData->dnFluxViscNumdWL_fI,
	             *const dnFluxViscNumdWR_fIL = FDATAL->NFluxData->dnFluxViscNumdWR_fI,
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

	if (Boundary) {
		double *const W_fIL = malloc(NfnI*Nvar * sizeof *W_fIL); // free
		for (size_t i = 0; i < NfnI*Nvar; i++)
			W_fIL[i] = 0.5*(WL_fIL[i]+WR_fIL[i]);

		double **Q_fIL = malloc(d * sizeof *Q_fIL); // free
		for (size_t dim = 0; dim < d; dim++) {
			Q_fIL[dim] = malloc(NfnI*Nvar * sizeof *Q_fIL[dim]); // free
			for (size_t i = 0; i < NfnI*Nvar; i++)
				Q_fIL[dim][i] = 0.5*(GradWL_fIL[dim][i]+GradWR_fIL[dim][i]);
		}

		flux_viscous(NfnI,1,W_fIL,(double const *const *const) Q_fIL,FluxViscNum_fIL);
		if (imex_type == 'I') {
			jacobian_flux_viscous(NfnI,1,W_fIL,(double const *const *const) Q_fIL,dFluxViscNumdWL_fIL,dFluxViscNumdQL_fIL);

			dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdWL_fIL,dnFluxViscNumdWL_fIL);
			for (size_t dim = 0; dim < d; dim++)
				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdQL_fIL[dim],dnFluxViscNumdQL_fIL[dim]);

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
		}
		free(W_fIL);
		array_free2_d(d,Q_fIL);
	} else {
		double *const FluxViscL_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscL_fIL); // free
		double *const FluxViscR_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscR_fIL); // free

		flux_viscous(NfnI,1,WL_fIL,GradWL_fIL,FluxViscL_fIL);
		flux_viscous(NfnI,1,WR_fIL,GradWR_fIL,FluxViscR_fIL);

		for (size_t i = 0; i < NfnI*d*Neq; i++)
			FluxViscNum_fIL[i] = 0.5*(FluxViscL_fIL[i]+FluxViscR_fIL[i]);

		free(FluxViscL_fIL);
		free(FluxViscR_fIL);

		if (imex_type == 'I') {
			jacobian_flux_viscous(NfnI,1,WL_fIL,(double const *const *const) GradWL_fIL,dFluxViscNumdWL_fIL,dFluxViscNumdQL_fIL);
			jacobian_flux_viscous(NfnI,1,WR_fIL,(double const *const *const) GradWR_fIL,dFluxViscNumdWR_fIL,dFluxViscNumdQR_fIL);

			dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdWL_fIL,dnFluxViscNumdWL_fIL);
			dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdWR_fIL,dnFluxViscNumdWR_fIL);
			for (size_t dim = 0; dim < d; dim++) {
				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdQL_fIL[dim],dnFluxViscNumdQL_fIL[dim]);
				dot_with_normal(NfnI,Neq*Nvar,n_fIL,dFluxViscNumdQR_fIL[dim],dnFluxViscNumdQR_fIL[dim]);
			}

			for (size_t i = 0; i < NfnI*Neq*Nvar; i++) {
				dnFluxViscNumdWL_fIL[i] *= 0.5;
				dnFluxViscNumdWR_fIL[i] *= 0.5;
			}

			for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NfnI*Neq*Nvar; i++) {
				dnFluxViscNumdQL_fIL[dim][i] *= 0.5;
				dnFluxViscNumdQR_fIL[dim][i] *= 0.5;
			}}
		}
	}

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
		if (!((BC % BC_STEP_SC == BC_NOSLIP_T) ||
		      (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC))) {
			printf("%d\n",BC);
			EXIT_UNSUPPORTED;
			// Update the boundary treatment if the assumptions on the treatment for the linearization with respect to
			// Q are not valid.
		}
	}

	// Include the BC information in dnFluxViscNumdWL_fIL if on a boundary
	if (imex_type == 'I' && Boundary) {
		// This is only done for dWBdW, the contribution from dQd(W/Q) is assumed to have been properly treated above.
		unsigned int const BC = FACE->BC;

		double const *const XYZ_fIL    = FACE->XYZ_fI;
		double       *const dWBdWL_fIL = malloc(NfnI*Nvar*Nvar * sizeof *dWBdWL_fIL); // free

		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;

		BCdata->XYZ    = XYZ_fIL;
		BCdata->nL     = n_fIL;
		BCdata->WL     = WL_fIL;
		BCdata->GradWL = NULL;

		BCdata->dWBdWL = dWBdWL_fIL;

		if (BC % BC_STEP_SC == BC_NOSLIP_T)
			jacobian_boundary_NoSlip_Dirichlet(BCdata);
		else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC)
			jacobian_boundary_NoSlip_Adiabatic(BCdata);
		else
			EXIT_UNSUPPORTED;
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
	                   Neq      = d+2,
	                   Nvar     = d+2,
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
				eqMax  = d+2;
				varMax = d+2;
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

	unsigned int const d              = DB.d,
	                   Neq            = d+2,
	                   Nvar           = d+2,
	                   IndFType       = FDATA->IndFType,
	                   NfnI           = OPS[IndFType]->NfnI,
	                   *const nOrdLR  = OPS[IndFType]->nOrdLR;

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
	unsigned int const d    = DB.d,
	                   Neq  = d+2,
	                   Nvar = d+2;

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
			RHS   = FDATA->FACE->RHSIn;
		} else if (side == 'R') {
			FDATA = FDATAR;
			RHS   = FDATA->FACE->RHSOut;

			swap_FACE_orientation(FDATAR,imex_type,coef_type);
		} else {
			EXIT_UNSUPPORTED;
		}

		struct S_OPERATORS_F const *const *const OPS = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

		unsigned int const d            = DB.d,
		                   Neq          = d+2,
		                   Nvar         = d+2,
		                   P            = FDATA->P,
		                   Eclass       = FDATA->Eclass,
		                   Vf           = FDATA->Vf,
		                   f            = FDATA->f,
		                   SpOp         = FDATA->SpOp,
		                   IndFType     = FDATA->IndFType,
		                   NfnI         = OPS[IndFType]->NfnI,
		                   NvnS         = OPS[0]->NvnS,
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

		unsigned int const VfL = FDATAL->Vf,
		                   VfR = FDATAR->Vf,
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
			compute_LHS_FACE_Inviscid_Weak(NvnSL,NvnSL,NfnI,I_FF,nANumL_fI,ChiSL_fIL,IdnFdW,FACE->LHSInIn);

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

			compute_LHS_FACE_Inviscid_Weak(NvnSL,NvnSR,NfnI,I_FF,nANumR_fI,ChiSR_fIL,IdnFdW,FACE->LHSOutIn);

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

			compute_LHS_FACE_Inviscid_Weak(NvnSR,NvnSL,NfnI,I_FF,nANumL_fI,ChiSL_fIR,IdnFdW,FACE->LHSInOut);

			free(ChiSL_fIR);

			// LHSRR (Effect of (R)ight VOLUME on (R)ight VOLUME)
			double const *const ChiSR_fIR = OPSR[0]->ChiS_fI[VfR];

			compute_LHS_FACE_Inviscid_Weak(NvnSR,NvnSR,NfnI,I_FF,nANumR_fI,ChiSR_fIR,IdnFdW,FACE->LHSOutOut);

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
		eqMax  = d+2;
		varMax = d+2;
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
                         char const imex_type, char const FORM_MF1)
{
	/*
	 *	Purpose:
	 *		Finalize the RHS/LHS term contributions to Qhat from FACEs.
	 *
	 *	Comments:
	 *		Except for the "-ve" sign, the need to compute d terms instead of 1, and the storage of the use of different
	 *		arrays to be operated on and stored in, this function is identical to finalize_FACE_Inviscid_Weak.
	 *
	 *		When FORM_MF1 == 'S', the local contribution is added to the numerical flux before evaluating QhatF.
	 */

	if (!(imex_type == 'E' || imex_type == 'I'))
		EXIT_UNSUPPORTED;

	struct S_FACE const *const FACE = FDATAL->FACE;

	if (imex_type == 'E') {
		if (FORM_MF1 == 'S')
			correct_numerical_solution_strong(FDATAL,'E',side,FORM_MF1);

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

		unsigned int const d            = DB.d,
		                   Neq          = d+2,
		                   Nvar         = d+2,
		                   P            = FDATA->P,
		                   Eclass       = FDATA->Eclass,
		                   Vf           = FDATA->Vf,
		                   f            = FDATA->f,
		                   SpOp         = FDATA->SpOp,
		                   IndFType     = FDATA->IndFType,
		                   NfnI         = OPS[IndFType]->NfnI,
		                   NvnS         = OPS[0]->NvnS,
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
			if (FORM_MF1 == 'S')
				correct_numerical_solution_strong(FDATAL,'I',side,FORM_MF1);

			unsigned int const NvnSL = OPSL[0]->NvnS;

			double const *const *const dnSolNumdWL_fI = (double const *const *const) FDATAL->NFluxData->dnSolNumdWL_fI;

			// Qhat_WhatLL (Effect of (L)eft VOLUME on (L)eft VOLUME)
			I_FF   = OPSL[0]->I_Weak_FF[VfL];
			IdnFdW = malloc(NvnSL*NfnI * sizeof *IdnFdW); // free

			double const *const ChiSL_fIL = OPSL[0]->ChiS_fI[VfL];
			compute_LHS_QhatF_Weak(NvnSL,NvnSL,NfnI,I_FF,dnSolNumdWL_fI,ChiSL_fIL,IdnFdW,FACE->Qhat_WhatLL,Boundary);

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

			// Qhat_WhatRL (Effect of (R)ight VOLUME on (L)eft VOLUME)
			double *const ChiSR_fIL = malloc(NvnSR*NfnI * sizeof *ChiSR_fIL); // free

			ChiS_fI = OPSR[0]->ChiS_fI[VfR];
			for (size_t i = 0; i < NfnI; i++) {
			for (size_t j = 0; j < NvnSR; j++) {
				ChiSR_fIL[i*NvnSR+j] = ChiS_fI[nOrdRL[i]*NvnSR+j];
			}}

			compute_LHS_QhatF_Weak(NvnSL,NvnSR,NfnI,I_FF,dnSolNumdWR_fI,ChiSR_fIL,IdnFdW,FACE->Qhat_WhatRL,Boundary);

			free(ChiSR_fIL);
			free(IdnFdW);

			// Swap orientation of numerical flux Jacobian terms
			if (FORM_MF1 == 'S')
				correct_numerical_solution_strong(FDATAR,'I',side,FORM_MF1);
			swap_FACE_orientation(FDATAR,'I','Q');

			I_FF   = OPSR[0]->I_Weak_FF[VfR];
			IdnFdW = malloc(NvnSR*NfnI * sizeof *IdnFdW); // free

			// Qhat_WhatLR (Effect of (L)eft VOLUME on (R)ight VOLUME)
			double *const ChiSL_fIR = malloc(NvnSL*NfnI * sizeof *ChiSL_fIR); // free

			ChiS_fI = OPSL[0]->ChiS_fI[VfL];
			for (size_t i = 0; i < NfnI; i++) {
			for (size_t j = 0; j < NvnSL; j++) {
				ChiSL_fIR[i*NvnSL+j] = ChiS_fI[nOrdLR[i]*NvnSL+j];
			}}

			compute_LHS_QhatF_Weak(NvnSR,NvnSL,NfnI,I_FF,dnSolNumdWL_fI,ChiSL_fIR,IdnFdW,FACE->Qhat_WhatLR,Boundary);

			free(ChiSL_fIR);

			// Qhat_WhatRR (Effect of (R)ight VOLUME on (R)ight VOLUME)
			double const *const ChiSR_fIR = OPSR[0]->ChiS_fI[VfR];

			compute_LHS_QhatF_Weak(NvnSR,NvnSR,NfnI,I_FF,dnSolNumdWR_fI,ChiSR_fIR,IdnFdW,FACE->Qhat_WhatRR,Boundary);

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
	                   Neq  = d+2,
	                   Nvar = d+2;

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
	 *		right operator is given by Q_What (as opposed to ChiS_fI).
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

		// Q_WhatLL (Effect of (L)eft VOLUME on (L)eft VOLUME)
		I_FF   = OPSL[0]->I_Weak_FF[VfL];
		IdnFdQ = malloc(NvnSL*NfnI * sizeof *IdnFdQ); // free

		double const *const *const Q_What = (double const *const *const) FDATAL->Q_WhatLL;
		compute_LHS_FACE_Q_Weak(NvnSL,NvnSL,NfnI,I_FF,dnFluxViscNumdQL_fI,Q_What,IdnFdQ,FACE->LHSInIn,Boundary);
		array_free2_d(d,(double **) Q_What);
		free(IdnFdQ);
	} else if (side == 'R') {
		if (Boundary)
			EXIT_UNSUPPORTED;

		double const *const *      Q_What;
		double const *const *const dnFluxViscNumdQL_fI = (double const *const *const) NFluxData->dnFluxViscNumdQL_fI;
		double const *const *const dnFluxViscNumdQR_fI = (double const *const *const) NFluxData->dnFluxViscNumdQR_fI;

		unsigned int const NvnSL = OPSL[0]->NvnS,
		                   NvnSR = OPSR[0]->NvnS;

		I_FF   = OPSL[0]->I_Weak_FF[VfL];
		IdnFdQ = malloc(NvnSL*NfnI * sizeof *IdnFdQ); // free

		// Q_WhatRL (Effect of (R)ight VOLUME on (L)eft VOLUME)
		Q_What = (double const *const *const) FDATAL->Q_WhatRL;
		compute_LHS_FACE_Q_Weak(NvnSL,NvnSR,NfnI,I_FF,dnFluxViscNumdQR_fI,Q_What,IdnFdQ,FACE->LHSOutIn,Boundary);
		array_free2_d(d,(double **) Q_What);
		free(IdnFdQ);

		// Swap orientation of numerical flux Jacobian terms
		swap_FACE_orientation(FDATAR,'I','P');

		I_FF   = OPSR[0]->I_Weak_FF[VfR];
		IdnFdQ = malloc(NvnSR*NfnI * sizeof *IdnFdQ); // free

		// Q_WhatLR (Effect of (L)eft VOLUME on (R)ight VOLUME)
		Q_What = (double const *const *const) FDATAR->Q_WhatLR;
		compute_LHS_FACE_Q_Weak(NvnSR,NvnSL,NfnI,I_FF,dnFluxViscNumdQL_fI,Q_What,IdnFdQ,FACE->LHSInOut,Boundary);
		array_free2_d(d,(double **) Q_What);

		// Q_WhatRR (Effect of (R)ight VOLUME on (R)ight VOLUME)
		Q_What = (double const *const *const) FDATAR->Q_WhatRR;
		compute_LHS_FACE_Q_Weak(NvnSR,NvnSR,NfnI,I_FF,dnFluxViscNumdQR_fI,Q_What,IdnFdQ,FACE->LHSOutOut,Boundary);
		array_free2_d(d,(double **) Q_What);
		free(IdnFdQ);
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_LHSQF_Weak(struct S_FACE *const FACE)
{
	if (!DB.Viscous)
		return;

	unsigned int const d    = DB.d,
	                   Nvar = d+2,
	                   Neq  = d+2;

	struct S_VOLUME *const VL = FACE->VIn,
	                *const VR = FACE->VOut;

	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t varQ = 0; varQ < Nvar; varQ++) {
	for (size_t var = 0; var < Nvar; var++) {
		unsigned int const NvnSL = VL->NvnS;
		if (FACE->Boundary) {
			size_t const IndLHSQ = (eq*Nvar+varQ)*NvnSL*NvnSL,
			             IndQ_W  = (varQ*Nvar+var)*NvnSL*NvnSL,
			             IndLHS  = (eq*Nvar+var)*NvnSL*NvnSL;
			for (size_t dim = 0; dim < d; dim++)
				mm_d(CBRM,CBNT,CBNT,NvnSL,NvnSL,NvnSL,1.0,1.0,&VL->LHSQ[dim][IndLHSQ],&FACE->Qhat_WhatLL[dim][IndQ_W],&VL->LHS[IndLHS]);
		} else {
			// dQhat/dWhat is block diagonal for this term
			if (var != varQ)
				continue;

			size_t const Indeqvar = (eq*Nvar+var);
			size_t IndLHS, IndLHSQ;

			unsigned int const NvnSR = VR->NvnS;
			for (size_t dim = 0; dim < d; dim++) {
				IndLHS = Indeqvar*NvnSL*NvnSL;
				mm_d(CBRM,CBNT,CBNT,NvnSL,NvnSL,NvnSL,1.0,1.0,&VL->LHSQ[dim][IndLHS],FACE->Qhat_WhatLL[dim],&VL->LHS[IndLHS]);

				IndLHS = Indeqvar*NvnSL*NvnSR;
				IndLHSQ = Indeqvar*NvnSL*NvnSL;
				mm_d(CBRM,CBNT,CBNT,NvnSL,NvnSR,NvnSL,1.0,1.0,&VL->LHSQ[dim][IndLHSQ],FACE->Qhat_WhatRL[dim],&FACE->LHSOutIn[IndLHS]);

				IndLHSQ = Indeqvar*NvnSR*NvnSR;
				mm_d(CBRM,CBNT,CBNT,NvnSR,NvnSL,NvnSR,1.0,1.0,&VR->LHSQ[dim][IndLHSQ],FACE->Qhat_WhatLR[dim],&FACE->LHSInOut[IndLHS]);

				IndLHS = Indeqvar*NvnSR*NvnSR;
				mm_d(CBRM,CBNT,CBNT,NvnSR,NvnSR,NvnSR,1.0,1.0,&VR->LHSQ[dim][IndLHS],FACE->Qhat_WhatRR[dim],&VR->LHS[IndLHS]);
			}
		}
	}}}
}
