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

#include "array_swap.h"
#include "exact_solutions.h"
#include "variable_functions.h"
#include "boundary_conditions.h"
#include "jacobian_boundary_conditions.h"
#include "fluxes_inviscid.h"
#include "jacobian_fluxes_inviscid.h"

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
// VOLUME functions
// **************************************************************************************************** //

void init_ops_VOLUME(struct S_OPERATORS_V *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	// Initialize DB Parameters
	unsigned int ***SF_BE = DB.SF_BE;

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
		OPS->D_Weak  = ELEMENT->Ds_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT->Is_Weak_VV[P][P][0];

		OPS->ChiS_vI_SF = ELEMENT_SF->ChiS_vIs[P][P][0];
		OPS->D_Weak_SF  = ELEMENT_SF->Ds_Weak_VV[P][P][0];
		OPS->I_Weak_SF  = ELEMENT_SF->Is_Weak_VV[P][P][0];

		OPS->D_Weak_sp = ELEMENT->Ds_Weak_VV_sp[P][P][0];
	} else {
		OPS->NvnI    = ELEMENT->NvnIc[P];
		OPS->NvnI_SF = ELEMENT_SF->NvnIc[P];

		OPS->ChiS_vI    = ELEMENT->ChiS_vIc[P][P][0];
		OPS->D_Weak  = ELEMENT->Dc_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT->Ic_Weak_VV[P][P][0];

		OPS->ChiS_vI_SF = ELEMENT_SF->ChiS_vIc[P][P][0];
		OPS->D_Weak_SF  = ELEMENT_SF->Dc_Weak_VV[P][P][0];
		OPS->I_Weak_SF  = ELEMENT_SF->Ic_Weak_VV[P][P][0];

		OPS->D_Weak_sp = ELEMENT->Dc_Weak_VV_sp[P][P][0];
	}
}

void init_VDATA(struct S_VDATA *VDATA, struct S_VOLUME *VOLUME)
{
	VDATA->VOLUME = VOLUME;

	VDATA->P      = VOLUME->P;
	VDATA->Eclass = VOLUME->Eclass;
	init_ops_VOLUME(VDATA->OPS[0],VOLUME,0);
	if (VOLUME->type == WEDGE)
		init_ops_VOLUME(VDATA->OPS[1],VOLUME,1);
}

void compute_W_vI(struct S_VDATA *VDATA, double *W_vI)
{
	/*
	 *	Purpose:
	 *		Interpolate VOLUME What coefficients to VOLUME cubature nodes.
	 *
	 *	Comments:
	 *		Various options are available for computing W_vI:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
	 */

	unsigned int Collocated = DB.Collocated;

	struct S_VOLUME *VOLUME = VDATA->VOLUME;

	if (Collocated) {
		EXIT_UNSUPPORTED; // Set W_vI = VOLUME->What in calling function.
	} else {
		unsigned int d        = DB.d,
		             Nvar     = DB.Nvar,
		             ***SF_BE = DB.SF_BE;

		unsigned int P      = VDATA->P,
		             Eclass = VDATA->Eclass;

		struct S_OPERATORS_V **OPS = VDATA->OPS;

		unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX];
		double       *OP[DMAX];
		if (Eclass == C_TP && SF_BE[P][0][0]) {
			get_sf_parameters(OPS[0]->NvnS_SF,OPS[0]->NvnI_SF,OPS[0]->ChiS_vI_SF,0,0,NULL,NIn,NOut,OP,d,3,Eclass);

			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 0;

			sf_apply_d(VOLUME->What,W_vI,NIn,NOut,Nvar,OP,Diag,d);
		} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
			get_sf_parameters(OPS[0]->NvnS_SF,OPS[0]->NvnI_SF,OPS[0]->ChiS_vI_SF,
			                  OPS[1]->NvnS_SF,OPS[1]->NvnI_SF,OPS[1]->ChiS_vI_SF,NIn,NOut,OP,d,3,Eclass);

			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 0;
			Diag[1] = 2; // TRI operator (dim = 0) takes care of dim = 1 as well.

			sf_apply_d(VOLUME->What,W_vI,NIn,NOut,Nvar,OP,Diag,d);
		} else {
			mm_CTN_d(OPS[0]->NvnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_vI,VOLUME->What,W_vI);
		}
	}
}

void convert_between_rp(const unsigned int Nn, const unsigned int Nrc, const double *C, double *Ap, double *Ar,
                        const char *conv_type)
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

	unsigned int d = DB.d;

	if (strstr(conv_type,"FluxToRef")) {
		memset(Ar,0.0,Nn*Nrc*d * sizeof *Ar);
		unsigned int IndAr, IndAp, IndC;

		for (size_t col = 0; col < Nrc; col++) {
		for (size_t dim1 = 0; dim1 < d; dim1++) {
			IndAr = (Nrc*dim1+col)*Nn;
			for (size_t dim2 = 0; dim2 < d; dim2++) {
				IndAp = (col*d+dim2)*Nn;
				IndC  = (dim1*d+dim2)*Nn;
				for (size_t n = 0; n < Nn; n++)
					Ar[IndAr+n] += Ap[IndAp+n]*C[IndC+n];
			}
		}}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Inviscid_Weak(const unsigned int Nrc, const double *Ar_vI, double *RLHS, const char *term_type,
                                   struct S_VDATA *VDATA)
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

	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
	             ***SF_BE   = DB.SF_BE;

	struct S_OPERATORS_V **OPS = VDATA->OPS;

	unsigned int P      = VDATA->P,
	             Eclass = VDATA->Eclass;

	unsigned int NvnS = OPS[0]->NvnS,
	             NvnI = OPS[0]->NvnI;

	if (strstr(term_type,"RHS")) {
		memset(RLHS,0.0,NvnS*Nrc * sizeof *RLHS);
		if (Eclass == C_TP && SF_BE[P][0][0]) {
			unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX];
			double       *OP[DMAX];

			double *DAr = malloc(NvnS*Nrc * sizeof *DAr); // free
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
			double       *OP[DMAX], *OP0, *OP1;

			double *DAr = malloc(NvnS*Nrc * sizeof *DAr); // free
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
			double *DAr = malloc(NvnS*Nrc * sizeof *DAr); // free
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
	} else if (strstr(term_type,"LHS")) {
		unsigned int Neq  = DB.Neq,
		             Nvar = DB.Nvar;

		const double *Ar_vI_ptr[d];
		for (size_t dim = 0; dim < d; dim++)
			Ar_vI_ptr[dim] = &Ar_vI[Nvar*Neq*NvnI*dim];

		if (Collocated) {
			if (Eclass == C_TP || Eclass == C_WEDGE) {
				struct S_OpCSR **D = OPS[0]->D_Weak_sp;

				memset(RLHS,0.0,NvnS*Nrc * sizeof *RLHS);
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					unsigned int IndAr  = (eq*Nvar+var)*NvnI,
								 IndLHS = (eq*Nvar+var)*NvnS*NvnS;

					for (size_t dim = 0; dim < d; dim++) {
						unsigned int *rowIndex = D[dim]->rowIndex,
						             *columns  = D[dim]->columns;
						double       *values   = D[dim]->values;

						for (size_t i = 0; i < NvnS; i++) {
							unsigned int IndD = i*NvnI;
							for (size_t Indj = *rowIndex, IndjMax = *(++rowIndex); Indj < IndjMax; Indj++) {
								RLHS[IndLHS+IndD+(*columns)] += (*values++)*Ar_vI_ptr[dim][IndAr+(*columns)];
								columns++;
							}
						}
					}
				}}
			} else {
				double **D = OPS[0]->D_Weak;

				memset(RLHS,0.0,NvnS*Nrc * sizeof *RLHS);
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					unsigned int IndAr  = (eq*Nvar+var)*NvnI,
								 IndLHS = (eq*Nvar+var)*NvnS*NvnS;

					for (size_t dim = 0; dim < d; dim++) {
						for (size_t i = 0; i < NvnS; i++) {
							unsigned int IndD = i*NvnI;
							for (size_t j = 0; j < NvnI; j++)
								RLHS[IndLHS+IndD+j] += D[dim][IndD+j]*Ar_vI_ptr[dim][IndAr+j];
						}
					}
				}}
			}
		} else {
			double **D = OPS[0]->D_Weak;

			double *DAr_vI = malloc(NvnS*NvnI * sizeof *DAr_vI); // free
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				unsigned int IndAr = (eq*Nvar+var)*NvnI;

				memset(DAr_vI,0.0,NvnS*NvnI * sizeof *DAr_vI);
				for (size_t dim = 0; dim < d; dim++) {
					for (size_t i = 0; i < NvnS; i++) {
						unsigned int IndD = i*NvnI;
						for (size_t j = 0; j < NvnI; j++)
							DAr_vI[IndD+j] += D[dim][IndD+j]*Ar_vI_ptr[dim][IndAr+j];
					}
				}

				unsigned int IndLHS = (eq*Nvar+var)*NvnS*NvnS;
				mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,0.0,DAr_vI,OPS[0]->ChiS_vI,&RLHS[IndLHS]);
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
						const double *Ar_vI_ptr2 = &Ar_vI_ptr[dim][(eq*Nvar+var)*NvnI+i];
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


// **************************************************************************************************** //
// FACE functions
// **************************************************************************************************** //

// Remove what is not necessary in init_ops_FACE (ToBeDeleted)
void init_ops_FACE(struct S_OPERATORS_F *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                   const unsigned int IndClass)
{
	// Initialize DB Parameters
	unsigned int ***SF_BE = DB.SF_BE;

	// Standard datatypes
	unsigned int PV, PF, Vtype, Eclass, FtypeInt, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS, *ELEMENT_FACE;

	// silence
	ELEMENT_OPS = NULL;

	PV       = VOLUME->P;
	PF       = FACE->P;
	Vtype    = VOLUME->type;
	Eclass   = VOLUME->Eclass;

	FtypeInt    = FACE->typeInt;
	IndOrdInOut = FACE->IndOrdInOut;
	IndOrdOutIn = FACE->IndOrdOutIn;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACE = get_ELEMENT_FACE(Vtype,IndClass);
	if ((Eclass == C_TP && SF_BE[PF][0][1]) || (Eclass == C_WEDGE && SF_BE[PF][1][1]))
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else
		ELEMENT_OPS = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[PV];
	OPS->NvnS_SF = ELEMENT_OPS->NvnS[PV];

	OPS->NfnS      = ELEMENT_FACE->NvnS[PF];
	OPS->ChiS_fS   = ELEMENT->ChiS_fS[PV][PF];
	if (FtypeInt == 's') {
		// Straight FACE Integration
		OPS->NfnI    = ELEMENT->NfnIs[PF][IndClass];
		OPS->NfnI_SF = ELEMENT_OPS->NfnIs[PF][0];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIs[PF];

		OPS->ChiS_fI   = ELEMENT_OPS->ChiS_fIs[PV][PF];
		OPS->ChiS_vI   = ELEMENT_OPS->ChiS_vIs[PV][PF];
		OPS->I_Weak_FF = ELEMENT_OPS->Is_Weak_FF[PV][PF];
		OPS->I_Weak_VV = ELEMENT_OPS->Is_Weak_VV[PV][PF];

		OPS->ChiS_fI_sp   = ELEMENT->ChiS_fIs_sp[PV][PF];
		OPS->I_Weak_FF_sp = ELEMENT->Is_Weak_FF_sp[PV][PF];

		OPS->GfS_fI    = ELEMENT->GfS_fIs[PF][PF];

		OPS->nOrdInOut = ELEMENT_FACE->nOrd_fIs[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fIs[PF][IndOrdOutIn];
	} else {
		// Curved FACE Integration
		OPS->NfnI    = ELEMENT->NfnIc[PF][IndClass];
		OPS->NfnI_SF = ELEMENT_OPS->NfnIc[PF][0];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIc[PF];

		OPS->ChiS_fI   = ELEMENT_OPS->ChiS_fIc[PV][PF];
		OPS->ChiS_vI   = ELEMENT_OPS->ChiS_vIc[PV][PF];
		OPS->I_Weak_FF = ELEMENT_OPS->Ic_Weak_FF[PV][PF];
		OPS->I_Weak_VV = ELEMENT_OPS->Ic_Weak_VV[PV][PF];

		OPS->ChiS_fI_sp   = ELEMENT->ChiS_fIc_sp[PV][PF];
		OPS->I_Weak_FF_sp = ELEMENT->Ic_Weak_FF_sp[PV][PF];

		OPS->GfS_fI    = ELEMENT->GfS_fIc[PF][PF];

		OPS->nOrdInOut = ELEMENT_FACE->nOrd_fIc[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fIc[PF][IndOrdOutIn];
	}
}

void init_FDATA(struct S_FDATA *FDATA, struct S_FACE *FACE, const char side)
{
	// Initialize DB Parameters
	unsigned int Collocated = DB.Collocated;

	FDATA->FACE = FACE;

	FDATA->P = FACE->P;

	if (side == 'L') {
		FDATA->VOLUME = FACE->VIn;
		FDATA->Vf     = FACE->VfIn;
	} else if (side == 'R') {
		FDATA->VOLUME = FACE->VOut;
		FDATA->Vf     = FACE->VfOut;
	} else {
		EXIT_UNSUPPORTED;
	}

	FDATA->f    = (FDATA->Vf)/NFREFMAX;
	FDATA->SpOp = Collocated && ((FDATA->Vf) % NFREFMAX == 0 && FDATA->VOLUME->P == FDATA->P);

	FDATA->Eclass = FDATA->VOLUME->Eclass;
	FDATA->IndFType = get_IndFType(FDATA->Eclass,FDATA->f);

	init_ops_FACE(FDATA->OPS[0],FDATA->VOLUME,FACE,0);
	if (FDATA->VOLUME->type == WEDGE || FDATA->VOLUME->type == PYR)
		// Needed for sum factorized operators and alternate FACE operators (TRIs/QUADs)
		init_ops_FACE(FDATA->OPS[1],FDATA->VOLUME,FACE,1);
}

void compute_W_fI(const struct S_FDATA *FDATA, double *W_fI)
{
	/*
	 *	Purpose:
	 *		Interpolate VOLUME What coefficients to FACE cubature nodes.
	 *
	 *	Comments:
	 *		Various options are available for computing W_fI:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
	 */

	unsigned int d            = DB.d,
	             Nvar         = DB.Nvar,
	             *VFPartUnity = DB.VFPartUnity,
	             ***SF_BE     = DB.SF_BE;

	struct S_OPERATORS_F **OPS   = FDATA->OPS;
	struct S_VOLUME      *VOLUME = FDATA->VOLUME;

	unsigned int P      = FDATA->P,
	             Eclass = FDATA->Eclass,
	             Vf     = FDATA->Vf,
	             f      = FDATA->f,
	             SpOp   = FDATA->SpOp;

	unsigned int NIn[DMAX], NOut[DMAX], Diag[DMAX], NOut0, NOut1;
	double       *OP[DMAX], **OP0, **OP1;

	if (Eclass == C_TP && SF_BE[P][0][1]) {
		get_sf_parametersF(OPS[0]->NvnS_SF,OPS[0]->NvnI_SF,OPS[0]->ChiS_vI,
		                   OPS[0]->NvnS_SF,OPS[0]->NfnI_SF,OPS[0]->ChiS_fI,NIn,NOut,OP,d,Vf,C_TP);

		if (SpOp) {
			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 2;
			Diag[f/2] = 0;
		} else {
			for (size_t dim = 0; dim < d; dim++)
				Diag[dim] = 0;
		}

		sf_apply_d(VOLUME->What,W_fI,NIn,NOut,Nvar,OP,Diag,d);
	} else if (Eclass == C_WEDGE && SF_BE[P][1][1]) {
		if (f < 3) { OP0   = OPS[0]->ChiS_fI, OP1   = OPS[1]->ChiS_vI;
		             NOut0 = OPS[0]->NfnI_SF, NOut1 = OPS[1]->NvnI_SF;
		} else {     OP0   = OPS[0]->ChiS_vI, OP1   = OPS[1]->ChiS_fI;
		             NOut0 = OPS[0]->NvnI_SF, NOut1 = OPS[1]->NfnI_SF; }
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

		sf_apply_d(VOLUME->What,W_fI,NIn,NOut,Nvar,OP,Diag,d);
	} else if (((SpOp && (Eclass == C_TP || Eclass == C_WEDGE)) || (VFPartUnity[Eclass])) && DB.AllowSparseFACE) {
		mm_CTN_CSR_d(OPS[0]->NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI_sp[Vf],VOLUME->What,W_fI);
	} else {
		mm_CTN_d(OPS[0]->NfnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_fI[Vf],VOLUME->What,W_fI);
	}
}

void compute_WR_fIL(const struct S_FDATA *FDATA, const double *WL_fIL, double *WR_fIL)
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

	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	struct S_OPERATORS_F **OPS = FDATA->OPS;
	struct S_FACE        *FACE = FDATA->FACE;

	unsigned int IndFType = FDATA->IndFType;

	unsigned int BC      = FACE->BC,
	             NfnI    = OPS[IndFType]->NfnI,
	             *nOrdRL = OPS[IndFType]->nOrdOutIn;

	double       *XYZ_fIL = FACE->XYZ_fI,
	             *n_fIL   = FACE->n_fI;

	if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACE
		compute_W_fI(FDATA,WR_fIL);
		array_rearrange_d(NfnI,Nvar,nOrdRL,'C',WR_fIL);
	} else { // Boundary FACE
		if (BC % BC_STEP_SC == BC_RIEMANN) {
			boundary_Riemann(NfnI,1,XYZ_fIL,WL_fIL,NULL,WR_fIL,n_fIL,d);
		} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
			bool ExactSlipWall = 0;
			if (ExactSlipWall) {
				double *UR_fIL = malloc(NVAR3D*NfnI * sizeof *UR_fIL); // free

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
		} else {
			EXIT_UNSUPPORTED;
		}
	}
}

void compute_numerical_flux(const struct S_FDATA *FDATA, struct S_NumericalFlux *NFluxData, const char imex_type)
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

	unsigned int d    = DB.d,
	             Neq  = DB.Neq,
	             Nvar = DB.Nvar;

	const double *WL_fIL = NFluxData->WL_fIL,
	             *WR_fIL = NFluxData->WR_fIL;

	double       *nFluxNum_fIL     = NFluxData->nFluxNum_fIL,
	             *dnFluxNumdWL_fIL = NFluxData->dnFluxNumdWL_fIL,
	             *dnFluxNumdWR_fIL = NFluxData->dnFluxNumdWR_fIL;

	struct S_OPERATORS_F **OPS = FDATA->OPS;
	struct S_FACE        *FACE = FDATA->FACE;

	unsigned int IndFType = FDATA->IndFType;

	unsigned int Boundary = FACE->Boundary,
	             NfnI     = OPS[IndFType]->NfnI;

	double       *n_fIL = FACE->n_fI;

	switch (DB.InviscidFluxType) {
	case FLUX_LF:
		flux_LF(NfnI,1,WL_fIL,WR_fIL,nFluxNum_fIL,n_fIL,d,Neq);
		if (imex_type == 'I') {
			jacobian_flux_LF(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWL_fIL,n_fIL,d,Neq,'L');
			if (!Boundary)
				jacobian_flux_LF(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWR_fIL,n_fIL,d,Neq,'R');
		}
		break;
	case FLUX_ROE:
		flux_Roe(NfnI,1,WL_fIL,WR_fIL,nFluxNum_fIL,n_fIL,d,Neq);
		if (imex_type == 'I') {
			jacobian_flux_Roe(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWL_fIL,n_fIL,d,Neq,'L');
			if (!Boundary)
				jacobian_flux_Roe(NfnI,1,WL_fIL,WR_fIL,dnFluxNumdWR_fIL,n_fIL,d,Neq,'R');
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	// Include the BC information in dnFluxNumWL_fIL if on a boundary
	if (Boundary) {
		unsigned int BC       = FACE->BC;
		double       *XYZ_fIL = FACE->XYZ_fI;

		dWRdWL_fIL = malloc(NfnI*Nvar*Nvar * sizeof *dWRdWL_fIL); // free
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
		else
			EXIT_UNSUPPORTED;

		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			InddnFdWL = (eq*Nvar+var)*NfnI;

			for (i = 0; i < Nvar; i++) {
				InddnFdWR = (eq*Neq+i)*NfnI;
				InddWRdWL = (var*Nvar+i)*NfnI;
				for (n = 0; n < NfnI; n++)
					dnFluxNumdWL_fIL[InddnFdWL+n] += dnFluxNumdWR_fIL[InddnFdWR+n]*dWRdWL_fIL[InddWRdWL+n];
			}
		}}

		free(dWRdWL_fIL);
	}
}
