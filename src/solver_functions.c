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
