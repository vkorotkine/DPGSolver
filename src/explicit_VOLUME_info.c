// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "explicit_VOLUME_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_OpCSR.h"

#include "element_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h"
#include "fluxes_inviscid.h"
#include "array_print.h"

#undef I // No complex variables used here

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS term.
 *
 *	Comments:
 *		Much of the explicit_VOLUME_info and implicit_VOLUME_info functions are identical, consider combining them
 *		(ToBeDeleted).
 *		Certain multiplications can be avoided when computing either Fr from F or when computing RHS terms based on the
 *		sparsity of the flux Jacobian. Test performance improvement if these terms are neglected (ToBeModified).
 *		Vectorization does not improve performance based on preliminary testing. This is likely a result of the large
 *		memory allocation/deallocation overhead which is required for the adaptive code. For non-adaptive versions of
 *		the code, this flexibility is not required and comparison may yield favourable results for the vectorized code.
 *		However, the ultimate goal of the code is only to run in the adaptive setting, thus this may not be worth
 *		pursuing. (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnI, NvnS, NvnS_SF, NvnI_SF;
	double       *ChiS_vI, **D_Weak, *I_Weak;

	struct S_OpCSR **D_Weak_sp;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass);
static void compute_VOLUME_RHS_EFE(void);
static void compute_VOLUMEVec_RHS_EFE(void);

void explicit_VOLUME_info(void)
{
	// Initialize DB Parameters
	unsigned int EFE        = DB.EFE,
	             Vectorized = DB.Vectorized;

	if (EFE) {
		switch (Vectorized) {
		case 0:
			compute_VOLUME_RHS_EFE();
			break;
		default:
			compute_VOLUMEVec_RHS_EFE();
			break;
		}
	} else {
		;
	}
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	// Initialize DB Parameters
	unsigned int ***SF_BE = DB.SF_BE;

	// Standard datatypes
	unsigned int P, type, curved, Eclass;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	// silence
	ELEMENT_OPS = NULL;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;
	Eclass = VOLUME->Eclass;

	ELEMENT = get_ELEMENT_type(type);
	if ((Eclass == C_TP && SF_BE[P][0][0]) || (Eclass == C_WEDGE && SF_BE[P][1][0]))
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else
		ELEMENT_OPS = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[P];
	OPS->NvnS_SF = ELEMENT_OPS->NvnS[P];
	if (!curved) {
		OPS->NvnI    = ELEMENT->NvnIs[P];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIs[P];

		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIs[P][P][0];
		OPS->D_Weak  = ELEMENT_OPS->Ds_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT_OPS->Is_Weak_VV[P][P][0];

		OPS->D_Weak_sp = ELEMENT->Ds_Weak_VV_sp[P][P][0];
	} else {
		OPS->NvnI    = ELEMENT->NvnIc[P];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIc[P];

		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIc[P][P][0];
		OPS->D_Weak  = ELEMENT_OPS->Dc_Weak_VV[P][P][0];
		OPS->I_Weak  = ELEMENT_OPS->Ic_Weak_VV[P][P][0];

		OPS->D_Weak_sp = ELEMENT->Dc_Weak_VV_sp[P][P][0];
	}
}

static void compute_VOLUME_RHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form = DB.Form;
	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
				 Nvar       = DB.Nvar,
				 Neq        = DB.Neq,
	             ***SF_BE   = DB.SF_BE;

	// Standard datatypes
	unsigned int i, eq, dim1, dim2, P,
	             IndFr, IndF, IndC, IndRHS, Eclass,
	             NvnI, NvnS, NvnI_SF[2], NvnS_SF[2], NIn[3], NOut[3], Diag[3];
	double       *W_vI, *F_vI, *Fr_vI, *C_vI, *RHS, *DFr, **D, *I, *OP[3], *OP0, *OP1;

	struct S_OpCSR **D_sp;

	struct S_OPERATORS *OPS[2];
	struct S_VOLUME    *VOLUME;

	// silence
	NvnI_SF[1] = 0;
	NvnS_SF[1] = 0;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(Form,"Weak")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			P = VOLUME->P;
//printf("VOLUME: %d\n",VOLUME->indexg);
			// Obtain operators
			init_ops(OPS[0],VOLUME,0);
			if (VOLUME->type == WEDGE)
				init_ops(OPS[1],VOLUME,1);

			Eclass = VOLUME->Eclass;

			// Obtain W_vI
			NvnI = OPS[0]->NvnI;
			if (Collocated) {
				W_vI = VOLUME->What;
			} else {
				W_vI = malloc(NvnI*Nvar * sizeof *W_vI); // free

				if (Eclass == C_TP && SF_BE[P][0][0]) {
					for (i = 0; i < 1; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,0,0,NULL,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;

					sf_apply_d(VOLUME->What,W_vI,NIn,NOut,Nvar,OP,Diag,d);
				} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
					for (i = 0; i < 2; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,
					                  NvnS_SF[1],NvnI_SF[1],OPS[1]->ChiS_vI,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;
					Diag[1] = 2;

					sf_apply_d(VOLUME->What,W_vI,NIn,NOut,Nvar,OP,Diag,d);
				} else {
					mm_CTN_d(NvnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_vI,VOLUME->What,W_vI);
				}
			}
/*
if (VOLUME->indexg == 0) {
printf("eVi: %d\n",VOLUME->indexg);
//array_print_d(OPS[0]->NvnS,Neq,VOLUME->What,'C');
array_print_d(NvnI,Neq,W_vI,'C');
}
*/
			// Compute Flux in reference space
			F_vI = malloc(NvnI*d*Neq * sizeof *F_vI); // free
			flux_inviscid(NvnI,1,W_vI,F_vI,d,Neq);

//array_print_d(NvnI,Neq*d,F_vI,'C');

			if (!Collocated)
				free(W_vI);

			C_vI = VOLUME->C_vI;

			Fr_vI = calloc(NvnI*d*Neq , sizeof *Fr_vI); // free
			for (eq = 0; eq < Neq; eq++) {
			for (dim1 = 0; dim1 < d; dim1++) {
			for (dim2 = 0; dim2 < d; dim2++) {
				IndFr = (eq*d+dim1)*NvnI;
				IndF  = (eq*d+dim2)*NvnI;
				IndC  = (dim1*d+dim2)*NvnI;
				for (i = 0; i < NvnI; i++)
					Fr_vI[IndFr+i] += F_vI[IndF+i]*C_vI[IndC+i];
			}}}
			free(F_vI);

//array_print_d(NvnI,d*d,C_vI,'C');
//array_print_d(NvnI,Neq*d,Fr_vI,'C');

for (eq = 0; eq < Neq; eq++) {
//array_print_d(NvnI,d,&Fr_vI[NvnI*d*eq],'C');
}

			// Compute RHS terms
			NvnS = OPS[0]->NvnS;

			if (VOLUME->RHS)
				free(VOLUME->RHS);
			RHS = calloc(NvnS*Neq , sizeof *RHS); // keep (requires external free)
			VOLUME->RHS = RHS;

			DFr = malloc(NvnS * sizeof *DFr); // free
			if (Eclass == C_TP && SF_BE[P][0][0]) {
				for (i = 0; i < 1; i++) {
					NvnS_SF[i] = OPS[i]->NvnS_SF;
					NvnI_SF[i] = OPS[i]->NvnI_SF;
				}

				I = OPS[0]->I_Weak;
				D = OPS[0]->D_Weak;

				for (dim1 = 0; dim1 < d; dim1++) {
					get_sf_parameters(NvnI_SF[0],NvnS_SF[0],I,NvnI_SF[0],NvnS_SF[0],D[0],NIn,NOut,OP,d,dim1,Eclass);

					if (Collocated) {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 2;
						Diag[dim1] = 0;
					} else {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 0;
					}

					for (eq = 0; eq < Neq; eq++) {
						sf_apply_d(&Fr_vI[(eq*d+dim1)*NvnI],DFr,NIn,NOut,1,OP,Diag,d);

						IndRHS = eq*NvnS;
						for (i = 0; i < NvnS; i++)
							RHS[IndRHS+i] += DFr[i];
					}
				}
			} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
				for (i = 0; i < 2; i++) {
					NvnS_SF[i] = OPS[i]->NvnS_SF;
					NvnI_SF[i] = OPS[i]->NvnI_SF;
				}

				for (dim1 = 0; dim1 < d; dim1++) {
					if (dim1 < 2) OP0 = OPS[0]->D_Weak[dim1], OP1 = OPS[1]->I_Weak;
					else          OP0 = OPS[0]->I_Weak,       OP1 = OPS[1]->D_Weak[0];
					get_sf_parameters(NvnI_SF[0],NvnS_SF[0],OP0,NvnI_SF[1],NvnS_SF[1],OP1,NIn,NOut,OP,d,3,Eclass);

					if (Collocated) {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 2;
						if (dim1 < 2)
							Diag[0] = 0;
						else
							Diag[dim1] = 0;
					} else {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 0;
						Diag[1] = 2;
					}

					for (eq = 0; eq < Neq; eq++) {
						sf_apply_d(&Fr_vI[(eq*d+dim1)*NvnI],DFr,NIn,NOut,1,OP,Diag,d);

						IndRHS = eq*NvnS;
						for (i = 0; i < NvnS; i++)
							RHS[IndRHS+i] += DFr[i];
					}
				}
			} else if (Collocated && (Eclass == C_TP || Eclass == C_WEDGE)) {
				D_sp = OPS[0]->D_Weak_sp;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_CSR_d(NvnS,1,NvnI,D_sp[dim1],&Fr_vI[(eq*d+dim1)*NvnI],DFr);

					IndRHS = eq*NvnS;
					for (i = 0; i < NvnS; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
			} else {
				D = OPS[0]->D_Weak;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_d(NvnS,1,NvnI,D[dim1],&Fr_vI[(eq*d+dim1)*NvnI],DFr);

					IndRHS = eq*NvnS;
					for (i = 0; i < NvnS; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
			}
			free(DFr);
			free(Fr_vI);
//printf("%d %d %d\n",NvnS,NvnI,SF_BE[P][0][0]);
//array_print_d(NvnS,Nvar,RHS,'C');
//exit(1);
		}
	} else if (strstr(Form,"Strong")) {
		printf("Exiting: Implement the strong form in compute_VOLUME_RHS_EFE.\n"), exit(1);
	}
//exit(1);

	for (i = 0; i < 2; i++)
		free(OPS[i]);
}

static void compute_VOLUMEVec_RHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form = DB.Form;
	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
	             Nvar       = DB.Nvar,
	             Neq        = DB.Neq,
	             NP         = DB.NP,
	             NECgrp     = DB.NECgrp,
	             *NVgrp     = DB.NVgrp,
	             ***SF_BE   = DB.SF_BE;

	// Standard datatypes
	unsigned int i, j, k, eq, dim1, dim2, P, curved, Eclass, iMax, jMax,
	             IndVgrp, IndC, IndF, IndFr, IndRHS, NvnI, NvnS, NV,
	             NvnS_SF[2], NvnI_SF[2], NIn[3], NOut[3], Diag[3];
	double       *What_vS, *W_vI, **What_vS_ptr, *What, *WhatVec, *F_vI, *Fr_vI,
	             **C_vI_ptr, *C_vIl, *C_vI, *C_vIVec,
	             **RHS_ptr, *RHSl, *RHS, *RHSVec, *DFr, **D, *I, *OP[3], *OP0, *OP1;

	struct S_OpCSR **D_sp;

	struct S_OPERATORS *OPS[2];
	struct S_VOLUME    *VOLUME;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	What_vS_ptr = malloc(Nvar * sizeof *What_vS_ptr); // free
	C_vI_ptr    = malloc(d*d  * sizeof *C_vI_ptr);    // free
	RHS_ptr     = malloc(Neq  * sizeof *RHS_ptr);     // free

	if (strstr(Form,"Weak")) {
		for (Eclass = 0; Eclass < NECgrp; Eclass++) {
		for (P = 0; P < NP; P++) {
		for (curved = 0; curved < 2; curved++) {
			IndVgrp = Eclass*NP*2 + P*2 + curved;
			VOLUME = DB.Vgrp[IndVgrp];

			if (VOLUME == NULL)
				continue;

			NV = NVgrp[IndVgrp];

			// Obtain operators
			init_ops(OPS[0],VOLUME,0);
			if (VOLUME->type == WEDGE)
				init_ops(OPS[1],VOLUME,1);

			// Obtain Vectorized W_vI
			NvnI = OPS[0]->NvnI;
			NvnS = OPS[0]->NvnS;
			W_vI = malloc(NvnI*NV*Nvar * sizeof *W_vI); // free
			C_vI = malloc(NvnI*NV*d*d  * sizeof *C_vI); // free

			if (Collocated) {
				for (i = 0; i < Nvar; i++)
					What_vS_ptr[i] = &W_vI[i*NV*NvnI];
			} else {
				What_vS = malloc(NvnS*NV*Nvar * sizeof *What_vS); // free

				for (i = 0; i < Nvar; i++)
					What_vS_ptr[i] = &What_vS[i*NV*NvnS];
			}
			for (i = 0, iMax = d*d; i < iMax; i++)
				C_vI_ptr[i] = &C_vI[i*NV*NvnI];

			for (i = 0, iMax = NV; i < iMax; i++) {
				if (i)
					VOLUME = VOLUME->grpnext;

				What = VOLUME->What;
				for (j = 0; j < Nvar; j++) {
					WhatVec = What_vS_ptr[j];
					for (k = 0; k < NvnS; k++)
						*WhatVec++ = *What++;
					What_vS_ptr[j] = WhatVec;
				}

				C_vIl = VOLUME->C_vI;
				for (j = 0, jMax = d*d; j < jMax; j++) {
					C_vIVec = C_vI_ptr[j];
					for (k = 0; k < NvnI; k++)
						*C_vIVec++ = *C_vIl++;
					C_vI_ptr[j] = C_vIVec;
				}
			}

			if (!Collocated) {
				if (Eclass == C_TP && SF_BE[P][0][0]) {
					for (i = 0; i < 1; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,0,0,NULL,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;

					sf_apply_d(What_vS,W_vI,NIn,NOut,NV*Nvar,OP,Diag,d);
				} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
					for (i = 0; i < 2; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,
					                  NvnS_SF[1],NvnI_SF[1],OPS[1]->ChiS_vI,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;

					sf_apply_d(What_vS,W_vI,NIn,NOut,NV*Nvar,OP,Diag,d);
				} else {
					mm_CTN_d(NvnI,NV*Nvar,NvnS,OPS[0]->ChiS_vI,What_vS,W_vI);
				}

				free(What_vS);
			}

//array_print_d(NvnI*NV,Nvar,W_vI,'C');

			// Compute Flux in reference space
			F_vI = malloc(NvnI*NV*d*Neq * sizeof *F_vI); // free
			flux_inviscid(NvnI,NV,W_vI,F_vI,d,Neq);
			free(W_vI);

//array_print_d(NvnI*NV,Nvar*d,F_vI,'C');

			Fr_vI = calloc(NvnI*NV*d*Neq , sizeof *Fr_vI); // free
			for (eq = 0; eq < Neq; eq++) {
			for (dim1 = 0; dim1 < d; dim1++) {
			for (dim2 = 0; dim2 < d; dim2++) {
				IndFr = (eq*d+dim1)*NvnI*NV;
				IndF  = (eq*d+dim2)*NvnI*NV;
				IndC  = (dim1*d+dim2)*NvnI*NV;
				for (i = 0, iMax = NvnI*NV; i < iMax; i++)
					Fr_vI[IndFr+i] += F_vI[IndF+i]*C_vI[IndC+i];
			}}}
			free(F_vI);
			free(C_vI);

//array_print_d(NvnI*NV,Nvar*d,Fr_vI,'C');

			// Compute RHS terms
			DFr = malloc(NvnS*NV * sizeof *DFr); // free
			RHS = calloc(NvnS*NV*Neq , sizeof *RHS); // free
			if (Eclass == C_TP && SF_BE[P][0][0]) {
				for (i = 0; i < 1; i++) {
					NvnS_SF[i] = OPS[i]->NvnS_SF;
					NvnI_SF[i] = OPS[i]->NvnI_SF;
				}

				I = OPS[0]->I_Weak;
				D = OPS[0]->D_Weak;

				for (dim1 = 0; dim1 < d; dim1++) {
					get_sf_parameters(NvnI_SF[0],NvnS_SF[0],I,NvnI_SF[0],NvnS_SF[0],D[0],NIn,NOut,OP,d,dim1,Eclass);

					if (Collocated) {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 2;
						Diag[dim1] = 0;
					} else {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 0;
					}

					for (eq = 0; eq < Neq; eq++) {
						sf_apply_d(&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr,NIn,NOut,NV,OP,Diag,d);

						IndRHS = eq*NvnS*NV;
						for (i = 0, iMax = NvnS*NV; i < iMax; i++)
							RHS[IndRHS+i] += DFr[i];
					}
				}
			} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
				for (i = 0; i < 2; i++) {
					NvnS_SF[i] = OPS[i]->NvnS_SF;
					NvnI_SF[i] = OPS[i]->NvnI_SF;
				}

				for (dim1 = 0; dim1 < d; dim1++) {
					if (dim1 < 2) OP0 = OPS[0]->D_Weak[dim1], OP1 = OPS[1]->I_Weak;
					else          OP0 = OPS[0]->I_Weak,       OP1 = OPS[1]->D_Weak[0];
					get_sf_parameters(NvnI_SF[0],NvnS_SF[0],OP0,NvnI_SF[1],NvnS_SF[1],OP1,NIn,NOut,OP,d,3,Eclass);

					if (Collocated) {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 2;
						if (dim1 < 2)
							Diag[0] = 0;
						else
							Diag[dim1] = 0;
					} else {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 0;
					}

					for (eq = 0; eq < Neq; eq++) {
						sf_apply_d(&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr,NIn,NOut,NV,OP,Diag,d);

						IndRHS = eq*NvnS*NV;
						for (i = 0, iMax = NvnS*NV; i < iMax; i++)
							RHS[IndRHS+i] += DFr[i];
					}
				}
			} else if (Collocated && (Eclass == C_TP || Eclass == C_WEDGE)) {
				D_sp = OPS[0]->D_Weak_sp;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_CSR_d(NvnS,NV,NvnI,D_sp[dim1],&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr);

					IndRHS = eq*NvnS*NV;
					for (i = 0, iMax = NvnS*NV; i < iMax; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
			} else {
				DFr = malloc(NvnS*NV * sizeof *DFr); // free
				D = OPS[0]->D_Weak;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_d(NvnS,NV,NvnI,D[dim1],&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr);

					IndRHS = eq*NvnS*NV;
					for (i = 0, iMax = NvnS*NV; i < iMax; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
			}
			free(DFr);
			free(Fr_vI);
//array_print_d(NvnS*NV,Neq,RHS,'C');

			// Copy values of Global RHS into local VOLUME RHS arrays
			// Profile and, if expensive, keep global RHS array
			for (i = 0; i < Neq; i++)
				RHS_ptr[i] = &RHS[i*NV*NvnS];

			for (VOLUME = DB.Vgrp[IndVgrp]; VOLUME; VOLUME = VOLUME->grpnext) {
				VOLUME->RHS = malloc(NvnS*Neq * sizeof *RHSl); // keep
				RHSl = VOLUME->RHS;

				for (eq = 0; eq < Neq; eq++) {
					RHSVec = RHS_ptr[eq];
					for (i = 0; i < NvnS; i++)
						*RHSl++ = *RHSVec++;
					RHS_ptr[eq] = RHSVec;
				}
//array_print_d(NvnS,Neq,VOLUME->RHS,'C');
			}
			free(RHS);
//exit(1);

		}}}
	}

	free(What_vS_ptr);
	free(C_vI_ptr);
	free(RHS_ptr);

	for (i = 0; i < 2; i++)
		free(OPS[i]);
}
