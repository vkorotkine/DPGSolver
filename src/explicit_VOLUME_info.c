// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS term.
 *
 *	Comments:
 *		Need to add in support for sum factorization when working (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnI, NvnS;
	double       *ChiS_vI, **D_Weak;
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
	// Standard datatypes
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	// silence
	ELEMENT_OPS = NULL;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);
	if (type == LINE || type == QUAD || type == HEX || type == WEDGE)
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else if (type == TRI || type == TET || type == PYR)
		ELEMENT_OPS = ELEMENT;
	
	OPS->NvnS = ELEMENT_OPS->NvnS[P];
	if (!curved) {
		OPS->NvnI = ELEMENT_OPS->NvnIs[P];

		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIs[P];
		OPS->D_Weak  = ELEMENT_OPS->Ds_Weak[P];
	} else {
		OPS->NvnI = ELEMENT_OPS->NvnIc[P];

		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIc[P];
		OPS->D_Weak  = ELEMENT_OPS->Dc_Weak[P];
	}
}

static void compute_VOLUME_RHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form = DB.Form;
	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
				 Nvar       = DB.Nvar,
				 Neq        = DB.Neq;

	// Standard datatypes
	unsigned int i, eq, dim1, dim2,
	             IndFr, IndF, IndC, IndRHS,
	             NnI, NvnS;
	double       *W_vI, *F_vI, *Fr_vI, *C_vI, *RHS, *DFr, **D;

	struct S_OPERATORS *OPS[2];
	struct S_VOLUME    *VOLUME;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(Form,"Weak") != NULL) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			// Obtain operators
			init_ops(OPS[0],VOLUME,0);
			if (VOLUME->type == WEDGE)
				init_ops(OPS[1],VOLUME,1);

			// Obtain W_vI
			NnI = OPS[0]->NvnI;
			if (Collocated) {
				W_vI = VOLUME->What;
			} else {
				W_vI = malloc(NnI*Nvar * sizeof *W_vI); // free
				mm_CTN_d(NnI,Nvar,OPS[0]->NvnS,OPS[0]->ChiS_vI,VOLUME->What,W_vI);
			}

			// Compute Flux in reference space
			F_vI = malloc(NnI*d*Neq * sizeof *F_vI); // free
			flux_inviscid(NnI,1,W_vI,F_vI,d,Neq);

			if (!Collocated)
				free(W_vI);

			C_vI = VOLUME->C_vI;

			Fr_vI = calloc(NnI*d*Neq , sizeof *Fr_vI); // free
			for (eq = 0; eq < Neq; eq++) {
			for (dim1 = 0; dim1 < d; dim1++) {
			for (dim2 = 0; dim2 < d; dim2++) {
				IndFr = (eq*d+dim1)*NnI;
				IndF  = (eq*d+dim2)*NnI;
				IndC  = (dim1*d+dim2)*NnI;
				for (i = 0; i < NnI; i++)
					Fr_vI[IndFr+i] += F_vI[IndF+i]*C_vI[IndC+i];
			}}}
			free(F_vI);

			// Compute RHS terms
			NvnS = OPS[0]->NvnS;

// VOLUME->RHS should be freed as soon as it is no longer needed (outside of this function)
			RHS = VOLUME->RHS;
			RHS = malloc(NvnS*Nvar * sizeof *RHS); // keep (requires external free)
			if (0 && VOLUME->Eclass == C_TP) {
				; // update this with sum factorization
			} else if (1 || VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR) {
				DFr = malloc(NvnS * sizeof *DFr); // free
				D = OPS[0]->D_Weak;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_d(NvnS,1,NnI,D[dim1],&Fr_vI[(eq*d+dim1)*NnI],DFr);
					IndRHS = eq*NnI;
					for (i = 0; i < NnI; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
				free(DFr);
			} else if (VOLUME->Eclass == C_WEDGE) {
				; // update this with sum factorization
			}
			free(Fr_vI);
		}
	} else if (strstr(Form,"Strong") != NULL) {
		printf("Exiting: Implement the strong form in compute_VOLUME_RHS_EFE.\n"), exit(1);
	}

	for (i = 0; i < 2; i++)
		free(OPS[i]);
}

static void compute_VOLUMEVec_RHS_EFE(void)
{

}
