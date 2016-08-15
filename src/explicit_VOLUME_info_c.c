// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "explicit_VOLUME_info_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "fluxes_inviscid_c.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Identical to explicit_VOLUME_info using complex variables (for complex step verification).
 *
 *	Comments:
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

void explicit_VOLUME_info_c(void)
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
			printf("Error: Vectorized version not yet converted.\n"), EXIT_MSG;
//			compute_VOLUMEVec_RHS_EFE();
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
	P = IndClass;
	ELEMENT_OPS = NULL;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);
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
				 Neq        = DB.Neq;

	// Standard datatypes
	unsigned int   i, eq, dim1, dim2,
	               IndFr, IndF, IndC, IndRHS,
	               NvnI, NvnS;
	double         *C_vI, **D;
	double complex *W_vI, *F_vI, *Fr_vI, *RHS, *DFr;

	struct S_OPERATORS *OPS[2];
	struct S_VOLUME    *VOLUME;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(Form,"Weak")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			// Obtain operators
			init_ops(OPS[0],VOLUME,0);
			if (VOLUME->type == WEDGE)
				init_ops(OPS[1],VOLUME,1);

			// Obtain W_vI
			NvnI = OPS[0]->NvnI;
			if (Collocated) {
				W_vI = VOLUME->What_c;
			} else {
				W_vI = malloc(NvnI*Nvar * sizeof *W_vI); // free

				mm_dcc(CBCM,CBT,CBNT,NvnI,Nvar,OPS[0]->NvnS,1.0,OPS[0]->ChiS_vI,VOLUME->What_c,W_vI);
			}

			// Compute Flux in reference space
			F_vI = malloc(NvnI*d*Neq * sizeof *F_vI); // free
			flux_inviscid_c(NvnI,1,W_vI,F_vI,d,Neq);

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

			// Compute RHS terms
			NvnS = OPS[0]->NvnS;

			if (VOLUME->RHS_c)
				free(VOLUME->RHS_c);
			RHS = calloc(NvnS*Neq , sizeof *RHS); // keep (requires external free)
			VOLUME->RHS_c = RHS;

			DFr = malloc(NvnS * sizeof *DFr); // free

			D = OPS[0]->D_Weak;

			for (eq = 0; eq < Neq; eq++) {
			for (dim1 = 0; dim1 < d; dim1++) {
				mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnI,1.0,D[dim1],&Fr_vI[(eq*d+dim1)*NvnI],DFr);

				IndRHS = eq*NvnS;
				for (i = 0; i < NvnS; i++)
					RHS[IndRHS+i] += DFr[i];
			}}

			free(DFr);
			free(Fr_vI);
		}
	} else if (strstr(Form,"Strong")) {
		printf("Exiting: Implement the strong form in compute_VOLUME_RHS_EFE.\n"), exit(1);
	}

	for (i = 0; i < 2; i++)
		free(OPS[i]);
}
