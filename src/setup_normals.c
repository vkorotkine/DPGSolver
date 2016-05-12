#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
//#include <string.h>

#include "database.h"
//#include "parameters.h"
//#include "functions.h"

//#include "petscsys.h"

/*
 *	Purpose:
 *		Set up normals at integration nodes on element facets.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
 */

/*
static void init_ops(const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     unsigned int *NvnC, unsigned int *NvnI, unsigned int *NfnI,
                     double *I_vC_vI, double **I_vC_fI, double *nr);

void setup_normals(struct S_FACET *FACET)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int fn, fnMax, curved,
	             IndfIn,
	             NvnC0, NvnI0, NfnI0, NvnC1, NvnI1, NfnI1, NnI;
	double       *I_vC_vI0, **I_vC_fI0, *nr0,
	             *I_vC_vI1, **I_vC_fI1, *nr1,
	             *C_fI, *C_vC, *nrIn, *n;

	struct S_VOLUME *VOLUMEIn, *VOLUMEOut;

	VOLUMEIn  = FACET->VOLUMEIn;
	VOLUMEOut = FACET->VOLUMEOut;

//  IndfIn = 
	curved = FACET->curved;

	init_ops(VOLUMEIn,FACET,&NvnC0,&NvnI0,&NfnI0,&I_vC_vI0,&I_vC_fI0,&nr0);
	if (VOLUMEIn->type == WEDGE)
		init_ops(VOLUMEIn,FACET,&NvnC1,&NvnI1,&NfnI1,&I_vC_vI1,&I_vC_fI1,&nr1);

	C_vC = VOLUMEIn->C_vC;
	if (VOLUMEIn->Eclass == C_TP) {
		printf("Add in support for C_TP in setup_normals.\n");
		exit(1);
	} else if (VOLUMEIn->Eclass == C_SI || VOLUMEIn->Eclass == C_PYR) {
		C_fI = malloc(NvnI0*d*d * sizeof *C_fI); // free

		mm_CTN_d(NfnI0,d*d,NvnC0,I_vC_fI0[IndfIn],C_vC,C_fI);

		NnI = NfnI0;
	} else if (VOLUMEIn->Eclass == C_WEDGE) {
		printf("Add in support for C_WEDGE in setup_normals.\n");
		exit(1);
	}
	nrIn = &nr0[IndfIn*3];

	// Store a single normal on straight FACETs
	if (!curved) fnMax = 1;
	else         fnMax = NnI;

	n = calloc(fnMax*d , sizeof *n); // keep
	for (fn = 0; fn < fnMax; fn++) {
		for (dim1 = 0, dim1 < d; dim1++) {
		for (dim2 = 0, dim2 < d; dim2++) {
			n[fn*d+dim1] += nrIn[dim2]*C_fI[NnI*(dim1+d*dim2)+fn];
		}}
	}

array_print_d(fnMax,d,n,'R');
exit(1);



	free(C_vI);
}

static void init_ops(const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     unsigned int *NvnC, unsigned int *NvnI, unsigned int *NfnI,
                     double *I_vC_vI, double **I_vC_fI, double *nr)
{
	// Standard datatypes
	unsigned int P, Vtype, Vcurved, FtypeInt;

	P        = FACET->P;
	Vtype    = VOLUME->type;
	Vcurved  = VOLUME->curved;
	FtypeInt = FACET->typeInt;

	ELEMENT = get_ELEMENT_type(Vtype);
	if (type == LINE || type == QUAD || type == HEX || type == WEDGE)
		ELEMENT_Ops = ELEMENT->ELEMENTclass[IndClass];
	else if (type == TRI || type == TET || type == PYR)
		ELEMENT_Ops = ELEMENT;

	if (!Vcurved) {
		// Straight VOLUME
		*NvnC = ELEMENT_Ops->NvnCs[P];
		if (FtypeInt == 's') {
			// Straight FACET Integration
			*NvnI = ELEMENT_Ops->NvnIs[P];
			*NfnI = ELEMENT_Ops->NfnIs[P];

			*I_vC_vI = ELEMENT_Ops->I_vCs_vIs[P];
			*I_vC_fI = ELEMENT_Ops->I_vCs_vfs[P];
		} else {
			// Curved FACET Integration
			*NvnI = ELEMENT_Ops->NvnIc[P];
			*NfnI = ELEMENT_Ops->NfnIc[P];

			*I_vC_vI = ELEMENT_Ops->I_vCs_vIc[P];
			*I_vC_fI = ELEMENT_Ops->I_vCs_vfc[P];
		}
	} else {
		// Curved VOLUME
		*NvnC = ELEMENT_Ops->NvnCc[P];
		if (FtypeInt == 's') {
			// Straight FACET Integration
			*NvnI = ELEMENT_Ops->NvnIs[P];
			*NfnI = ELEMENT_Ops->NfnIs[P];

			*I_vC_vI = ELEMENT_Ops->I_vCc_vIs[P];
			*I_vC_fI = ELEMENT_Ops->I_vCc_vfs[P];
		} else {
			// Curved FACET Integration
			*NvnI = ELEMENT_Ops->NvnIc[P];
			*NfnI = ELEMENT_Ops->NfnIc[P];

			*I_vC_vI = ELEMENT_Ops->I_vCc_vIc[P];
			*I_vC_fI = ELEMENT_Ops->I_vCc_vfc[P];
		}
	}
	*nr = ELEMENT->nr;
}
*/
