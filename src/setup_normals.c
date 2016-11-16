// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_normals.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Set up normals at integration nodes on element faces.
 *
 *	Comments:
 *		The unit normal is stored so that the correct normal velocities are computed in the numerical flux functions.
 *
 *	Notation:
 *		n_fI     : Physical unit normal vector evaluated at the (f)ace (I)ntegration nodes.
 *		detJF_fI : Area element evaluated at the (f)ace (I)ntegration nodes.
 *
 *	References:
 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
 */

struct S_OPERATORS {
	unsigned int NvnC, NvnI, NfnI, NfnS;
	double       *I_vC_vI, **I_vC_fI, **I_vC_fS, *nr;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                     const unsigned int IndClass);

void setup_normals(struct S_FACE *FACE)
{
	// Initialize DB Parameters
	unsigned int d        = DB.d,
	             Adapt    = DB.Adapt,
	             NfrefMax = DB.NfrefMax;

	// Standard datatypes
	unsigned int fn, fnMax, curved, dim, dim1, dim2,
	             VfIn, Eclass, IndFType, fIn,
	             NvnC0, NfnI0, NfnS0, Nn;
	double       nSum, nSum2,
	             *C_fI, *C_fS, *C_vC, *nrIn, *n_fI, *n_fS, *detJF_fI, *detJF_fS;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VIn;

	OPS = malloc(sizeof *OPS); // free

	VIn  = FACE->VIn;

	VfIn = FACE->VfIn;
	curved = FACE->curved;

	fIn = VfIn/NfrefMax;

	Eclass = get_Eclass(VIn->type);
	IndFType = get_IndFType(Eclass,fIn);

	init_ops(OPS,VIn,FACE,IndFType);

	NvnC0 = OPS->NvnC;
	NfnI0 = OPS->NfnI;
	NfnS0 = OPS->NfnS;

	C_vC = VIn->C_vC;

	if (VfIn % NFREFMAX != 0)
		printf("Error: VfIn should be h-conforming in setup_normals.\n"), exit(1);

	nrIn = &(OPS->nr[fIn*d]);
	switch (Adapt) {
	default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in setup_normals.\n"), exit(1);
		C_fS = malloc(NfnS0*d*d * sizeof *C_fS); // free
		mm_CTN_d(NfnS0,d*d,NvnC0,OPS->I_vC_fS[VfIn],C_vC,C_fS);

		Nn = NfnS0;

		// Store a single normal on straight FACEs
		if (!curved) fnMax = 1;
		else         fnMax = Nn;

		n_fS     = calloc(fnMax*d , sizeof *n_fS);     // keep
		detJF_fS = calloc(fnMax   , sizeof *detJF_fS); // keep
		for (fn = 0; fn < fnMax; fn++) {
			for (dim1 = 0; dim1 < d; dim1++) {
			for (dim2 = 0; dim2 < d; dim2++) {
				n_fS[fn*d+dim1] += nrIn[dim2]*C_fS[Nn*(dim1+d*dim2)+fn];
			}}
		}
		free(C_fS);

		for (fn = 0; fn < fnMax; fn++) {
			nSum2 = 0;
			for (dim = 0; dim < d; dim++)
				nSum2 += pow(n_fS[fn*d+dim],2.0);

			nSum = sqrt(nSum2);
			detJF_fS[fn] = nSum;
			for (dim = 0; dim < d; dim++)
				n_fS[fn*d+dim] /= nSum;
		}

		if (FACE->n_fS)
			free(FACE->n_fS);
		FACE->n_fS = n_fS;

		if (FACE->detJF_fS)
			free(FACE->detJF_fS);
		FACE->detJF_fS = detJF_fS;
		break;
case ADAPT_P: // ToBeModified
case ADAPT_H:
case ADAPT_HP:
	case ADAPT_0:
		C_fI = malloc(NfnI0*d*d * sizeof *C_fI); // free
		mm_CTN_d(NfnI0,d*d,NvnC0,OPS->I_vC_fI[VfIn],C_vC,C_fI);

		Nn = NfnI0;

		// Potentially store a single normal on straight FACEs (ToBeModified)
//		if (!curved) fnMax = 1;
//		else         fnMax = Nn;
		fnMax = Nn;

		n_fI     = calloc(fnMax*d , sizeof *n_fI);     // keep
		detJF_fI = calloc(fnMax   , sizeof *detJF_fI); // keep
		for (fn = 0; fn < fnMax; fn++) {
			for (dim1 = 0; dim1 < d; dim1++) {
			for (dim2 = 0; dim2 < d; dim2++) {
				n_fI[fn*d+dim1] += nrIn[dim2]*C_fI[Nn*(dim1+d*dim2)+fn];
			}}
		}
		free(C_fI);

		for (fn = 0; fn < fnMax; fn++) {
			nSum2 = 0;
			for (dim = 0; dim < d; dim++)
				nSum2 += pow(n_fI[fn*d+dim],2.0);

			nSum = sqrt(nSum2);
			detJF_fI[fn] = nSum;
			for (dim = 0; dim < d; dim++)
				n_fI[fn*d+dim] /= nSum;
		}

		if (FACE->n_fI)
			free(FACE->n_fI);
		FACE->n_fI = n_fI;

		if (FACE->detJF_fI)
			free(FACE->detJF_fI);
		FACE->detJF_fI = detJF_fI;
		break;
	}

//printf("%d %d %d\n",FACE->indexg,VfIn,IndFType);
//array_print_d(fnMax,d,n,'R');

	free(OPS);
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                     const unsigned int IndClass)
{
	// Standard datatypes
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	PV       = VOLUME->P;
	PF       = FACE->P;
	Vtype    = VOLUME->type;
	Vcurved  = VOLUME->curved;
	FtypeInt = FACE->typeInt;

	ELEMENT     = get_ELEMENT_type(Vtype);
	ELEMENT_OPS = ELEMENT;

	OPS->NfnS    = ELEMENT_OPS->NfnS[PF][IndClass];
	if (!Vcurved) {
		// Straight VOLUME
		OPS->NvnC = ELEMENT_OPS->NvnCs[PV];
		OPS->NvnI = ELEMENT_OPS->NvnIs[PV];

		OPS->I_vC_fS = ELEMENT_OPS->I_vCs_fS[PV][PF];
		if (FtypeInt == 's') {
			// Straight FACE Integration
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCs_vIs[PV][PV][0];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCs_fIs[PV][PF];
		} else {
			// Curved FACE Integration
			OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCs_vIc[PV][PV][0];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCs_fIc[PV][PF];
		}
	} else {
		// Curved VOLUME
		OPS->NvnC = ELEMENT_OPS->NvnCc[PV];
		OPS->NvnI = ELEMENT_OPS->NvnIc[PV];

		OPS->I_vC_fS = ELEMENT_OPS->I_vCc_fS[PV][PF];
		if (FtypeInt == 's') {
			// Straight FACE Integration
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCc_vIs[PV][PV][0];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCc_fIs[PV][PF];
		} else {
			// Curved FACE Integration
			OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];

			OPS->I_vC_vI = ELEMENT_OPS->I_vCc_vIc[PV][PV][0];
			OPS->I_vC_fI = ELEMENT_OPS->I_vCc_fIc[PV][PF];
		}
	}
	OPS->nr = ELEMENT->nr;
}
