#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up geometric factors.
 *
 *	Comments:
 *		Only the curl-form of the cofactor matrix terms is used following the analysis of Kopriva(2006)
 *		The second set of operators is only initialized and used for WEDGE elements.
 *		this may be slightly less readable. Profile and decide. (ToBeModified)
 *
 *		Investigation into the importance of proper treatment of the metric terms is still required. For this reason,
 *		the cofactor matrix terms are first projected to the basis of order PC and then projected to the VOLUME and
 *		FACET integration nodes. After the code is working, modify this by trying different combinations, such as
 *		directly computing C_vI, C_fI and make conclusions. (ToBeDeleted)
 *		detJ is computed directly at the VOLUME integration nodes. Thus, likely don't need PJs/PJc in setup_parameters
 *		=> Delete later (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 *		Kopriva(2006)-Metric_Identities_and_the_Discontinuous_Spectral_Element_Method_on_Curvilinear_Meshes
 */

static void init_ops(const struct S_VOLUME *VOLUME, const unsigned int IndClass,
                     unsigned int *NvnG, unsigned int *NvnC, unsigned int *NvnI,
                     double **IC, double **I_vG_vC, double **I_vG_vI, double **I_vC_vI,
                     double ***D_vG_vC, double ***D_vG_vI, double ***D_vC_vC);

void setup_geom_factors(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, row, col,
	             NvnG0, NvnC0, NvnI0, NvnG1, NvnC1, NvnI1;
	double       *IC0, *I_vG_vC0, *I_vG_vI0, *I_vC_vI0, **D_vG_vC0, **D_vG_vI0, **D_vC_vC0,
	             *IC1, *I_vG_vC1, *I_vG_vI1, *I_vC_vI1, **D_vG_vC1, **D_vG_vI1, **D_vC_vC1,
	             *XYZ, *J_vI, *J_vC, *detJV_vI, *C_vC, *C_vI;

	// silence
	detJV_vI = NULL;
	C_vC     = NULL;
	C_vI     = NULL;

	// Obtain operators
	init_ops(VOLUME,0,&NvnG0,&NvnC0,&NvnI0,&IC0,&I_vG_vC0,&I_vG_vI0,&I_vC_vI0,&D_vG_vC0,&D_vG_vI0,&D_vC_vC0);
	if (VOLUME->type == WEDGE)
		init_ops(VOLUME,1,&NvnG1,&NvnC1,&NvnI1,&IC1,&I_vG_vC1,&I_vG_vI1,&I_vC_vI1,&D_vG_vC1,&D_vG_vI1,&D_vC_vC1);

	XYZ = VOLUME->XYZ;
	if (VOLUME->Eclass == C_TP) {
		printf("Add in support for C_TP in setup_geom_factors.\n");
		exit(1);
	} else if (VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR) {
		J_vI     = malloc(NvnI0*d*d * sizeof *J_vI);     // free
		J_vC     = malloc(NvnC0*d*d * sizeof *J_vC);     // free
		detJV_vI = malloc(NvnI0     * sizeof *detJV_vI); // keep
		C_vC     = malloc(NvnC0*d*d * sizeof *C_vC);     // keep (free after setup_normals)
		C_vI     = malloc(NvnI0*d*d * sizeof *C_vI);     // keep

		for (row = 0; row < d; row++) {
		for (col = 0; col < d; col++) {
			mm_CTN_d(NvnI0,1,NvnG0,D_vG_vI0[col],&XYZ[NvnG0*row],&J_vI[NvnI0*(row+d*col)]);
			mm_CTN_d(NvnC0,1,NvnG0,D_vG_vC0[col],&XYZ[NvnG0*row],&J_vC[NvnC0*(row+d*col)]);
		}}

		if (d == 1) {
			for (n = 0; n < NvnI0; n++) {
				detJV_vI[n] = J_vI[n];
			}

			for (n = 0; n < NvnC0; n++) {
				C_vC[n] = 1.0;
			}
		} else if (d == 2) {
			for (n = 0; n < NvnI0; n++) {
				detJV_vI[n] =   J_vI[NvnI0*(0+d*0)+n]*J_vI[NvnI0*(1+d*1)+n]
				              - J_vI[NvnI0*(0+d*1)+n]*J_vI[NvnI0*(1+d*0)+n];
			}

			for (n = 0; n < NvnC0; n++) {
				C_vC[NvnC0*(0+d*0)+n] = pow(-1.0,1.0+1.0)*J_vC[NvnC0*(1+d*1)+n];
				C_vC[NvnC0*(0+d*1)+n] = pow(-1.0,1.0+2.0)*J_vC[NvnC0*(1+d*0)+n];
				C_vC[NvnC0*(1+d*0)+n] = pow(-1.0,2.0+1.0)*J_vC[NvnC0*(0+d*1)+n];
				C_vC[NvnC0*(1+d*1)+n] = pow(-1.0,2.0+2.0)*J_vC[NvnC0*(0+d*0)+n];
			}

		} else if (d == 3) {
			for (n = 0; n < NvnI0; n++) {
				detJV_vI[n] =   J_vI[NvnI0*(0+d*0)+n]*(  J_vI[NvnI0*(1+d*1)+n]*J_vI[NvnI0*(2+d*2)+n]
				                                       - J_vI[NvnI0*(1+d*2)+n]*J_vI[NvnI0*(2+d*1)+n])
				              - J_vI[NvnI0*(0+d*1)+n]*(  J_vI[NvnI0*(1+d*0)+n]*J_vI[NvnI0*(2+d*2)+n]
				                                       - J_vI[NvnI0*(1+d*2)+n]*J_vI[NvnI0*(2+d*0)+n])
				              + J_vI[NvnI0*(0+d*2)+n]*(  J_vI[NvnI0*(1+d*0)+n]*J_vI[NvnI0*(2+d*1)+n]
				                                       - J_vI[NvnI0*(1+d*1)+n]*J_vI[NvnI0*(2+d*0)+n]);
			}

			// standard form
			for (n = 0; n < NvnC0; n++) {
				C_vC[NvnC0*(0+d*0)+n] = pow(-1.0,1.0+1.0)*(  J_vC[NvnC0*(1+d*1)+n]*J_vC[NvnC0*(2+d*2)+n]
				                                           - J_vC[NvnC0*(1+d*2)+n]*J_vC[NvnC0*(2+d*1)+n]);
				C_vC[NvnC0*(0+d*1)+n] = pow(-1.0,1.0+2.0)*(  J_vC[NvnC0*(1+d*2)+n]*J_vC[NvnC0*(2+d*0)+n]
				                                           - J_vC[NvnC0*(1+d*0)+n]*J_vC[NvnC0*(2+d*2)+n]);
				C_vC[NvnC0*(0+d*2)+n] = pow(-1.0,1.0+3.0)*(  J_vC[NvnC0*(1+d*0)+n]*J_vC[NvnC0*(2+d*1)+n]
				                                           - J_vC[NvnC0*(1+d*1)+n]*J_vC[NvnC0*(2+d*0)+n]);
				C_vC[NvnC0*(1+d*0)+n] = pow(-1.0,2.0+1.0)*(  J_vC[NvnC0*(0+d*2)+n]*J_vC[NvnC0*(2+d*1)+n]
				                                           - J_vC[NvnC0*(0+d*1)+n]*J_vC[NvnC0*(2+d*2)+n]);
				C_vC[NvnC0*(1+d*1)+n] = pow(-1.0,2.0+2.0)*(  J_vC[NvnC0*(0+d*0)+n]*J_vC[NvnC0*(2+d*2)+n]
				                                           - J_vC[NvnC0*(0+d*2)+n]*J_vC[NvnC0*(2+d*0)+n]);
				C_vC[NvnC0*(1+d*2)+n] = pow(-1.0,2.0+3.0)*(  J_vC[NvnC0*(0+d*1)+n]*J_vC[NvnC0*(2+d*0)+n]
				                                           - J_vC[NvnC0*(0+d*0)+n]*J_vC[NvnC0*(2+d*1)+n]);
				C_vC[NvnC0*(2+d*0)+n] = pow(-1.0,3.0+1.0)*(  J_vC[NvnC0*(1+d*1)+n]*J_vC[NvnC0*(1+d*2)+n]
				                                           - J_vC[NvnC0*(1+d*2)+n]*J_vC[NvnC0*(1+d*1)+n]);
				C_vC[NvnC0*(2+d*1)+n] = pow(-1.0,3.0+2.0)*(  J_vC[NvnC0*(1+d*2)+n]*J_vC[NvnC0*(1+d*0)+n]
				                                           - J_vC[NvnC0*(1+d*0)+n]*J_vC[NvnC0*(1+d*2)+n]);
				C_vC[NvnC0*(2+d*2)+n] = pow(-1.0,3.0+3.0)*(  J_vC[NvnC0*(1+d*0)+n]*J_vC[NvnC0*(1+d*1)+n]
				                                           - J_vC[NvnC0*(1+d*1)+n]*J_vC[NvnC0*(1+d*0)+n]);
			}

			// check correctness of cofactor terms by looking at normal vectors.
			// implement curl-form and check that values are approximately equal
		}
		mm_CTN_d(NvnI0,d*d,NvnC0,I_vC_vI0,C_vC,C_vI);

/*
array_print_d(NvnG0,d,XYZ,'C');
for (dim = 0; dim < d; dim++) {
	array_print_d(NvnC0,d,&J_vC[NvnC0*(d*dim)],'C');
}
*/


		free(J_vI);
		free(J_vC);

//printf("Exiting setup_geom_factors.\n"), exit(1);
	} else if (VOLUME->Eclass == C_WEDGE) {
		printf("Add in support for C_TP in setup_geom_factors.\n");
		exit(1);
	}
	VOLUME->detJV_vI = detJV_vI;
	VOLUME->C_vC     = C_vC;
	VOLUME->C_vI     = C_vI;
}

static void init_ops(const struct S_VOLUME *VOLUME, const unsigned int IndClass,
                     unsigned int *NvnG, unsigned int *NvnC, unsigned int *NvnI,
                     double **IC, double **I_vG_vC, double **I_vG_vI, double **I_vC_vI,
                     double ***D_vG_vC, double ***D_vG_vI, double ***D_vC_vC)
{
	// Standard datatypes
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT, *ELEMENT_Ops;

	// silence
	ELEMENT_Ops = NULL;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);
	if (type == LINE || type == QUAD || type == HEX || type == WEDGE)
		ELEMENT_Ops = ELEMENT->ELEMENTclass[IndClass];
	else if (type == TRI || type == TET || type == PYR)
		ELEMENT_Ops = ELEMENT;

	if (!curved) {
		*NvnG = ELEMENT_Ops->NvnGs[0];
		*NvnC = ELEMENT_Ops->NvnCs[P];
		*NvnI = ELEMENT_Ops->NvnIs[P];

		*IC      = ELEMENT_Ops->ICs[P];
		*I_vG_vC = ELEMENT_Ops->I_vGs_vCs[P];
		*I_vG_vI = ELEMENT_Ops->I_vGs_vIs[P];
		*I_vC_vI = ELEMENT_Ops->I_vCs_vIs[P];
		*D_vG_vC = ELEMENT_Ops->D_vGs_vCs[P];
		*D_vG_vI = ELEMENT_Ops->D_vGs_vIs[P];
		*D_vC_vC = ELEMENT_Ops->D_vCs_vCs[P];
	} else {
		*NvnG = ELEMENT_Ops->NvnGc[P];
		*NvnC = ELEMENT_Ops->NvnCc[P];
		*NvnI = ELEMENT_Ops->NvnIc[P];

		*IC      = ELEMENT_Ops->ICc[P];
		*I_vG_vC = ELEMENT_Ops->I_vGc_vCc[P];
		*I_vG_vI = ELEMENT_Ops->I_vGc_vIc[P];
		*I_vC_vI = ELEMENT_Ops->I_vCc_vIc[P];
		*D_vG_vC = ELEMENT_Ops->D_vGc_vCc[P];
		*D_vG_vI = ELEMENT_Ops->D_vGc_vIc[P];
		*D_vC_vC = ELEMENT_Ops->D_vCc_vCc[P];
	}
}
