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

struct S_OPERATORS {
	unsigned int NvnG, NvnC, NvnI;
	double       *IC, *I_vG_vC, *I_vG_vI, *I_vC_vI,
	             **D_vG_vC, **D_vG_vI, **D_vC_vC;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass);

void setup_geom_factors(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int i, n, row, col,
	             NvnG0, NvnC0, NvnI0, NvnG1, NvnC1, NvnI1;
	double       *XYZ, *J_vI, *J_vC, *detJV_vI, *C_vC, *C_vI;

	struct S_OPERATORS *OPS[2];

	// silence
	detJV_vI = NULL;
	C_vC     = NULL;
	C_vI     = NULL;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	// Obtain operators
	init_ops(OPS[0],VOLUME,0);
	if (VOLUME->type == WEDGE)
		init_ops(OPS[1],VOLUME,1);

	NvnG0 = OPS[0]->NvnG;
	NvnC0 = OPS[0]->NvnC;
	NvnI0 = OPS[0]->NvnI;

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
			mm_CTN_d(NvnI0,1,NvnG0,OPS[0]->D_vG_vI[col],&XYZ[NvnG0*row],&J_vI[NvnI0*(d*row+col)]);
			mm_CTN_d(NvnC0,1,NvnG0,OPS[0]->D_vG_vC[col],&XYZ[NvnG0*row],&J_vC[NvnC0*(d*row+col)]);
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
				detJV_vI[n] =   J_vI[NvnI0*(d*0+0)+n]*J_vI[NvnI0*(d*1+1)+n]
				              - J_vI[NvnI0*(d*0+1)+n]*J_vI[NvnI0*(d*1+0)+n];
			}

			for (n = 0; n < NvnC0; n++) {
				C_vC[NvnC0*(d*0+0)+n] =  J_vC[NvnC0*(d*1+1)+n];
				C_vC[NvnC0*(d*0+1)+n] = -J_vC[NvnC0*(d*1+0)+n];
				C_vC[NvnC0*(d*1+0)+n] = -J_vC[NvnC0*(d*0+1)+n];
				C_vC[NvnC0*(d*1+1)+n] =  J_vC[NvnC0*(d*0+0)+n];
			}

		} else if (d == 3) {
			for (n = 0; n < NvnI0; n++) {
				detJV_vI[n] =   J_vI[NvnI0*(d*0+0)+n]*(  J_vI[NvnI0*(d*1+1)+n]*J_vI[NvnI0*(d*2+2)+n]
				                                       - J_vI[NvnI0*(d*1+2)+n]*J_vI[NvnI0*(d*2+1)+n])
				              - J_vI[NvnI0*(d*0+1)+n]*(  J_vI[NvnI0*(d*1+0)+n]*J_vI[NvnI0*(d*2+2)+n]
				                                       - J_vI[NvnI0*(d*1+2)+n]*J_vI[NvnI0*(d*2+0)+n])
				              + J_vI[NvnI0*(d*0+2)+n]*(  J_vI[NvnI0*(d*1+0)+n]*J_vI[NvnI0*(d*2+1)+n]
				                                       - J_vI[NvnI0*(d*1+1)+n]*J_vI[NvnI0*(d*2+0)+n]);
			}

/*
			// standard form
			for (n = 0; n < NvnC0; n++) {
				C_vC[NvnC0*(d*0+0)+n] = J_vC[NvnC0*(d*1+1)+n]*J_vC[NvnC0*(d*2+2)+n]
				                       -J_vC[NvnC0*(d*1+2)+n]*J_vC[NvnC0*(d*2+1)+n];
				C_vC[NvnC0*(d*0+1)+n] = J_vC[NvnC0*(d*1+2)+n]*J_vC[NvnC0*(d*2+0)+n]
				                       -J_vC[NvnC0*(d*1+0)+n]*J_vC[NvnC0*(d*2+2)+n];
				C_vC[NvnC0*(d*0+2)+n] = J_vC[NvnC0*(d*1+0)+n]*J_vC[NvnC0*(d*2+1)+n]
				                       -J_vC[NvnC0*(d*1+1)+n]*J_vC[NvnC0*(d*2+0)+n];
				C_vC[NvnC0*(d*1+0)+n] = J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*2+1)+n]
				                       -J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*2+2)+n];
				C_vC[NvnC0*(d*1+1)+n] = J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*2+2)+n]
				                       -J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*2+0)+n];
				C_vC[NvnC0*(d*1+2)+n] = J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*2+0)+n]
				                       -J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*2+1)+n];
				C_vC[NvnC0*(d*2+0)+n] = J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*1+2)+n]
				                       -J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*1+1)+n];
				C_vC[NvnC0*(d*2+1)+n] = J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*1+0)+n]
				                       -J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*1+2)+n];
				C_vC[NvnC0*(d*2+2)+n] = J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*1+1)+n]
				                       -J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*1+0)+n];
			}
*/

			// curl form
			unsigned int dim, IndCurl, IndXYZJ, IndCurlXYZJ;
			unsigned int OrderCurl[12] = {1, 2, 2, 1, 2, 0, 0, 2, 0, 1, 1, 0};
			double *XYZJ_vC, *XYZ_vC, *tmp_vC, *CurlXYZJ_vC;

			XYZ_vC = malloc(d*NvnC0 * sizeof *XYZ_vC); // free
			XYZJ_vC = malloc(6*3*NvnC0 * sizeof *XYZJ_vC); // free

			mm_CTN_d(NvnC0,d,NvnG0,OPS[0]->I_vG_vC,XYZ,XYZ_vC);


			IndXYZJ = 0;
			IndCurl = 0;
			for (i = 0; i < 6; i++) {
				for (dim = 0; dim < d; dim++) {
					for (n = 0; n < NvnC0; n++) {
						XYZJ_vC[NvnC0*IndXYZJ+n] = XYZ_vC[NvnC0*(OrderCurl[IndCurl])+n]
						                         * J_vC[NvnC0*(OrderCurl[IndCurl+1]*d+dim)+n];
					}
					IndXYZJ++;
				}
				IndCurl += 2;
			}
			free(XYZ_vC);

			tmp_vC = malloc(NvnC0 * sizeof *tmp_vC); // free
			CurlXYZJ_vC = malloc(6*3*NvnC0 * sizeof *CurlXYZJ_vC); // free

			IndXYZJ = 0;
			IndCurlXYZJ = 0;
			for (i = 0; i < 6; i++) {
				IndCurl = 0;
				for (dim = 0; dim < d; dim++) {
					mm_CTN_d(NvnC0,1,NvnC0,OPS[0]->D_vC_vC[OrderCurl[IndCurl+0]],
					         &XYZJ_vC[NvnC0*(i*3+OrderCurl[IndCurl+2])],&CurlXYZJ_vC[NvnC0*IndCurlXYZJ]);
					mm_CTN_d(NvnC0,1,NvnC0,OPS[0]->D_vC_vC[OrderCurl[IndCurl+1]],
					         &XYZJ_vC[NvnC0*(i*3+OrderCurl[IndCurl+2+1])],tmp_vC);
					for (n = 0; n < NvnC0; n++)
						CurlXYZJ_vC[NvnC0*IndCurlXYZJ+n] -= tmp_vC[n];

					IndCurlXYZJ++;
					IndCurl += 4;
				}
				IndXYZJ += 3;
			}
			free(XYZJ_vC);
			free(tmp_vC);

//array_print_d(NvnC0,9,C_vC,'C');

			IndCurlXYZJ = 0;
			for (row = 0; row < d; row++) {
				for (col = 0; col < d; col++) {
					for (n = 0; n < NvnC0; n++) {
						C_vC[NvnC0*(row*d+col)+n] =
							0.5*(CurlXYZJ_vC[NvnC0*IndCurlXYZJ+n]-CurlXYZJ_vC[NvnC0*(IndCurlXYZJ+3)+n]);
					}

					IndCurlXYZJ++;
				}
				IndCurlXYZJ += 3;
			}

//array_print_d(NvnC0,9,C_vC,'C');

			free(CurlXYZJ_vC);
		}
		mm_CTN_d(NvnI0,d*d,NvnC0,OPS[0]->I_vC_vI,C_vC,C_vI);

/*
array_print_d(NvnG0,d,XYZ,'C');
for (int dim = 0; dim < d; dim++) {
	array_print_d(NvnC0,d,&J_vC[NvnC0*(d*dim)],'C');
}
*/

		free(J_vI);
		free(J_vC);

//printf("Exiting setup_geom_factors.\n"), exit(1);
	} else if (VOLUME->Eclass == C_WEDGE) {
		NvnG1 = OPS[1]->NvnG;
		NvnC1 = OPS[1]->NvnC;
		NvnI1 = OPS[1]->NvnI;
		printf("Add in support for C_TP in setup_geom_factors.\n");
		exit(1);
	}
	VOLUME->detJV_vI = detJV_vI;
	VOLUME->C_vC     = C_vC;
	VOLUME->C_vI     = C_vI;

	for (i = 0; i < 2; i++)
		free(OPS[i]);
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

	if (!curved) {
		OPS->NvnG = ELEMENT_OPS->NvnGs[0];
		OPS->NvnC = ELEMENT_OPS->NvnCs[P];
		OPS->NvnI = ELEMENT_OPS->NvnIs[P];

		OPS->IC      = ELEMENT_OPS->ICs[P];
		OPS->I_vG_vC = ELEMENT_OPS->I_vGs_vCs[P];
		OPS->I_vG_vI = ELEMENT_OPS->I_vGs_vIs[P];
		OPS->I_vC_vI = ELEMENT_OPS->I_vCs_vIs[P];
		OPS->D_vG_vC = ELEMENT_OPS->D_vGs_vCs[P];
		OPS->D_vG_vI = ELEMENT_OPS->D_vGs_vIs[P];
		OPS->D_vC_vC = ELEMENT_OPS->D_vCs_vCs[P];
	} else {
		OPS->NvnG = ELEMENT_OPS->NvnGc[P];
		OPS->NvnC = ELEMENT_OPS->NvnCc[P];
		OPS->NvnI = ELEMENT_OPS->NvnIc[P];

		OPS->IC      = ELEMENT_OPS->ICc[P];
		OPS->I_vG_vC = ELEMENT_OPS->I_vGc_vCc[P];
		OPS->I_vG_vI = ELEMENT_OPS->I_vGc_vIc[P];
		OPS->I_vC_vI = ELEMENT_OPS->I_vCc_vIc[P];
		OPS->D_vG_vC = ELEMENT_OPS->D_vGc_vCc[P];
		OPS->D_vG_vI = ELEMENT_OPS->D_vGc_vIc[P];
		OPS->D_vC_vC = ELEMENT_OPS->D_vCc_vCc[P];
	}
}
