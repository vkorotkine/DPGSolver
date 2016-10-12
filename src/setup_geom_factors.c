// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_geom_factors.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Set up geometric factors.
 *
 *	Comments:
 *		Only the curl-form of the cofactor matrix terms is used following the analysis of Kopriva(2006).
 *		C is stored sequentially as C = [ C11 C21 C31 C12 C22 C32 C13 C23 C33 ]. This is done so as to reduce memory
 *		stride while computing reference flux terms as [fr]_(1xd) = [f]_(1xd)*[C]_(dxd).
 *		J is stored sequentially as J = [ x_r x_s x_t y_r y_s y_t z_r z_s z_t ] (standard).
 *
 *		Investigation into the importance of proper treatment of the metric terms is still required. For this reason,
 *		the cofactor matrix terms are first projected to the basis of order PC and then projected to the VOLUME and
 *		FACET integration nodes. After the code is working, modify this by trying different combinations, such as
 *		directly computing C_vI, C_fI and make conclusions. (ToBeDeleted)
 *		detJ is computed directly at the VOLUME integration nodes. Thus, likely don't need PJs/PJc in setup_parameters
 *		=> Delete later (ToBeDeleted).
 *
 *	Notation:
 *		J    : (J)acobian of VOLUME coordinates (XYZ)
 *		C    : (C)ofactor terms of the Jacobian matrix
 *		detJ : (det)erminant of the Jacobian matrix
 *
 *	References:
 *		Kopriva(2006)-Metric_Identities_and_the_Discontinuous_Spectral_Element_Method_on_Curvilinear_Meshes
 */

struct S_OPERATORS {
	unsigned int NvnG, NvnC, NvnI, NfnI;
	double       *IC, *I_vG_vC, *I_vG_vI, *I_vC_vI,
	             **D_vG_vC, **D_vG_vI, **D_vC_vC, ***D_vG_fI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndFType);

void setup_geom_factors(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int i, n, row, col,
	             NvnG0, NvnC0, NvnI0;
	double       *XYZ, *J_vI, *J_vC, *detJV_vI, *C_vC, *C_vI;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	// Obtain operators
	init_ops(OPS,VOLUME,NULL,0);

	NvnG0 = OPS->NvnG;
	NvnC0 = OPS->NvnC;
	NvnI0 = OPS->NvnI;

	XYZ = VOLUME->XYZ;

	J_vI     = malloc(NvnI0*d*d * sizeof *J_vI);     // free
	J_vC     = malloc(NvnC0*d*d * sizeof *J_vC);     // free
	detJV_vI = malloc(NvnI0     * sizeof *detJV_vI); // keep
	C_vC     = malloc(NvnC0*d*d * sizeof *C_vC);     // keep (free after setup_normals)
	C_vI     = malloc(NvnI0*d*d * sizeof *C_vI);     // keep

	for (row = 0; row < d; row++) {
	for (col = 0; col < d; col++) {
		mm_CTN_d(NvnI0,1,NvnG0,OPS->D_vG_vI[col],&XYZ[NvnG0*row],&J_vI[NvnI0*(d*row+col)]);
		mm_CTN_d(NvnC0,1,NvnG0,OPS->D_vG_vC[col],&XYZ[NvnG0*row],&J_vC[NvnC0*(d*row+col)]);
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
			C_vC[NvnC0*(0+d*0)+n] =  J_vC[NvnC0*(d*1+1)+n]; // C11
			C_vC[NvnC0*(1+d*0)+n] = -J_vC[NvnC0*(d*0+1)+n]; // C21
			C_vC[NvnC0*(0+d*1)+n] = -J_vC[NvnC0*(d*1+0)+n]; // C12
			C_vC[NvnC0*(1+d*1)+n] =  J_vC[NvnC0*(d*0+0)+n]; // C22
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
			C_vC[NvnC0*(0+d*0)+n] = J_vC[NvnC0*(d*1+1)+n]*J_vC[NvnC0*(d*2+2)+n]
			                       -J_vC[NvnC0*(d*1+2)+n]*J_vC[NvnC0*(d*2+1)+n]; // C11
			C_vC[NvnC0*(0+d*1)+n] = J_vC[NvnC0*(d*1+2)+n]*J_vC[NvnC0*(d*2+0)+n]
			                       -J_vC[NvnC0*(d*1+0)+n]*J_vC[NvnC0*(d*2+2)+n]; // C12
			C_vC[NvnC0*(0+d*2)+n] = J_vC[NvnC0*(d*1+0)+n]*J_vC[NvnC0*(d*2+1)+n]
			                       -J_vC[NvnC0*(d*1+1)+n]*J_vC[NvnC0*(d*2+0)+n]; // C13
			C_vC[NvnC0*(1+d*0)+n] = J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*2+1)+n]
			                       -J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*2+2)+n]; // C21
			C_vC[NvnC0*(1+d*1)+n] = J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*2+2)+n]
			                       -J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*2+0)+n]; // C22
			C_vC[NvnC0*(1+d*2)+n] = J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*2+0)+n]
			                       -J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*2+1)+n]; // C23
			C_vC[NvnC0*(2+d*0)+n] = J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*1+2)+n]
			                       -J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*1+1)+n]; // C31
			C_vC[NvnC0*(2+d*1)+n] = J_vC[NvnC0*(d*0+2)+n]*J_vC[NvnC0*(d*1+0)+n]
			                       -J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*1+2)+n]; // C32
			C_vC[NvnC0*(2+d*2)+n] = J_vC[NvnC0*(d*0+0)+n]*J_vC[NvnC0*(d*1+1)+n]
			                       -J_vC[NvnC0*(d*0+1)+n]*J_vC[NvnC0*(d*1+0)+n]; // C33
		}
*/

		// curl form
		unsigned int dim, IndCurl, IndXYZJ, IndCurlXYZJ;
		unsigned int OrderCurl[12] = {1, 2, 2, 1, 2, 0, 0, 2, 0, 1, 1, 0};
		double *XYZJ_vC, *XYZ_vC, *tmp_vC, *CurlXYZJ_vC;

		XYZ_vC = malloc(d*NvnC0 * sizeof *XYZ_vC); // free
		XYZJ_vC = malloc(6*3*NvnC0 * sizeof *XYZJ_vC); // free

		mm_CTN_d(NvnC0,d,NvnG0,OPS->I_vG_vC,XYZ,XYZ_vC);


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
				mm_CTN_d(NvnC0,1,NvnC0,OPS->D_vC_vC[OrderCurl[IndCurl+0]],
				         &XYZJ_vC[NvnC0*(i*3+OrderCurl[IndCurl+2])],&CurlXYZJ_vC[NvnC0*IndCurlXYZJ]);
				mm_CTN_d(NvnC0,1,NvnC0,OPS->D_vC_vC[OrderCurl[IndCurl+1]],
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
					C_vC[NvnC0*(row+d*col)+n] =
						0.5*(CurlXYZJ_vC[NvnC0*IndCurlXYZJ+n]-CurlXYZJ_vC[NvnC0*(IndCurlXYZJ+3)+n]);
				}
				IndCurlXYZJ++;
			}
			IndCurlXYZJ += 3;
		}

//array_print_d(NvnC0,9,C_vC,'C');

		free(CurlXYZJ_vC);
	}

	mm_CTN_d(NvnI0,d*d,NvnC0,OPS->I_vC_vI,C_vC,C_vI);

	free(J_vI);
	free(J_vC);
	free(OPS);

//printf("Exiting setup_geom_factors.\n"), exit(1);

	VOLUME->detJV_vI = detJV_vI;
	if (VOLUME->C_vC)
		free(VOLUME->C_vC);
	VOLUME->C_vC     = C_vC;
	VOLUME->C_vI     = C_vI;

}

static void compute_detJV(const unsigned int Nn, double *J, double *detJV)
{
	/*
	 *	Comments:
	 *		Consider implementing the symmetric conservative form from:
	 *		Abe(2016)-Conservative_high-order_flux-reconstruction_schemes_on_moving_and_deforming_grids
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n;

	if (d == 1) {
		for (n = 0; n < Nn; n++) {
			detJV[n] = J[n];
		}
	} else if (d == 2) {
		for (n = 0; n < Nn; n++) {
			detJV[n] =   J[Nn*(d*0+0)+n]*J[Nn*(d*1+1)+n]
			           - J[Nn*(d*0+1)+n]*J[Nn*(d*1+0)+n];
		}
	} else if (d == 3) {
		for (n = 0; n < Nn; n++) {
			detJV[n] =   J[Nn*(d*0+0)+n]*(  J[Nn*(d*1+1)+n]*J[Nn*(d*2+2)+n]
			                              - J[Nn*(d*1+2)+n]*J[Nn*(d*2+1)+n])
			           - J[Nn*(d*0+1)+n]*(  J[Nn*(d*1+0)+n]*J[Nn*(d*2+2)+n]
			                              - J[Nn*(d*1+2)+n]*J[Nn*(d*2+0)+n])
			           + J[Nn*(d*0+2)+n]*(  J[Nn*(d*1+0)+n]*J[Nn*(d*2+1)+n]
			                              - J[Nn*(d*1+1)+n]*J[Nn*(d*2+0)+n]);
		}
	}
}

void setup_geom_factors_highorder(struct S_FACET *FACET)
{
	/*
	 *	Purpose:
	 *		Compute detJV_fI from each VOLUME.
	 *
	 *	Comments:
	 *		detJV_fI is only used for computing gradients at FACET nodes and thus not required for systems of 1st order
	 *		equations (such as when solving only the Euler equations).
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int row, col, IndFType,
	             NvnG, NfnI,
	             Vf, f, Eclass;
	double       *XYZ, *J_fI, *detJV_fI;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	// Obtain operators
	VOLUME = FACET->VIn;
	Vf     = FACET->VfIn;
	f      = Vf/NFREFMAX;

	Eclass = get_Eclass(VOLUME->type);
	IndFType = get_IndFType(Eclass,f);

	init_ops(OPS,VOLUME,FACET,IndFType);

	NvnG = OPS->NvnG;
	NfnI = OPS->NfnI;

	XYZ = VOLUME->XYZ;

	J_fI = malloc(NfnI*d*d * sizeof *J_fI); // free
	for (row = 0; row < d; row++) {
	for (col = 0; col < d; col++) {
		mm_CTN_d(NfnI,1,NvnG,OPS->D_vG_fI[Vf][col],&XYZ[NvnG*row],&J_fI[NfnI*(d*row+col)]);
	}}

	FACET->detJVIn_fI = malloc(NfnI * sizeof *detJV_fI); // keep
	compute_detJV(NfnI,J_fI,FACET->detJVIn_fI);

	if (!FACET->Boundary) {
		VOLUME = FACET->VOut;
		Vf     = FACET->VfOut;
		f      = Vf/NFREFMAX;

		init_ops(OPS,VOLUME,FACET,IndFType);

		NvnG = OPS->NvnG;

		XYZ = VOLUME->XYZ;

		for (row = 0; row < d; row++) {
		for (col = 0; col < d; col++) {
			mm_CTN_d(NfnI,1,NvnG,OPS->D_vG_fI[Vf][col],&XYZ[NvnG*row],&J_fI[NfnI*(d*row+col)]);
		}}

		FACET->detJVOut_fI = malloc(NfnI * sizeof *detJV_fI); // keep
		compute_detJV(NfnI,J_fI,FACET->detJVOut_fI);
	}
	free(J_fI);

	free(OPS);
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndFType)
{
	// Standard datatypes
	unsigned int P, type, curved;

	struct S_ELEMENT *ELEMENT_OPS;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT_OPS = get_ELEMENT_type(type);

	if (!curved) {
		OPS->NvnG = ELEMENT_OPS->NvnGs[1];
		OPS->NvnC = ELEMENT_OPS->NvnCs[P];
		OPS->NvnI = ELEMENT_OPS->NvnIs[P];

		OPS->I_vG_vC = ELEMENT_OPS->I_vGs_vCs[1][P][0];
		OPS->I_vC_vI = ELEMENT_OPS->I_vCs_vIs[P][P][0];
		OPS->D_vG_vC = ELEMENT_OPS->D_vGs_vCs[1][P][0];
		OPS->D_vG_vI = ELEMENT_OPS->D_vGs_vIs[1][P][0];
		OPS->D_vC_vC = ELEMENT_OPS->D_vCs_vCs[P][P][0];
	} else {
		OPS->NvnG = ELEMENT_OPS->NvnGc[P];
		OPS->NvnC = ELEMENT_OPS->NvnCc[P];
		OPS->NvnI = ELEMENT_OPS->NvnIc[P];

		OPS->I_vG_vC = ELEMENT_OPS->I_vGc_vCc[P][P][0];
		OPS->I_vC_vI = ELEMENT_OPS->I_vCc_vIc[P][P][0];
		OPS->D_vG_vC = ELEMENT_OPS->D_vGc_vCc[P][P][0];
		OPS->D_vG_vI = ELEMENT_OPS->D_vGc_vIc[P][P][0];
		OPS->D_vC_vC = ELEMENT_OPS->D_vCc_vCc[P][P][0];
	}

	if (FACET) {
		unsigned int PF, PV, FtypeInt;

		PV = P;
		PF = FACET->P;

		FtypeInt = FACET->typeInt;

		if (!curved) {
			// Straight VOLUME
			if (FtypeInt == 's') {
				// Straight FACET Integration
				OPS->NfnI    = ELEMENT_OPS->NfnIs[PF][IndFType];
				OPS->D_vG_fI = ELEMENT_OPS->D_vGs_fIs[1][PF];
			} else {
				// Curved FACET Integration
				OPS->NfnI    = ELEMENT_OPS->NfnIc[PF][IndFType];
				OPS->D_vG_fI = ELEMENT_OPS->D_vGs_fIc[1][PF];
			}
		} else {
			// Curved VOLUME
			if (FtypeInt == 's') {
				OPS->NfnI    = ELEMENT_OPS->NfnIs[PF][IndFType];
				OPS->D_vG_fI = ELEMENT_OPS->D_vGc_fIs[PV][PF];
			} else {
				OPS->NfnI    = ELEMENT_OPS->NfnIc[PF][IndFType];
				OPS->D_vG_fI = ELEMENT_OPS->D_vGc_fIc[PV][PF];
			}
		}
	}
}
