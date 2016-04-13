#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

#include "mkl.h"

// DON'T FORGET TO CHANGE COMMENTS (ToBeDeleted)

/*
 *	Purpose:
 *		Set up operators to be used throughout the code.
 *
 *	Comments:
 *		Currently, only tensor-product operators are assembled using the sum factorized 1D operators. (ToBeModified)
 *		Intuitively, the collapsed tensor-product elements seem to be less efficient than those developed for other
 *		element types.
 *
 *	Notation:
 *		Theta_[] : Angles for conversion between reference and facet coordinates (Zwanenburg(2016): Table 14)
 *		           Options: eta
 *		                    zeta
 *
 *	References:
 *		Zwanenburg(2016)_Equivalence between the Energy Stable Flux Reconstruction and Filtered Discontinuous Galerkin
 *		                 Schemes
 *		Karniadakis(1999)_Spectral hp Element Methods for CFD
 *		Hesthaven(2000)_Stable Spectral Methods on Tetrahedral Elements
 *		Stiller(2008)_Factorization Techniques for Nodal Spectral Elements in Curved Domains
*/

void setup_operators(void)
{
	// Initialize DB Parameters
	unsigned int d          = DB.d,
	             NP         = DB.NP,
	             PMax       = DB.PMax,
	             PGs        = DB.PGs,
	             *PGc       = DB.PGc,
	             **PCs      = DB.PCs,
	             **PCc      = DB.PCc,
	             **PJs      = DB.PJs,
	             **PJc      = DB.PJc,
	             *PF        = DB.PF,
	             **PFrs     = DB.PFrs,
	             **PFrc     = DB.PFrc,
	             **PIfs     = DB.PIfs,
	             **PIfc     = DB.PIfc,
	             **PIvs     = DB.PIvs,
	             **PIvc     = DB.PIvc,
	             PR         = DB.PR,
	             PP         = DB.PP,
	             Restart    = DB.Restart,
	             EFE        = DB.EFE,
	             Collocated = DB.Collocated,
	             **SF_BE    = DB.SF_BE;

	char         *BasisType     = DB.BasisType,
	             ***NodeTypeS   = DB.NodeTypeS,
	             ***NodeTypeF   = DB.NodeTypeF,
	             ***NodeTypeFrs = DB.NodeTypeFrs,
	             ***NodeTypeFrc = DB.NodeTypeFrc,
	             ***NodeTypeIfs = DB.NodeTypeIfs,
	             ***NodeTypeIfc = DB.NodeTypeIfc,
	             ***NodeTypeIvs = DB.NodeTypeIvs,
	             ***NodeTypeIvc = DB.NodeTypeIvc;

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int i, j, count, dE, f, P,
	             Nf,
	             *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnS, *NvnF, *NvnFrs, *NvnFrc, *NvnIs, *NvnIc, *NvnP,
	             *NfnGc, *NfnIs, *NfnIc,
	             *ToReturn, *dummyCon,
	             **Con_rst_vP;
	double       Theta_eta[6], Theta_zeta[6],
	             *nr,*dummyw, *dummyrst,
	             **rst_vGs, **rst_vGc, **rst_vCs, **rst_vCc, **rst_vJs, **rst_vJc, **rst_vS, **rst_vF, **rst_vFrs, **rst_vFrc,
	             **rst_vIs, **rst_vIc, **rst_vP, **wvIs, **wvIc,
	             ***rst_fGc, ***rst_fIs, ***rst_fIc, **wfIs, **wfIc,
	             *IGs,
	             *ChiRefGs_vGs,
	             *ChiGs_vGs,
	             *ChiRefInvGs_vGs,
	             *ChiInvGs_vGs,
	             *TGs,
	             **ChiRefGs_vGc,
	             **ChiGs_vGc,
	             **I_vGs_vGc;

	struct S_ELEMENT *ELEMENT;

	ToReturn = malloc(4 * sizeof *ToReturn); // free

	// LINE and TP Class
	ELEMENT = DB.ELEMENT; while(ELEMENT->type != LINE) ELEMENT = ELEMENT->next;
	dE = ELEMENT->d;

	// Set up Normals
	Nf = ELEMENT->Nf;

	Theta_eta[0]  = 0.; Theta_eta[1]  = 0.;
	Theta_zeta[0] = PI; Theta_zeta[1] = 0.;

	nr = ELEMENT->nr;
	for (f = 0; f < Nf; f++) {
		nr[f*d+0] = cos(Theta_eta[f])*cos(Theta_zeta[f]);
		if (d > 1) nr[f*d+1] = cos(Theta_eta[f])*sin(Theta_zeta[f]);
		if (d > 2) nr[f*d+2] = -sin(Theta_zeta[f]);
	}

	// VOLUME Nodes

	/* NOTE: NEED TO REMOVE SOME OF THESE FROM THE ELEMENT STRUCTURE, THEY ARE NOT ALL NECESSARY. FOR THOSE REMOVED,
	 *       ALLOCATE MEMORY AND FREE IN setup_operators. FROM NOW ON, ONLY INCLUDE OPERATORS IN THE ELEMENT STRUCTURE
	 *       IF IT IS NEEDED SOMEWHERE IN THE CODE.
	 */

	rst_vGs  = ELEMENT->rst_vGs ; NvnGs  = ELEMENT->NvnGs;
	rst_vGc  = ELEMENT->rst_vGc ; NvnGc  = ELEMENT->NvnGc;
	rst_vCs  = ELEMENT->rst_vCs ; NvnCs  = ELEMENT->NvnCs;
	rst_vCc  = ELEMENT->rst_vCc ; NvnCc  = ELEMENT->NvnCc;
	rst_vJs  = ELEMENT->rst_vJs ; NvnJs  = ELEMENT->NvnJs;
	rst_vJc  = ELEMENT->rst_vJc ; NvnJc  = ELEMENT->NvnJc;
	rst_vS   = ELEMENT->rst_vS  ; NvnS   = ELEMENT->NvnS;
	rst_vF   = ELEMENT->rst_vF  ; NvnF   = ELEMENT->NvnF;
	rst_vFrs = ELEMENT->rst_vFrs; NvnFrs = ELEMENT->NvnFrs;
	rst_vFrc = ELEMENT->rst_vFrc; NvnFrc = ELEMENT->NvnFrc;
	rst_vIs  = ELEMENT->rst_vIs ; NvnIs  = ELEMENT->NvnIs ; wvIs = ELEMENT->wvIs;
	rst_vIc  = ELEMENT->rst_vIc ; NvnIc  = ELEMENT->NvnIc ; wvIc = ELEMENT->wvIc;
	rst_vP   = ELEMENT->rst_vP  ; NvnP   = ELEMENT->NvnP  ; Con_rst_vP = ELEMENT->Con_rst_vP;

	// FACET Nodes
	rst_fGc = ELEMENT->rst_fGc; NfnGc = ELEMENT->NfnGc;
	rst_fIs = ELEMENT->rst_fIs; NfnIs = ELEMENT->NfnIs; wfIs = ELEMENT->wfIs;
	rst_fIc = ELEMENT->rst_fIc; NfnIc = ELEMENT->NfnIc; wfIc = ELEMENT->wfIc;

	// Preliminary Operators
	ChiRefGs_vGc = malloc(NP * sizeof *ChiRefGs_vGc); // tbd
	ChiGs_vGc = malloc(NP * sizeof *ChiGs_vGc); // tbd

	// Operators
	I_vGs_vGc = ELEMENT->I_vGs_vGc;




	// VOLUME Nodes (Order Independent)
	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst_vGs[0],&dummyw,&dummyCon,&NvnGs[0],ToReturn,PGs,dE,"GLL"); // tbd

	ToReturn[2] = 1;
	cubature_TP(&rst_vP[0],&dummyw,&Con_rst_vP[0],&NvnP[0],ToReturn,PP,d,"ES"); // tbd

	// Preliminary Operators
	IGs = identity_d(NvnGs[0]); // tbd

	ChiRefGs_vGs = basis_TP(PGs,rst_vGs[0],NvnGs[0],dE); // tbd

	if (strstr(BasisType,"Modal") != NULL) {
		ChiGs_vGs = ChiRefGs_vGs;
	} else if (strstr(BasisType,"Nodal") != NULL) {
		ChiGs_vGs = IGs;
	}

	ChiRefInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiRefGs_vGs,IGs); // tbd

	ChiInvGs_vGs    = inverse_d(NvnGs[0],NvnGs[0],ChiGs_vGs,IGs); // tbd

/// Add in CBLAS_LAYOUT parameter to mm_Alloc ///
	TGs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
	                 NvnGs[0],NvnGs[0],NvnGs[0],1.0,ChiRefInvGs_vGs,ChiGs_vGs); // tbd


/*
array_print_d(NvnGs[0],NvnGs[0],ChiRefGs_vGs,'R');
array_print_d(NvnGs[0],NvnGs[0],IGs,'R');
array_print_d(NvnGs[0],NvnGs[0],ChiRefInvGs_vGs,'R');
array_print_d(NvnGs[0],NvnGs[0],TGs,'R');
*/

    for (P = 0; P <= PMax; P++) {

		// VOLUME Nodes (Order dependent)
		ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
		cubature_TP(&rst_vGc[P],&dummyw,&dummyCon,&NvnGc[P],ToReturn,PGc[P]   ,dE,"GLL"); // tbd
		cubature_TP(&rst_vCs[P],&dummyw,&dummyCon,&NvnCs[P],ToReturn,PCs[P][0],dE,"GLL"); // tbd
		cubature_TP(&rst_vCc[P],&dummyw,&dummyCon,&NvnCc[P],ToReturn,PCc[P][0],dE,"GLL"); // tbd
		cubature_TP(&rst_vJs[P],&dummyw,&dummyCon,&NvnCs[P],ToReturn,PJs[P][0],dE,"GLL"); // tbd
		cubature_TP(&rst_vJc[P],&dummyw,&dummyCon,&NvnCc[P],ToReturn,PJc[P][0],dE,"GLL"); // tbd

		cubature_TP(&rst_vS[P]  ,&dummyw,&dummyCon,&NvnS[P]  ,ToReturn,P         ,dE,NodeTypeS[P][0]); // tbd
		cubature_TP(&rst_vF[P]  ,&dummyw,&dummyCon,&NvnF[P]  ,ToReturn,PF[P]     ,dE,NodeTypeF[P][0]); // tbd
		cubature_TP(&rst_vFrs[P],&dummyw,&dummyCon,&NvnFrs[P],ToReturn,PFrs[P][0],dE,NodeTypeFrs[P][0]); // tbd
		cubature_TP(&rst_vFrc[P],&dummyw,&dummyCon,&NvnFrc[P],ToReturn,PFrc[P][0],dE,NodeTypeFrc[P][0]); // tbd

		ToReturn[1] = 1;
		cubature_TP(&rst_vIs[P],&wvIs[P],&dummyCon,&NvnIs[P],ToReturn,PIvs[P][0],dE,NodeTypeIvs[P][0]); // tbd
		cubature_TP(&rst_vIc[P],&wvIc[P],&dummyCon,&NvnIc[P],ToReturn,PIvc[P][0],dE,NodeTypeIvc[P][0]); // tbd

		// FACET Nodes (Order dependent in the general case => redundant here)
		wfIs[P] = malloc(1 * sizeof **wfIs); // tbd
		wfIc[P] = malloc(1 * sizeof **wfIc); // tbd

		rst_fGc[P] = malloc(Nf * sizeof **rst_fGc); // tbd
		rst_fIs[P] = malloc(Nf * sizeof **rst_fIs); // tbd
		rst_fIc[P] = malloc(Nf * sizeof **rst_fIc); // tbd

		NfnGc[P] = 1;
		NfnIs[P] = 1; wfIs[P][0] = 1;
		NfnIc[P] = 1; wfIc[P][0] = 1;

		for (f = 0; f < Nf; f++) {
			rst_fGc[P][f] = malloc(NfnGc[P]*dE * sizeof ***rst_fGc); // tbd
			rst_fIs[P][f] = malloc(NfnIs[P]*dE * sizeof ***rst_fIs); // tbd
			rst_fIc[P][f] = malloc(NfnIc[P]*dE * sizeof ***rst_fIc); // tbd

			rst_fGc[P][f][0] = pow(-1,f+1);
			rst_fIs[P][f][0] = pow(-1,f+1);
			rst_fIc[P][f][0] = pow(-1,f+1);
		}

		// Preliminary Operators
		ChiRefGs_vGc[P] = basis_TP(PGs,rst_vGc[P],NvnGc[P],dE); // tbd

		ChiGs_vGc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                          NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vGc[P],TGs); // tbd



		// Operators
		I_vGs_vGc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                          NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vGc[P],ChiInvGs_vGs); // keep

//array_print_d(NvnGc[P],NvnGs[0],I_vGs_vGc[P],'R');

	}



	// QUAD
	if (is_ELEMENT_present(QUAD)) {
		ELEMENT = get_ELEMENT_type(QUAD);

		dE = ELEMENT->d;

		NvnGs  = ELEMENT->NvnGs;

		// VOLUME Nodes (Order Independent)
		ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
		cubature_TP(&dummyrst,&dummyw,&dummyCon,&NvnGs[0],ToReturn,PGs,dE,"GLL"); // tbd

	}

	// HEX
	if (is_ELEMENT_present(HEX)) {
		ELEMENT = get_ELEMENT_type(HEX);

		dE = ELEMENT->d;

		NvnGs  = ELEMENT->NvnGs;

		// VOLUME Nodes (Order Independent)
		ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
		cubature_TP(&dummyrst,&dummyw,&dummyCon,&NvnGs[0],ToReturn,PGs,dE,"GLL"); // tbd
	}








	free(ToReturn);

}
