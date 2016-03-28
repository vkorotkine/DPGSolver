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

void setup_operators()
{
	// Initialize DB Parameters
	int  d          = DB.d,
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
	     **SF_BE    = DB.SF_BE,

	     Testing    = DB.Testing;

	char *BasisType     = DB.BasisType,
	     ***NodeTypeS   = DB.NodeTypeS,
	     ***NodeTypeF   = DB.NodeTypeF,
	     ***NodeTypeFrs = DB.NodeTypeFrs,
	     ***NodeTypeFrc = DB.NodeTypeFrc,
	     ***NodeTypeIfs = DB.NodeTypeIfs,
	     ***NodeTypeIfc = DB.NodeTypeIfc,
	     ***NodeTypeIvs = DB.NodeTypeIvs,
	     ***NodeTypeIvc = DB.NodeTypeIvc;

	int  PrintTesting = 0, MPIrank = DB.MPIrank;

	// Standard datatypes
	int    i, j, count, dE, f, P,
	       Nf,
	       *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnS, *NvnF, *NvnFrs, *NvnFrc, *NvnIs, *NvnIc, *NvnP,
	       *NfnGc, *NfnIs, *NfnIc,
	       *ToReturn, *dummyCon,
	       **Con_xir_vP;
	double Theta_eta[6], Theta_zeta[6],
	       *nr,*dummyW, *dummyxir,
	       **xir_vGs, **xir_vGc, **xir_vCs, **xir_vCc, **xir_vJs, **xir_vJc, **xir_vS, **xir_vF, **xir_vFrs, **xir_vFrc,
	       **xir_vIs, **xir_vIc, **xir_vP, **WvIs, **WvIc,
	       ***xir_fGc, ***xir_fIs, ***xir_fIc, **WfIs, **WfIc,
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
	xir_vGs  = ELEMENT->xir_vGs ; NvnGs  = ELEMENT->NvnGs;
	xir_vGc  = ELEMENT->xir_vGc ; NvnGc  = ELEMENT->NvnGc;
	xir_vCs  = ELEMENT->xir_vCs ; NvnCs  = ELEMENT->NvnCs;
	xir_vCc  = ELEMENT->xir_vCc ; NvnCc  = ELEMENT->NvnCc;
	xir_vJs  = ELEMENT->xir_vJs ; NvnJs  = ELEMENT->NvnJs;
	xir_vJc  = ELEMENT->xir_vJc ; NvnJc  = ELEMENT->NvnJc;
	xir_vS   = ELEMENT->xir_vS  ; NvnS   = ELEMENT->NvnS;
	xir_vF   = ELEMENT->xir_vF  ; NvnF   = ELEMENT->NvnF;
	xir_vFrs = ELEMENT->xir_vFrs; NvnFrs = ELEMENT->NvnFrs;
	xir_vFrc = ELEMENT->xir_vFrc; NvnFrc = ELEMENT->NvnFrc;
	xir_vIs  = ELEMENT->xir_vIs ; NvnIs  = ELEMENT->NvnIs ; WvIs = ELEMENT->WvIs;
	xir_vIc  = ELEMENT->xir_vIc ; NvnIc  = ELEMENT->NvnIc ; WvIc = ELEMENT->WvIc;
	xir_vP   = ELEMENT->xir_vP  ; NvnP   = ELEMENT->NvnP  ; Con_xir_vP = ELEMENT->Con_xir_vP;

	// FACET Nodes
	xir_fGc = ELEMENT->xir_fGc; NfnGc = ELEMENT->NfnGc;
	xir_fIs = ELEMENT->xir_fIs; NfnIs = ELEMENT->NfnIs; WfIs = ELEMENT->WfIs;
	xir_fIc = ELEMENT->xir_fIc; NfnIc = ELEMENT->NfnIc; WfIc = ELEMENT->WfIc;

	// Preliminary Operators
	ChiRefGs_vGc = malloc(NP * sizeof *ChiRefGs_vGc); // tbd
	ChiGs_vGc = malloc(NP * sizeof *ChiGs_vGc); // tbd

	// Operators
	I_vGs_vGc = ELEMENT->I_vGs_vGc;




	// VOLUME Nodes (Order Independent)
	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&xir_vGs[0],&dummyW,&dummyCon,&NvnGs[0],ToReturn,PGs,dE,"GLL");

	ToReturn[2] = 1;
	cubature_TP(&xir_vP[0],&dummyW,&Con_xir_vP[0],&NvnP[0],ToReturn,PP,d,"ES");

	// Preliminary Operators
	IGs = identity_d(NvnGs[0]); // tbd

	ChiRefGs_vGs = basis_TP(PGs,xir_vGs[0],NvnGs[0],dE); // tbd

	if (strstr(BasisType,"Modal") != NULL) {
		ChiGs_vGs = ChiRefGs_vGs;
	} else if (strstr(BasisType,"Nodal") != NULL) {
		ChiGs_vGs = IGs;
	}

	ChiRefInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiRefGs_vGs,IGs); // tbd

	ChiInvGs_vGs    = inverse_d(NvnGs[0],NvnGs[0],ChiGs_vGs,IGs); // tbd

	TGs = mm_d(CblasNoTrans,CblasNoTrans,NvnGs[0],NvnGs[0],NvnGs[0],1.0,ChiRefInvGs_vGs,ChiGs_vGs); // tbd


/*
array_print_d(NvnGs[0],NvnGs[0],ChiRefGs_vGs);
array_print_d(NvnGs[0],NvnGs[0],IGs);
array_print_d(NvnGs[0],NvnGs[0],ChiRefInvGs_vGs);
array_print_d(NvnGs[0],NvnGs[0],TGs);
*/


    for (P = 0; P <= PMax; P++) {

		// VOLUME Nodes (Order dependent)
		ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
		cubature_TP(&xir_vGc[P],&dummyW,&dummyCon,&NvnGc[P],ToReturn,PGc[P]   ,dE,"GLL");
		cubature_TP(&xir_vCs[P],&dummyW,&dummyCon,&NvnCs[P],ToReturn,PCs[P][0],dE,"GLL");
		cubature_TP(&xir_vCc[P],&dummyW,&dummyCon,&NvnCc[P],ToReturn,PCc[P][0],dE,"GLL");
		cubature_TP(&xir_vJs[P],&dummyW,&dummyCon,&NvnCs[P],ToReturn,PJs[P][0],dE,"GLL");
		cubature_TP(&xir_vJc[P],&dummyW,&dummyCon,&NvnCc[P],ToReturn,PJc[P][0],dE,"GLL");

		cubature_TP(&xir_vS[P]  ,&dummyW,&dummyCon,&NvnS[P]  ,ToReturn,P         ,dE,NodeTypeS[P][0]);
		cubature_TP(&xir_vF[P]  ,&dummyW,&dummyCon,&NvnF[P]  ,ToReturn,PF[P]     ,dE,NodeTypeF[P][0]);
		cubature_TP(&xir_vFrs[P],&dummyW,&dummyCon,&NvnFrs[P],ToReturn,PFrs[P][0],dE,NodeTypeFrs[P][0]);
		cubature_TP(&xir_vFrc[P],&dummyW,&dummyCon,&NvnFrc[P],ToReturn,PFrc[P][0],dE,NodeTypeFrc[P][0]);

		ToReturn[1] = 1;
		cubature_TP(&xir_vIs[P],&WvIs[P],&dummyCon,&NvnIs[P],ToReturn,PIvs[P][0],dE,NodeTypeIvs[P][0]);
		cubature_TP(&xir_vIc[P],&WvIc[P],&dummyCon,&NvnIc[P],ToReturn,PIvc[P][0],dE,NodeTypeIvc[P][0]);

		// FACET Nodes (Order dependent in the general case => redundant here)
		WfIs[P] = malloc(1 * sizeof **WfIs); // tbd
		WfIc[P] = malloc(1 * sizeof **WfIc); // tbd

		xir_fGc[P] = malloc(Nf * sizeof **xir_fGc); // tbd
		xir_fIs[P] = malloc(Nf * sizeof **xir_fIs); // tbd
		xir_fIc[P] = malloc(Nf * sizeof **xir_fIc); // tbd

		NfnGc[P] = 1;
		NfnIs[P] = 1; WfIs[P][0] = 1;
		NfnIc[P] = 1; WfIc[P][0] = 1;

		for (f = 0; f < Nf; f++) {
			xir_fGc[P][f] = malloc(NfnGc[P]*dE * sizeof ***xir_fGc); // tbd
			xir_fIs[P][f] = malloc(NfnIs[P]*dE * sizeof ***xir_fIs); // tbd
			xir_fIc[P][f] = malloc(NfnIc[P]*dE * sizeof ***xir_fIc); // tbd

			xir_fGc[P][f][0] = pow(-1,f+1);
			xir_fIs[P][f][0] = pow(-1,f+1);
			xir_fIc[P][f][0] = pow(-1,f+1);
		}

		// Preliminary Operators
		ChiRefGs_vGc[P] = basis_TP(PGs,xir_vGc[P],NvnGc[P],dE); // tbd

		ChiGs_vGc[P] = mm_d(CblasNoTrans,CblasNoTrans,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vGc[P],TGs); // tbd



		// Operators
		I_vGs_vGc[P] = mm_d(CblasNoTrans,CblasNoTrans,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vGc[P],ChiInvGs_vGs); // keep

array_print_d(NvnGc[P],NvnGs[0],I_vGs_vGc[P]);

	}



	// QUAD
	if (is_ELEMENT_present(QUAD)) {
		ELEMENT = get_ELEMENT_type(QUAD);

		dE = ELEMENT->d;

		NvnGs  = ELEMENT->NvnGs;

		// VOLUME Nodes (Order Independent)
		ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
		cubature_TP(&dummyxir,&dummyW,&dummyCon,&NvnGs[0],ToReturn,PGs,dE,"GLL");

	}

	// HEX
	if (is_ELEMENT_present(HEX)) {
		ELEMENT = get_ELEMENT_type(HEX);

		dE = ELEMENT->d;

		NvnGs  = ELEMENT->NvnGs;

		// VOLUME Nodes (Order Independent)
		ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
		cubature_TP(&dummyxir,&dummyW,&dummyCon,&NvnGs[0],ToReturn,PGs,dE,"GLL");
	}








	free(ToReturn);

}
