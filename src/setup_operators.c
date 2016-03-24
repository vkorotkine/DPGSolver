#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"

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
	int    i, j, count, dim, f, P,
	       Nf,
	       *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnS, *NvnF, *NvnFrs, *NvnFrc, *NvnIs, *NvnIc, *NvnP,
	       *NfnGc, *NfnIs, *NfnIc,
	       *ToReturn, *dummyi,
	       **Con_xir_vP;
	double Theta_eta[6], Theta_zeta[6],
	       *nr,*dummyd,
	       **xir_vGs, **xir_vGc, **xir_vCs, **xir_vCc, **xir_vJs, **xir_vJc, **xir_vS, **xir_vF, **xir_vFrs, **xir_vFrc,
	       **xir_vIs, **xir_vIc, **xir_vP, **WvIs, **WvIc,
	       ***xir_fGc, ***xir_fIs, ***xir_fIc, **WfIs, **WfIc;

	struct S_ELEMENT *ELEMENT;


	// Tensor-Product Operators
	ELEMENT = DB.ELEMENT; while(ELEMENT->type != LINE) ELEMENT = ELEMENT->next;
	dim = ELEMENT->d;

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
	ToReturn = malloc(4 * sizeof *ToReturn); // free

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





	// VOLUME Nodes (Order Independent)
	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&xir_vGs[0],&dummyd,&dummyi,&NvnGs[0],ToReturn,PGs,dim,"GLL");

	ToReturn[2] = 1;
	cubature_TP(&xir_vP[0],&dummyd,&Con_xir_vP[0],&NvnP[0],ToReturn,PP,d,"ES");
// array_print_d(NvnP[0],d,xir_vP[0]);
// array_print_i(pow(PP,d),pow(2,d),Con_xir_vP[0]);

    for (P = 0; P <= PMax; P++) {

		// VOLUME Nodes (Order dependent)
		ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
		cubature_TP(&xir_vGc[P],&dummyd,&dummyi,&NvnGc[P],ToReturn,PGc[P]   ,dim,"GLL");
		cubature_TP(&xir_vCs[P],&dummyd,&dummyi,&NvnCs[P],ToReturn,PCs[P][0],dim,"GLL");
		cubature_TP(&xir_vCc[P],&dummyd,&dummyi,&NvnCc[P],ToReturn,PCc[P][0],dim,"GLL");
		cubature_TP(&xir_vJs[P],&dummyd,&dummyi,&NvnCs[P],ToReturn,PJs[P][0],dim,"GLL");
		cubature_TP(&xir_vJc[P],&dummyd,&dummyi,&NvnCc[P],ToReturn,PJc[P][0],dim,"GLL");

		cubature_TP(&xir_vS[P]  ,&dummyd,&dummyi,&NvnS[P]  ,ToReturn,P         ,dim,NodeTypeS[P][0]);
		cubature_TP(&xir_vF[P]  ,&dummyd,&dummyi,&NvnF[P]  ,ToReturn,PF[P]     ,dim,NodeTypeF[P][0]);
		cubature_TP(&xir_vFrs[P],&dummyd,&dummyi,&NvnFrs[P],ToReturn,PFrs[P][0],dim,NodeTypeFrs[P][0]);
		cubature_TP(&xir_vFrc[P],&dummyd,&dummyi,&NvnFrc[P],ToReturn,PFrc[P][0],dim,NodeTypeFrc[P][0]);

		ToReturn[1] = 1;
		cubature_TP(&xir_vIs[P],&WvIs[P],&dummyi,&NvnIs[P],ToReturn,PIvs[P][0],dim,NodeTypeIvs[P][0]);
		cubature_TP(&xir_vIc[P],&WvIc[P],&dummyi,&NvnIc[P],ToReturn,PIvc[P][0],dim,NodeTypeIvc[P][0]);

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
			xir_fGc[P][f] = malloc(NfnGc[P]*dim * sizeof ***xir_fGc); // tbd
			xir_fIs[P][f] = malloc(NfnIs[P]*dim * sizeof ***xir_fIs); // tbd
			xir_fIc[P][f] = malloc(NfnIc[P]*dim * sizeof ***xir_fIc); // tbd

			xir_fGc[P][f][0] = pow(-1,f+1);
			xir_fIs[P][f][0] = pow(-1,f+1);
			xir_fIc[P][f][0] = pow(-1,f+1);
		}

		// Preliminary Operators

	}

	free(ToReturn);


	// Implementing basis_TP

	int P_tmp, Nn_tmp, dim_tmp, d_tmp,
	    *dummyi_tmp, ToReturn_tmp[4];
	double *xir_tmp, *dummyd_tmp, *ChiRef_tmp, **GradChiRef_tmp;

	P_tmp = 2;
	dim_tmp = 2;

	ToReturn_tmp[0] = 1; ToReturn_tmp[1] = 0; ToReturn_tmp[2] = 0; ToReturn_tmp[3] = 1;
	cubature_TP(&xir_tmp,&dummyd_tmp,&dummyi_tmp,&Nn_tmp,ToReturn_tmp,P_tmp,dim_tmp,"GL");
	ChiRef_tmp = basis_TP(P_tmp,xir_tmp,Nn_tmp,dim_tmp);
	GradChiRef_tmp = grad_basis_TP(P_tmp,xir_tmp,Nn_tmp,dim_tmp);

	array_print_d(Nn_tmp,pow(P_tmp+1,dim_tmp),ChiRef_tmp);

	for (d_tmp = 0; d_tmp < dim_tmp; d_tmp++)
		array_print_d(Nn_tmp,pow(P_tmp+1,dim_tmp),GradChiRef_tmp[d_tmp]);

	free(xir_tmp);
	free(ChiRef_tmp);
	array_free2_d(dim_tmp,GradChiRef_tmp);






}
