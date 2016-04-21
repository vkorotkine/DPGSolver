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
 *		Standard (i.e. non sum factorized) operators are only computed for performance critical functions; in all other
 *		functions, the sum factorized operators are used even if they are less efficient than the application of the
 *		standard operator.
 *		Currently, only tensor-product operators are assembled using the sum factorized 1D operators. (ToBeModified)
 *		Intuitively, the collapsed tensor-product elements seem to be less efficient than those developed for other
 *		element types.
 *		Ensure that operators for hp refinement are only stored when refinement is enabled (ToBeDeleted).
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

typedef void (*cubature_tdef) (double **rst, double **w_vec, unsigned int *Nn,const unsigned int return_w,
                               const unsigned int P, const unsigned int d, const char *NodeType);
typedef double *(*basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                               const unsigned int d);
typedef double **(*grad_basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn,
                                     unsigned int *NbfOut, const unsigned int d);

static void select_functions(basis_tdef *basis, grad_basis_tdef *grad_basis, cubature_tdef *cubature,
                             const unsigned int type)
{
	printf("In select functions\n");
	if (type == LINE) {
		*basis      = basis_TP;
		*grad_basis = grad_basis_TP;
		*cubature   = cubature_TP;
	} else if (type == TRI) {
		printf("Add support: select functions (TRI)\n");
	} else if (type == TET) {
		printf("Add support: select functions (TET)\n");
	} else if (type == WEDGE) {
		printf("Add support: select functions (WEDGE)\n");
	} else if (type == PYR) {
		printf("Add support: select functions (PYR)\n");
	} else {
		printf("Error: Unsupported type in select_functions.\n"), exit(1);
	}
}

static void setup_ELEMENT_plotting(const unsigned int EType)
{
	// Initialize DB Parameters
	unsigned int PP = DB.PP;

	// Standard datatypes
	unsigned int NvnP, NE,
	             *connectivity, *types;
	double       *rst_vP;
	
	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	plotting_element_info(&rst_vP,&connectivity,&types,&NvnP,&NE,PP,EType); // free
	ELEMENT->connectivity  = connectivity;
	ELEMENT->connect_types = types;
	ELEMENT->connect_NE    = NE;
	ELEMENT->NvnP = NvnP;

	free(rst_vP);
}

static void setup_ELEMENT_normals(const unsigned int EType)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes 
	unsigned int f, Nf;
	double Theta_eta[6], Theta_zeta[6], *nr;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	Nf = ELEMENT->Nf;
	nr = ELEMENT->nr;

	if (EType == LINE) {
		Theta_eta[0] = 0.0; Theta_zeta[0] = PI ;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 0.0;
	} else if (EType == TRI) {
		Theta_eta[0] = 0.0; Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = 0.0; Theta_zeta[2] = 9.0/6.0*PI;
	} else if (EType == TET) {
		double Theta_e = atan(2.0*sqrt(2.0));

		Theta_eta[0] = Theta_e - PI/2.0; Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = Theta_e - PI/2.0; Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = Theta_e - PI/2.0; Theta_zeta[2] = 9.0/6.0*PI;
		Theta_eta[3] = PI/2.0          ; Theta_zeta[3] = 0.0       ;
	} else {
		printf("Add support setup_ELEMENT_normals\n");
		exit(1);
	}

	for (f = 0; f < Nf; f++) {
		nr[f*d+0] = cos(Theta_eta[f])*cos(Theta_zeta[f]);
		if (d > 1) nr[f*d+1] = cos(Theta_eta[f])*sin(Theta_zeta[f]);
		if (d > 2) nr[f*d+2] = -sin(Theta_zeta[f]);
	}
}

static void setup_ELEMENT_operators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, *NvnGc;
	double       **I_vGs_vGc, **I_vGs_vP, **I_vGc_vP;

	// Initialize DB Parameters
	unsigned int PMax        = DB.PMax,
	             PGs         = DB.PGs,
	             *PGc        = DB.PGc,
	             PP          = DB.PP;

	char         *BasisType  = DB.BasisType;

	// Standard datatypes
	unsigned int dE, P,
	             Nbf,
	             NvnP,
	             dummy_ui, *dummyPtr_ui[2];
	double       *rst_vP, *rst_vGs, *rst_vGc,
	             *IGs, *IGc,
	             *TGs, *TGc,
	             *ChiRefGs_vGs, *ChiGs_vGs, *ChiRefGc_vGc, *ChiGc_vGc,
	             *ChiRefInvGs_vGs, *ChiInvGs_vGs, *ChiRefInvGc_vGc, *ChiInvGc_vGc,
	             *ChiRefGs_vP, *ChiGs_vP, *ChiRefGs_vGc, *ChiGs_vGc, *ChiRefGc_vP, *ChiGc_vP,
	             *dummyPtr_d;

	struct S_ELEMENT *ELEMENT;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	ELEMENT = get_ELEMENT_type(EType);

	dE = ELEMENT->d;

	select_functions(&basis,&grad_basis,&cubature,EType);


	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnGc = ELEMENT->NvnGc;

	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGc_vP  = ELEMENT->I_vGc_vP;



	// VOLUME Nodes (Order Independent)
	plotting_element_info(&rst_vP,&dummyPtr_ui[0],&dummyPtr_ui[1],&NvnP,&dummy_ui,PP,LINE); // free
	free(dummyPtr_ui[0]);
	free(dummyPtr_ui[1]);

	cubature(&rst_vGs,&dummyPtr_d,&NvnGs[0],0,PGs,dE,"GLL"); // free

	// Preliminary Operators
	IGs = identity_d(NvnGs[0]); // free

	ChiRefGs_vGs = basis(PGs,rst_vGs,NvnGs[0],&Nbf,dE); // free
	

	if (strstr(BasisType,"Modal") != NULL) {
		ChiGs_vGs = ChiRefGs_vGs;
	} else if (strstr(BasisType,"Nodal") != NULL) {
		ChiGs_vGs = IGs;
	}

	ChiRefInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiRefGs_vGs,IGs); // free

	ChiInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiGs_vGs,IGs); // free

	TGs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
	                 NvnGs[0],NvnGs[0],NvnGs[0],1.0,ChiRefInvGs_vGs,ChiGs_vGs); // free

	ChiRefGs_vP = basis(PGs,rst_vP,NvnP,&Nbf,dE); // free

	ChiGs_vP = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,NvnGs[0],NvnGs[0],1.0,ChiRefGs_vP,TGs); // free

	// Returned Operators
	I_vGs_vP[0] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,NvnGs[0],NvnGs[0],1.0,ChiGs_vP,ChiInvGs_vGs); // keep

	free(rst_vGs);

	free(IGs);
	free(ChiRefGs_vGs);
	free(ChiRefInvGs_vGs);
	free(ChiRefGs_vP);
	free(ChiGs_vP);

	for (P = 0; P <= PMax; P++) {
		cubature(&rst_vGc,&dummyPtr_d,&NvnGc[P],0,PGc[P]   ,dE,"GLL"); // free

		// Preliminary Operators
		IGc = identity_d(NvnGc[P]); // free

		ChiRefGc_vGc = basis(PGc[P],rst_vGc,NvnGc[P],&Nbf,dE); // free

		if (strstr(BasisType,"Modal") != NULL) {
			ChiGc_vGc = ChiRefGc_vGc;
		} else if (strstr(BasisType,"Nodal") != NULL) {
			ChiGc_vGc = IGc;
		}

		ChiRefInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiRefGc_vGc,IGc); // free

		ChiInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiGc_vGc,IGc); // free

		TGc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                 NvnGc[P],NvnGc[P],NvnGc[P],1.0,ChiRefInvGc_vGc,ChiGc_vGc); // free


		ChiRefGs_vGc = basis(PGs,rst_vGc,NvnGc[P],&Nbf,dE); // free
		ChiRefGc_vP  = basis(PGc[P],rst_vP,NvnP,&Nbf,dE);   // free

		ChiGs_vGc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                       NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vGc,TGs); // free
		ChiGc_vP  = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                       NvnP,NvnGc[P],NvnGc[P],1.0,ChiRefGc_vP,TGc);      // free

		// Returned Operators
		I_vGs_vGc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                          NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vGc,ChiInvGs_vGs); // keep

		I_vGc_vP[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                         NvnP,NvnGc[P],NvnGc[P],1.0,ChiGc_vP,ChiInvGc_vGc);    // keep

		free(rst_vGc);

		free(ChiRefGs_vGc);
		free(ChiGs_vGc);

		free(IGc);
		free(ChiRefGc_vGc);
		free(ChiRefInvGc_vGc);
		free(ChiInvGc_vGc);
		free(TGc);
		free(ChiRefGc_vP);
		free(ChiGc_vP);
	}

	free(rst_vP);

	free(ChiInvGs_vGs);
	free(TGs);
}

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
	             EType, Nf, NE, Nbf,
	             *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnS, *NvnF, *NvnFrs, *NvnFrc, *NvnIs, *NvnIc, NvnP,
	             *NfnGc, *NfnIs, *NfnIc,
	             *connectivity, *types;
	double       Theta_eta[6], Theta_zeta[6],
	             *nr,*dummyw, *dummyrst,
	             *rst_vGs, *rst_vGc, **rst_vCs, **rst_vCc, **rst_vJs, **rst_vJc, **rst_vS, **rst_vF, **rst_vFrs, **rst_vFrc,
	             **rst_vIs, **rst_vIc, *rst_vP, **wvIs, **wvIc,
	             ***rst_fGc, ***rst_fIs, ***rst_fIc, **wfIs, **wfIc,
	             *IP, *IGs, *IGc,
	             *ChiRefGs_vGs, *ChiRefGc_vGc,
	             *ChiGs_vGs, *ChiGc_vGc,
	             *ChiRefInvGs_vGs, *ChiRefInvGc_vGc,
	             *ChiInvGs_vGs, *ChiInvGc_vGc,
	             *TGs, *TGc,
	             *ChiRefGs_vGc, *ChiRefGs_vP, *ChiRefGc_vP,
	             *ChiGs_vGc, *ChiGs_vP, *ChiGc_vP,
	             **I_vGs_vGc, **I_vGs_vP, **I_vGc_vP;

	struct S_ELEMENT *ELEMENT;

	// Function pointers
	basis_tdef      basis;
	grad_basis_tdef grad_basis;
	cubature_tdef   cubature;


	// LINE (Includes TP Class)
	EType = LINE;

	setup_ELEMENT_plotting(EType);
	setup_ELEMENT_normals(EType);
	setup_ELEMENT_operators(EType);

	// QUAD
	EType = QUAD;

	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_plotting(EType);
	}

	// HEX
	EType = HEX;

	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_plotting(EType);
	}

	// TRI
	EType = TRI;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);

// Testing
		ELEMENT = get_ELEMENT_type(EType);
array_print_d(6,d,ELEMENT->nr,'R');
	}

	// TET
	EType = TET;
	if (is_ELEMENT_present(EType)) {

	}



/*
	EType = LINE;

	select_functions(&basis,&grad_basis,&cubature,EType);

	ELEMENT = get_ELEMENT_type(LINE);

	dE = ELEMENT->d;

	// Plotting
	plotting_element_info(&rst_vP,&connectivity,&types,&NvnP,&NE,PP,LINE); // free
	free(connectivity);
	free(types);

	// Operators
	NvnGs = ELEMENT->NvnGs;
	NvnGc = ELEMENT->NvnGc;

	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGc_vP  = ELEMENT->I_vGc_vP;


	// Temporary arrays
	ChiGs_vGs       = NULL;
	ChiGc_vGc       = NULL;
*/

	// VOLUME Nodes
/*
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

	// FACET Nodes
	rst_fGc = ELEMENT->rst_fGc; NfnGc = ELEMENT->NfnGc;
	rst_fIs = ELEMENT->rst_fIs; NfnIs = ELEMENT->NfnIs; wfIs = ELEMENT->wfIs;
	rst_fIc = ELEMENT->rst_fIc; NfnIc = ELEMENT->NfnIc; wfIc = ELEMENT->wfIc;
*/

/*
	// VOLUME Nodes (Order Independent)
	cubature(&rst_vGs,&dummyw,&NvnGs[0],0,PGs,dE,"GLL"); // free

	// Preliminary Operators
	IGs = identity_d(NvnGs[0]); // free

	ChiRefGs_vGs = basis(PGs,rst_vGs,NvnGs[0],&Nbf,dE); // free
	

	if (strstr(BasisType,"Modal") != NULL) {
		ChiGs_vGs = ChiRefGs_vGs;
	} else if (strstr(BasisType,"Nodal") != NULL) {
		ChiGs_vGs = IGs;
	}

	ChiRefInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiRefGs_vGs,IGs); // free

	ChiInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiGs_vGs,IGs); // free

	TGs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
	                 NvnGs[0],NvnGs[0],NvnGs[0],1.0,ChiRefInvGs_vGs,ChiGs_vGs); // free

	ChiRefGs_vP = basis(PGs,rst_vP,NvnP,&Nbf,dE); // free

	ChiGs_vP = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,NvnGs[0],NvnGs[0],1.0,ChiRefGs_vP,TGs); // free

	I_vGs_vP[0] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,NvnGs[0],NvnGs[0],1.0,ChiGs_vP,ChiInvGs_vGs); // keep

	free(rst_vGs);

	free(IGs);
	free(ChiRefGs_vGs);
	free(ChiRefInvGs_vGs);
	free(ChiRefGs_vP);
	free(ChiGs_vP);

    for (P = 0; P <= PMax; P++) {

		// VOLUME Nodes (Order dependent)
		cubature(&rst_vGc,&dummyw,&NvnGc[P],0,PGc[P]   ,dE,"GLL"); // free
*/
/*
		cubature(&rst_vCs[P],&dummyw,&NvnCs[P],0,PCs[P][0],dE,"GLL"); // tbd
		cubature(&rst_vCc[P],&dummyw,&NvnCc[P],0,PCc[P][0],dE,"GLL"); // tbd
		cubature(&rst_vJs[P],&dummyw,&NvnCs[P],0,PJs[P][0],dE,"GLL"); // tbd
		cubature(&rst_vJc[P],&dummyw,&NvnCc[P],0,PJc[P][0],dE,"GLL"); // tbd

		cubature(&rst_vS[P]  ,&dummyw,&NvnS[P]  ,0,P         ,dE,NodeTypeS[P][0]);   // tbd
		cubature(&rst_vF[P]  ,&dummyw,&NvnF[P]  ,0,PF[P]     ,dE,NodeTypeF[P][0]);   // tbd
		cubature(&rst_vFrs[P],&dummyw,&NvnFrs[P],0,PFrs[P][0],dE,NodeTypeFrs[P][0]); // tbd
		cubature(&rst_vFrc[P],&dummyw,&NvnFrc[P],0,PFrc[P][0],dE,NodeTypeFrc[P][0]); // tbd

		cubature(&rst_vIs[P],&wvIs[P],&NvnIs[P],1,PIvs[P][0],dE,NodeTypeIvs[P][0]); // tbd
		cubature(&rst_vIc[P],&wvIc[P],&NvnIc[P],1,PIvc[P][0],dE,NodeTypeIvc[P][0]); // tbd

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
*/

/*
		// Preliminary Operators
		IGc = identity_d(NvnGc[P]); // free

		ChiRefGc_vGc = basis(PGc[P],rst_vGc,NvnGc[P],&Nbf,dE); // free

		if (strstr(BasisType,"Modal") != NULL) {
			ChiGc_vGc = ChiRefGc_vGc;
		} else if (strstr(BasisType,"Nodal") != NULL) {
			ChiGc_vGc = IGc;
		}

		ChiRefInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiRefGc_vGc,IGc); // free

		ChiInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiGc_vGc,IGc); // free

		TGc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                 NvnGc[P],NvnGc[P],NvnGc[P],1.0,ChiRefInvGc_vGc,ChiGc_vGc); // free


		ChiRefGs_vGc = basis(PGs,rst_vGc,NvnGc[P],&Nbf,dE); // free
		ChiRefGc_vP  = basis(PGc[P],rst_vP,NvnP,&Nbf,dE);   // free

		ChiGs_vGc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                       NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vGc,TGs); // free
		ChiGc_vP  = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                       NvnP,NvnGc[P],NvnGc[P],1.0,ChiRefGc_vP,TGc);      // free



		// Operators
		I_vGs_vGc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                          NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vGc,ChiInvGs_vGs); // keep

		I_vGc_vP[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		                         NvnP,NvnGc[P],NvnGc[P],1.0,ChiGc_vP,ChiInvGc_vGc);    // keep

		free(rst_vGc);

		free(ChiRefGs_vGc);
		free(ChiGs_vGc);

		free(IGc);
		free(ChiRefGc_vGc);
		free(ChiRefInvGc_vGc);
		free(ChiInvGc_vGc);
		free(TGc);
		free(ChiRefGc_vP);
		free(ChiGc_vP);
	}

	free(rst_vP);

	free(ChiInvGs_vGs);
	free(TGs);
*/
}
