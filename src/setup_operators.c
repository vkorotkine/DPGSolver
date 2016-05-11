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

typedef void (*cubature_tdef) (double **rst, double **w_vec, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                               const unsigned int return_w, const unsigned int P, const unsigned int d,
                               const char *NodeType);
typedef double *(*basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                               const unsigned int d);
typedef double **(*grad_basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn,
                                     unsigned int *NbfOut, const unsigned int d);

static void select_functions(basis_tdef *basis, grad_basis_tdef *grad_basis, cubature_tdef *cubature,
                             const unsigned int type)
{
	if (type == LINE) {
		*basis      = basis_TP;
		*grad_basis = grad_basis_TP;
		*cubature   = cubature_TP;
	} else if (type == TRI) {
		*basis      = basis_SI;
		*grad_basis = grad_basis_SI;
		*cubature   = cubature_TRI;
	} else if (type == TET) {
		*basis      = basis_SI;
		*grad_basis = grad_basis_SI;
		*cubature   = cubature_TET;
	} else if (type == WEDGE) {
		printf("Error: WEDGE elements use a combination of TRI and LINE basis functions/nodes.\n"), exit(1);
	} else if (type == PYR) {
		*basis      = basis_PYR;
		*grad_basis = grad_basis_PYR;
		*cubature   = cubature_PYR;
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
	} else if (EType == WEDGE) {
		Theta_eta[0] = 0.0;        Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = 0.0;        Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = 0.0;        Theta_zeta[2] = 9.0/6.0*PI;
		Theta_eta[3] = 1.0/2.0*PI; Theta_zeta[3] = 0.0;
		Theta_eta[4] = 3.0/2.0*PI; Theta_zeta[4] = 0.0;
	} else if (EType == PYR) {
		double Theta_e = atan(sqrt(2.0));

		Theta_eta[0] = Theta_e - PI/2.0; Theta_zeta[0] = 2.0/2.0*PI;
		Theta_eta[1] = Theta_e - PI/2.0; Theta_zeta[1] = 0.0/2.0*PI;
		Theta_eta[2] = Theta_e - PI/2.0; Theta_zeta[2] = 3.0/2.0*PI;
		Theta_eta[3] = Theta_e - PI/2.0; Theta_zeta[3] = 1.0/2.0*PI;
		Theta_eta[4] = PI/2.0          ; Theta_zeta[4] = 0.0       ;
	} else {
		printf("Add support setup_ELEMENT_normals\n"), exit(1);
	}

	for (f = 0; f < Nf; f++) {
		           nr[f*d+0] =  cos(Theta_eta[f])*cos(Theta_zeta[f]);
		if (d > 1) nr[f*d+1] =  cos(Theta_eta[f])*sin(Theta_zeta[f]);
		if (d > 2) nr[f*d+2] = -sin(Theta_eta[f]);
	}
}

static void setup_ELEMENT_operators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc;
	double       **ICs, **ICc,
	             **I_vGs_vP, **I_vGs_vGc, **I_vGs_vCs, **I_vGs_vJs,
	             **I_vGc_vP,              **I_vGc_vCc, **I_vGc_vJc,
	             ***D_vGs_vCs, ***D_vGs_vJs,
	             ***D_vGc_vCc, ***D_vGc_vJc,
				 ***D_vCs_vCs,
				 ***D_vCc_vCc;

	// Initialize DB Parameters
	unsigned int PMax        = DB.PMax,
	             PGs         = DB.PGs,
	             *PGc        = DB.PGc,
	             **PCs        = DB.PCs,
	             **PCc        = DB.PCc,
	             **PJs        = DB.PJs,
	             **PJc        = DB.PJc,
	             PP          = DB.PP;

	char         *BasisType  = DB.BasisType,
	             **NodeTypeG = DB.NodeTypeG;

	// Standard datatypes
	unsigned int dim, dE, P,
	             Nbf, Eclass,
	             NvnP,
	             dummy_ui, *dummyPtr_ui[2];
	double       *rst_vP, *rst_vGs, *rst_vGc, *rst_vCs, *rst_vCc, *rst_vJs, *rst_vJc,
	             *IGs, *IGc,
	             *TGs, *TGc, *TCs, *TCc,
	             *ChiRefGs_vGs, *ChiRefGc_vGc, *ChiRefCs_vCs, *ChiRefCc_vCc,
	             *ChiGs_vGs,    *ChiGc_vGc,    *ChiCs_vCs,    *ChiCc_vCc,
	             *ChiRefInvGs_vGs, *ChiRefInvGc_vGc, *ChiRefInvCs_vCs, *ChiRefInvCc_vCc,
	             *ChiInvGs_vGs,    *ChiInvGc_vGc,    *ChiInvCs_vCs,    *ChiInvCc_vCc,
	             *ChiRefGs_vP, *ChiRefGs_vGc, *ChiRefGs_vCs, *ChiRefGs_vJs,
	             *ChiRefGc_vP,                *ChiRefGc_vCc, *ChiRefGc_vJc,
	             *ChiGs_vP, *ChiGs_vGc, *ChiGs_vCs, *ChiGs_vJs,
	             *ChiGc_vP,             *ChiGc_vCc, *ChiGc_vJc,
				 **GradChiRefGs_vCs, **GradChiRefGs_vJs,
				 **GradChiRefGc_vCc, **GradChiRefGc_vJc,
				 **GradChiRefCs_vCs, **GradChiRefCc_vCc,
				 **GradChiGs_vCs, **GradChiGs_vJs,
				 **GradChiGc_vCc, **GradChiGc_vJc,
				 **GradChiCs_vCs, **GradChiCc_vCc,
	             *dummyPtr_d;

	struct S_ELEMENT *ELEMENT;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	// silence
	ChiGs_vGs = NULL; ChiGc_vGc = NULL;
	ChiCs_vCs = NULL; ChiCc_vCc = NULL;

	ELEMENT = get_ELEMENT_type(EType);

	// No need to consider the second Eclass as WEDGE basis functions will not be built
	Eclass = get_Eclass(EType);

	dE = ELEMENT->d;

	select_functions(&basis,&grad_basis,&cubature,EType);

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnGc = ELEMENT->NvnGc;
	NvnCs = ELEMENT->NvnCs;
	NvnCc = ELEMENT->NvnCc;
	NvnJs = ELEMENT->NvnJs;
	NvnJc = ELEMENT->NvnJc;

	ICs = ELEMENT->ICs;
	ICc = ELEMENT->ICc;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vCs = ELEMENT->I_vGs_vCs;
	I_vGs_vJs = ELEMENT->I_vGs_vJs;
	I_vGc_vP  = ELEMENT->I_vGc_vP;
	I_vGc_vCc = ELEMENT->I_vGc_vCc;
	I_vGc_vJc = ELEMENT->I_vGc_vJc;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vJs = ELEMENT->D_vGs_vJs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vJc = ELEMENT->D_vGc_vJc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	// Allocate memory for arrays with multiple levels of dereferencing
	GradChiGs_vCs = malloc(dE * sizeof *GradChiGs_vCs); // free
	GradChiGs_vJs = malloc(dE * sizeof *GradChiGs_vJs); // free
	GradChiGc_vCc = malloc(dE * sizeof *GradChiGc_vCc); // free
	GradChiGc_vJc = malloc(dE * sizeof *GradChiGc_vJc); // free
	GradChiCs_vCs = malloc(dE * sizeof *GradChiGs_vCs); // free
	GradChiCc_vCc = malloc(dE * sizeof *GradChiGc_vCc); // free

	// VOLUME Nodes (Order Independent)
	plotting_element_info(&rst_vP,&dummyPtr_ui[0],&dummyPtr_ui[1],&NvnP,&dummy_ui,PP,EType); // free
	free(dummyPtr_ui[0]);
	free(dummyPtr_ui[1]);

	cubature(&rst_vGs,&dummyPtr_d,&dummyPtr_ui[0],&NvnGs[0],&dummy_ui,0,PGs,dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free

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
		cubature(&rst_vGc,&dummyPtr_d,&dummyPtr_ui[0],&NvnGc[P],&dummy_ui,0,PGc[P]        ,dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
		cubature(&rst_vCs,&dummyPtr_d,&dummyPtr_ui[0],&NvnCs[P],&dummy_ui,0,PCs[P][Eclass],dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
		cubature(&rst_vCc,&dummyPtr_d,&dummyPtr_ui[0],&NvnCc[P],&dummy_ui,0,PCc[P][Eclass],dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
		cubature(&rst_vJs,&dummyPtr_d,&dummyPtr_ui[0],&NvnJs[P],&dummy_ui,0,PJs[P][Eclass],dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
		cubature(&rst_vJc,&dummyPtr_d,&dummyPtr_ui[0],&NvnJc[P],&dummy_ui,0,PJc[P][Eclass],dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free

		// Preliminary Operators
		IGc    = identity_d(NvnGc[P]); // free
		ICs[P] = identity_d(NvnCs[P]); // keep
		ICc[P] = identity_d(NvnCc[P]); // keep

		ChiRefGc_vGc = basis(PGc[P],rst_vGc,NvnGc[P],&Nbf,dE); // free
		ChiRefCs_vCs = basis(PCs[P][Eclass],rst_vCs,NvnCs[P],&Nbf,dE); // free
		ChiRefCc_vCc = basis(PCc[P][Eclass],rst_vCc,NvnCc[P],&Nbf,dE); // free

		if (strstr(BasisType,"Modal") != NULL) {
			ChiGc_vGc = ChiRefGc_vGc;
			ChiCs_vCs = ChiRefCs_vCs;
			ChiCc_vCc = ChiRefCc_vCc;
		} else if (strstr(BasisType,"Nodal") != NULL) {
			ChiGc_vGc = IGc;
			ChiCs_vCs = ICs[P];
			ChiCc_vCc = ICc[P];
		}

		ChiRefInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiRefGc_vGc,IGc);    // free
		ChiRefInvCs_vCs = inverse_d(NvnCs[P],NvnCs[P],ChiRefCs_vCs,ICs[P]); // free
		ChiRefInvCc_vCc = inverse_d(NvnCc[P],NvnCc[P],ChiRefCc_vCc,ICc[P]); // free

		ChiInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiGc_vGc,IGc);    // free
		ChiInvCs_vCs = inverse_d(NvnCs[P],NvnCs[P],ChiCs_vCs,ICs[P]); // free
		ChiInvCc_vCc = inverse_d(NvnCc[P],NvnCc[P],ChiCc_vCc,ICc[P]); // free

		TGc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnGc[P],NvnGc[P],NvnGc[P],1.0,ChiRefInvGc_vGc,ChiGc_vGc); // free
		TCs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnCs[P],NvnCs[P],1.0,ChiRefInvCs_vCs,ChiCs_vCs); // free
		TCc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnCc[P],NvnCc[P],1.0,ChiRefInvCc_vCc,ChiCc_vCc); // free

		ChiRefGs_vGc = basis(PGs,   rst_vGc,NvnGc[P],&Nbf,dE); // free
		ChiRefGs_vCs = basis(PGs,   rst_vCs,NvnCs[P],&Nbf,dE); // free
		ChiRefGs_vJs = basis(PGs,   rst_vJs,NvnJs[P],&Nbf,dE); // free
		ChiRefGc_vP  = basis(PGc[P],rst_vP, NvnP,    &Nbf,dE); // free
		ChiRefGc_vCc = basis(PGc[P],rst_vCc,NvnCc[P],&Nbf,dE); // free
		ChiRefGc_vJc = basis(PGc[P],rst_vJc,NvnJc[P],&Nbf,dE); // free

		ChiGs_vGc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vGc,TGs); // free
		ChiGs_vCs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vCs,TGs); // free
		ChiGs_vJs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJs[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vJs,TGs); // free
		ChiGc_vP  = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiRefGc_vP, TGc); // free
		ChiGc_vCc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vCc,TGc); // free
		ChiGc_vJc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vJc,TGc); // free

		GradChiRefGs_vCs = grad_basis(PGs,           rst_vCs,NvnCs[P],&Nbf,dE); // free
		GradChiRefGs_vJs = grad_basis(PGs,           rst_vJs,NvnJs[P],&Nbf,dE); // free
		GradChiRefGc_vCc = grad_basis(PGc[P],        rst_vCc,NvnCc[P],&Nbf,dE); // free
		GradChiRefGc_vJc = grad_basis(PGc[P],        rst_vJc,NvnJc[P],&Nbf,dE); // free
		GradChiRefCs_vCs = grad_basis(PCs[P][Eclass],rst_vCs,NvnCs[P],&Nbf,dE); // free
		GradChiRefCc_vCc = grad_basis(PCc[P][Eclass],rst_vCc,NvnCc[P],&Nbf,dE); // free

		for (dim = 0; dim < dE; dim++) {
			GradChiGs_vCs[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,GradChiRefGs_vCs[dim],TGs); // free
			GradChiGs_vJs[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJs[P],NvnGs[0],NvnGs[0],1.0,GradChiRefGs_vJs[dim],TGs); // free
			GradChiGc_vCc[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vCc[dim],TGc); // free
			GradChiGc_vJc[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vJc[dim],TGc); // free
			GradChiCs_vCs[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiRefCs_vCs[dim],TCs); // free
			GradChiCc_vCc[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiRefCc_vCc[dim],TCc); // free
		}

		// Returned Operators
		I_vGs_vGc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vGc,ChiInvGs_vGs); // keep
		I_vGs_vCs[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vCs,ChiInvGs_vGs); // keep
		I_vGs_vJs[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJs[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vJs,ChiInvGs_vGs); // keep
		I_vGc_vP[P]  = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiGc_vP, ChiInvGc_vGc); // keep
		I_vGc_vCc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vCc,ChiInvGc_vGc); // keep
		I_vGc_vJc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vJc,ChiInvGc_vGc); // keep

		for (dim = 0; dim < dE; dim++) {
			D_vGs_vCs[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,GradChiGs_vCs[dim],ChiInvGs_vGs); // keep
			D_vGs_vJs[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJs[P],NvnGs[0],NvnGs[0],1.0,GradChiGs_vJs[dim],ChiInvGs_vGs); // keep
			D_vGc_vCc[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vCc[dim],ChiInvGc_vGc); // keep
			D_vGc_vJc[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnJc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vJc[dim],ChiInvGc_vGc); // keep
			D_vCs_vCs[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiCs_vCs[dim],ChiInvCs_vCs); // keep
			D_vCc_vCc[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiCc_vCc[dim],ChiInvCc_vCc); // keep
		}

/*
if (P == 2 && EType == TRI) {
	array_print_d(NvnCc[P],NvnGc[P],GradChiRefGc_vCc[0],'R');
	array_print_d(NvnCc[P],NvnGc[P],GradChiRefGc_vCc[1],'R');
	exit(1);
}
*/

		free(rst_vGc);
		free(rst_vCs);
		free(rst_vCc);
		free(rst_vJs);
		free(rst_vJc);

		free(IGc);
		free(ChiRefGc_vGc);
		free(ChiRefCs_vCs);
		free(ChiRefCc_vCc);

		free(ChiRefInvGc_vGc);
		free(ChiRefInvCs_vCs);
		free(ChiRefInvCc_vCc);

		free(ChiInvGc_vGc);
		free(ChiInvCs_vCs);
		free(ChiInvCc_vCc);

		free(TGc);
		free(TCs);
		free(TCc);

		free(ChiRefGs_vGc);
		free(ChiRefGs_vCs);
		free(ChiRefGs_vJs);
		free(ChiRefGc_vP);
		free(ChiRefGc_vCc);
		free(ChiRefGc_vJc);

		free(ChiGs_vGc);
		free(ChiGs_vCs);
		free(ChiGs_vJs);
		free(ChiGc_vP);
		free(ChiGc_vCc);
		free(ChiGc_vJc);

		array_free2_d(dE,GradChiRefGs_vCs);
		array_free2_d(dE,GradChiRefGs_vJs);
		array_free2_d(dE,GradChiRefGc_vCc);
		array_free2_d(dE,GradChiRefGc_vJc);
		array_free2_d(dE,GradChiRefCs_vCs);
		array_free2_d(dE,GradChiRefCc_vCc);

		for (dim = 0; dim < dE; dim++) {
			free(GradChiGs_vCs[dim]);
			free(GradChiGs_vJs[dim]);
			free(GradChiGc_vCc[dim]);
			free(GradChiGc_vJc[dim]);
			free(GradChiCs_vCs[dim]);
			free(GradChiCc_vCc[dim]);
		}
	}

	free(rst_vP);

	free(ChiInvGs_vGs);
	free(TGs);

	free(GradChiGs_vCs);
	free(GradChiGs_vJs);
	free(GradChiGc_vCc);
	free(GradChiGc_vJc);
	free(GradChiCs_vCs);
	free(GradChiCc_vCc);
}

void setup_operators(void)
{
	// Initialize DB Parameters

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int EType;



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
		setup_ELEMENT_operators(EType);

// Testing
//		ELEMENT = get_ELEMENT_type(EType);
//array_print_d(6,d,ELEMENT->nr,'R');
	}

	// TET
	EType = TET;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);

	}

	// PYR
	EType = PYR;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);

	}

	// WEDGE
	EType = WEDGE;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_plotting(EType);

	}
}
