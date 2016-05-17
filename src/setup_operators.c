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

struct S_BCOORDS {
	unsigned int Nve,
 	             *NfnIs, *NfnIc;
	double       **BCoords_Is, **BCoords_Ic;
};

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
	switch(type) {
	case LINE:
	case QUAD:
	case HEX:
		*basis      = basis_TP;
		*grad_basis = grad_basis_TP;
		*cubature   = cubature_TP;
		break;
	case TRI:
		*basis      = basis_SI;
		*grad_basis = grad_basis_SI;
		*cubature   = cubature_TRI;
		break;
	case TET:
		*basis      = basis_SI;
		*grad_basis = grad_basis_SI;
		*cubature   = cubature_TET;
		break;
	case WEDGE:
		printf("Error: WEDGE elements use a combination of TRI and LINE basis functions/nodes.\n"), exit(1);
		break;
	case PYR:
		*basis      = basis_PYR;
		*grad_basis = grad_basis_PYR;
		*cubature   = cubature_PYR;
		break;
	default:
		printf("%d\n",type);
		printf("Error: Unsupported type in select_functions.\n"), exit(1);
		break;
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
		Theta_eta[0] = 0.0; Theta_zeta[0] = PI;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 0.0;
	} else if (EType == QUAD) {
		Theta_eta[0] = 0.0; Theta_zeta[0] = PI;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 0.0;
		Theta_eta[2] = 0.0; Theta_zeta[2] = 1.5*PI;
		Theta_eta[3] = 0.0; Theta_zeta[3] = 0.5*PI;
	} else if (EType == HEX) {
		Theta_eta[0] = 0.0;    Theta_zeta[0] = PI;
		Theta_eta[1] = 0.0;    Theta_zeta[1] = 0.0;
		Theta_eta[2] = 0.0;    Theta_zeta[2] = 1.5*PI;
		Theta_eta[3] = 0.0;    Theta_zeta[3] = 0.5*PI;
		Theta_eta[4] = 0.5*PI; Theta_zeta[4] = 0.0;
		Theta_eta[5] = 1.5*PI; Theta_zeta[5] = PI;
	} else if (EType == TRI) {
		Theta_eta[0] = 0.0; Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = 0.0; Theta_zeta[2] = 9.0/6.0*PI;
	} else if (EType == TET) {
		double Theta_e = atan(2.0*sqrt(2.0));

		Theta_eta[0] = Theta_e - 0.5*PI; Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = Theta_e - 0.5*PI; Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = Theta_e - 0.5*PI; Theta_zeta[2] = 9.0/6.0*PI;
		Theta_eta[3] = 0.5*PI;           Theta_zeta[3] = 0.0;
	} else if (EType == WEDGE) {
		Theta_eta[0] = 0.0;    Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = 0.0;    Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = 0.0;    Theta_zeta[2] = 9.0/6.0*PI;
		Theta_eta[3] = 0.5*PI; Theta_zeta[3] = 0.0;
		Theta_eta[4] = 1.5*PI; Theta_zeta[4] = 0.0;
	} else if (EType == PYR) {
		double Theta_e = atan(sqrt(2.0));

		Theta_eta[0] = Theta_e - 0.5*PI; Theta_zeta[0] = PI;
		Theta_eta[1] = Theta_e - 0.5*PI; Theta_zeta[1] = 0.0;
		Theta_eta[2] = Theta_e - 0.5*PI; Theta_zeta[2] = 1.5*PI;
		Theta_eta[3] = Theta_e - 0.5*PI; Theta_zeta[3] = 0.5*PI;
		Theta_eta[4] = 0.5*PI;           Theta_zeta[4] = 0.0;
	} else {
		printf("Add support setup_ELEMENT_normals\n"), exit(1);
	}

	for (f = 0; f < Nf; f++) {
		           nr[f*d+0] =  cos(Theta_eta[f])*cos(Theta_zeta[f]);
		if (d > 1) nr[f*d+1] =  cos(Theta_eta[f])*sin(Theta_zeta[f]);
		if (d > 2) nr[f*d+2] = -sin(Theta_eta[f]);
	}
}

static double *get_rst_vC(const struct S_ELEMENT *ELEMENT)
{
	unsigned int d, Nve, EType;
	double *rst_vC;

	d     = ELEMENT->d;
	Nve   = ELEMENT->Nve;
	EType = ELEMENT->type;

	rst_vC = malloc(Nve*d * sizeof *rst_vC); // keep (requires external free)

	switch (EType) {
	case LINE:
		rst_vC[0*Nve+0] = -1.0;
		rst_vC[0*Nve+1] =  1.0;
		break;
	case TRI:
		rst_vC[0*Nve+0] = -1.0; rst_vC[1*Nve+0] = -1.0/sqrt(3.0);
		rst_vC[0*Nve+1] =  1.0; rst_vC[1*Nve+1] = -1.0/sqrt(3.0);
		rst_vC[0*Nve+2] =  0.0; rst_vC[1*Nve+2] =  2.0/sqrt(3.0);
		break;
	case QUAD:
		rst_vC[0*Nve+0] = -1.0; rst_vC[1*Nve+0] = -1.0;
		rst_vC[0*Nve+1] =  1.0; rst_vC[1*Nve+1] = -1.0;
		rst_vC[0*Nve+2] =  1.0; rst_vC[1*Nve+2] =  1.0;
		rst_vC[0*Nve+3] = -1.0; rst_vC[1*Nve+3] =  1.0;
		break;
	case TET:
		rst_vC[0*Nve+0] = -1.0; rst_vC[1*Nve+0] = -1.0/sqrt(3.0); rst_vC[2*Nve+0] = -1.0/sqrt(6.0);
		rst_vC[0*Nve+1] =  1.0; rst_vC[1*Nve+1] = -1.0/sqrt(3.0); rst_vC[2*Nve+1] = -1.0/sqrt(6.0);
		rst_vC[0*Nve+2] =  0.0; rst_vC[1*Nve+2] =  2.0/sqrt(3.0); rst_vC[2*Nve+2] = -1.0/sqrt(6.0);
		rst_vC[0*Nve+3] =  0.0; rst_vC[1*Nve+3] =  0.0;           rst_vC[2*Nve+3] =  3.0/sqrt(6.0);
		break;
	case HEX:
		rst_vC[0*Nve+0] = -1.0; rst_vC[1*Nve+0] = -1.0; rst_vC[2*Nve+0] = -1.0;
		rst_vC[0*Nve+1] =  1.0; rst_vC[1*Nve+1] = -1.0; rst_vC[2*Nve+1] = -1.0;
		rst_vC[0*Nve+2] =  1.0; rst_vC[1*Nve+2] =  1.0; rst_vC[2*Nve+2] = -1.0;
		rst_vC[0*Nve+3] = -1.0; rst_vC[1*Nve+3] =  1.0; rst_vC[2*Nve+3] = -1.0;
		rst_vC[0*Nve+4] = -1.0; rst_vC[1*Nve+4] = -1.0; rst_vC[2*Nve+4] =  1.0;
		rst_vC[0*Nve+5] =  1.0; rst_vC[1*Nve+5] = -1.0; rst_vC[2*Nve+5] =  1.0;
		rst_vC[0*Nve+6] =  1.0; rst_vC[1*Nve+6] =  1.0; rst_vC[2*Nve+6] =  1.0;
		rst_vC[0*Nve+7] = -1.0; rst_vC[1*Nve+7] =  1.0; rst_vC[2*Nve+7] =  1.0;
		break;
	case WEDGE:
		rst_vC[0*Nve+0] = -1.0; rst_vC[1*Nve+0] = -1.0/sqrt(3.0); rst_vC[2*Nve+0] = -1.0;
		rst_vC[0*Nve+1] =  1.0; rst_vC[1*Nve+1] = -1.0/sqrt(3.0); rst_vC[2*Nve+1] = -1.0;
		rst_vC[0*Nve+2] =  0.0; rst_vC[1*Nve+2] =  2.0/sqrt(3.0); rst_vC[2*Nve+2] = -1.0;
		rst_vC[0*Nve+3] = -1.0; rst_vC[1*Nve+3] = -1.0/sqrt(3.0); rst_vC[2*Nve+3] =  1.0;
		rst_vC[0*Nve+4] =  1.0; rst_vC[1*Nve+4] = -1.0/sqrt(3.0); rst_vC[2*Nve+4] =  1.0;
		rst_vC[0*Nve+5] =  0.0; rst_vC[1*Nve+5] =  2.0/sqrt(3.0); rst_vC[2*Nve+5] =  1.0;
		break;
	case PYR:
		rst_vC[0*Nve+0] = -1.0; rst_vC[1*Nve+0] = -1.0; rst_vC[2*Nve+0] = -1.0/5.0*sqrt(2.0);
		rst_vC[0*Nve+1] =  1.0; rst_vC[1*Nve+1] = -1.0; rst_vC[2*Nve+1] = -1.0/5.0*sqrt(2.0);
		rst_vC[0*Nve+2] =  1.0; rst_vC[1*Nve+2] =  1.0; rst_vC[2*Nve+2] = -1.0/5.0*sqrt(2.0);
		rst_vC[0*Nve+3] = -1.0; rst_vC[1*Nve+3] =  1.0; rst_vC[2*Nve+3] = -1.0/5.0*sqrt(2.0);
		rst_vC[0*Nve+4] =  0.0; rst_vC[1*Nve+4] =  0.0; rst_vC[2*Nve+4] =  4.0/5.0*sqrt(2.0);
		break;
	default:
		printf("Error: Unsupported EType in get_rst_vC.\n"), exit(1);
		break;
	}
	return rst_vC;
}

static struct S_BCOORDS *get_BCoords_dEm1(const struct S_ELEMENT *ELEMENT, const unsigned int IndFType)
{
	// Initialize DB Parameters
	unsigned int NP   = DB.NP,
	             PMax = DB.PMax;

	// Standard datatypes
	unsigned int P, EType, Nve,
	             *NfnIs, *NfnIc;
	double       **BCoords_Is, **BCoords_Ic;

	struct S_BCOORDS *BCoords_dEm1;
	struct S_ELEMENT *ELEMENT_F;

	BCoords_dEm1 = malloc(sizeof *BCoords_dEm1);      // keep (requires external free)
	NfnIs        = malloc(NP * sizeof *(NfnIs));      // keep (requires external free)
	NfnIc        = malloc(NP * sizeof *(NfnIc));      // keep (requires external free)
	BCoords_Is   = malloc(NP * sizeof *(BCoords_Is)); // keep (requires external free)
	BCoords_Ic   = malloc(NP * sizeof *(BCoords_Ic)); // keep (requires external free)

	EType = ELEMENT->type;
	if (EType == LINE) {
// fix this to be consistent with treatment of other elements if possible (ToBeDeleted)
		ELEMENT_F = get_ELEMENT_type(POINT);

		Nve = 1;
		for (P = 0; P <= PMax; P++) {
			NfnIs[P] = 1;
			NfnIc[P] = 1;

			BCoords_Is[P] = malloc(1 * sizeof **(BCoords_Is));
			BCoords_Ic[P] = malloc(1 * sizeof **(BCoords_Ic));

			BCoords_Is[P][0] = 1.0;
			BCoords_Ic[P][0] = 1.0;
		}
	} else {
		if (EType == TRI || EType == QUAD) {
			ELEMENT_F = get_ELEMENT_type(LINE);
		} else if (EType == TET || (EType == WEDGE && IndFType == 1) || (EType == PYR && IndFType == 0)) {
			ELEMENT_F = get_ELEMENT_type(TRI);
		} else if (EType == HEX || (EType == WEDGE && IndFType == 0) || (EType == PYR && IndFType == 1)) {
			ELEMENT_F = get_ELEMENT_type(QUAD);
		} else {
			printf("Error: Unsupported EType/IndFType combination in get_BCoords_dEm1.\n"), exit(1);
		}

		// Initialize DB Parameters
		unsigned int **PIfs         = DB.PIfs,
		             **PIfc         = DB.PIfc;
		char         **NodeTypeG    = DB.NodeTypeG,
		             ***NodeTypeIfs = DB.NodeTypeIfs,
		             ***NodeTypeIfc = DB.NodeTypeIfc;

		// Standard datatypes
		unsigned int dE, Nbf,
		             NfnGs, EType, Eclass,
		             dummy_ui, *dummyPtr_ui;
		double       *rst_vGs, *rst_fIs, *rst_fIc,
		             *dummyPtr_d[2],
		             *IGs,
		             *ChiRefGs_vGs,
		             *ChiRefInvGs_vGs,
		             *ChiRefGs_fIs, *ChiRefGs_fIc;

		// Function pointers
		cubature_tdef   cubature;
		basis_tdef      basis;
		grad_basis_tdef grad_basis;

		EType = ELEMENT_F->type;
		select_functions(&basis,&grad_basis,&cubature,EType);

		Eclass = get_Eclass(EType);


		dE = ELEMENT_F->d;

		Nve = ELEMENT_F->Nve;

		// It is important to use the nodes corresponding to the VeF ordering
		rst_vGs = get_rst_vC(ELEMENT_F); // free
		cubature(&dummyPtr_d[0],&dummyPtr_d[1],&dummyPtr_ui,&NfnGs,&dummy_ui,0,1,dE,NodeTypeG[Eclass]); // free
		free(dummyPtr_ui);
		free(dummyPtr_d[0]);

		IGs             = identity_d(NfnGs);                       // free
		ChiRefGs_vGs    = basis(1,rst_vGs,NfnGs,&Nbf,dE);          // free
		ChiRefInvGs_vGs = inverse_d(NfnGs,NfnGs,ChiRefGs_vGs,IGs); // free

		free(rst_vGs);

		free(IGs);
		free(ChiRefGs_vGs);

		for (P = 0; P <= PMax; P++) {
			cubature(&rst_fIs,&dummyPtr_d[0],&dummyPtr_ui,&NfnIs[P],&dummy_ui,0,PIfs[P][Eclass],dE,NodeTypeIfs[P][Eclass]); free(dummyPtr_ui); // free
			cubature(&rst_fIc,&dummyPtr_d[0],&dummyPtr_ui,&NfnIc[P],&dummy_ui,0,PIfc[P][Eclass],dE,NodeTypeIfc[P][Eclass]); free(dummyPtr_ui); // free

			ChiRefGs_fIs = basis(1,rst_fIs,NfnIs[P],&Nbf,dE); // free
			ChiRefGs_fIc = basis(1,rst_fIc,NfnIc[P],&Nbf,dE); // free

			BCoords_Is[P] = mm_Alloc_d(CblasColMajor,CblasTrans,CblasTrans,NfnIs[P],NfnGs,NfnGs,1.0,ChiRefGs_fIs,ChiRefInvGs_vGs); // keep
			BCoords_Ic[P] = mm_Alloc_d(CblasColMajor,CblasTrans,CblasTrans,NfnIc[P],NfnGs,NfnGs,1.0,ChiRefGs_fIc,ChiRefInvGs_vGs); // keep

			free(rst_fIs);
			free(rst_fIc);

			free(ChiRefGs_fIs);
			free(ChiRefGs_fIc);
		}
		free(ChiRefInvGs_vGs);
	}

//printf("%d %d\n",ELEMENT->type,ELEMENT_F->type);

	BCoords_dEm1->Nve   = Nve;
	BCoords_dEm1->NfnIs = NfnIs;
	BCoords_dEm1->NfnIc = NfnIc;
	BCoords_dEm1->BCoords_Is = BCoords_Is;
	BCoords_dEm1->BCoords_Ic = BCoords_Ic;

	return BCoords_dEm1;
}

static void setup_ELEMENT_operators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, **NfnIs, **NfnIc;
	double       **ICs, **ICc,
	             **I_vGs_vP, **I_vGs_vGc, **I_vGs_vCs, **I_vGs_vIs, ****I_vGs_fIs, ****I_vGs_fIc,
	             **I_vGc_vP,              **I_vGc_vCc, **I_vGc_vIc, ****I_vGc_fIs, ****I_vGc_fIc,
	             **I_vCs_vIs, **I_vCs_vIc, ****I_vCs_fIs, ****I_vCs_fIc,
	             **I_vCc_vIs, **I_vCc_vIc, ****I_vCc_fIs, ****I_vCc_fIc,
	             ***D_vGs_vCs, ***D_vGs_vIs,
	             ***D_vGc_vCc, ***D_vGc_vIc,
	             ***D_vCs_vCs,
	             ***D_vCc_vCc;

	// Initialize DB Parameters
	unsigned int NfMax       = DB.NfMax,
	             NfveMax     = DB.NfveMax,
	             NveMax      = DB.NveMax,
	             NfrefMax    = DB.NfrefMax,
	             PMax        = DB.PMax,
	             PGs         = DB.PGs,
	             *PGc        = DB.PGc,
	             **PCs       = DB.PCs,
	             **PCc       = DB.PCc,
	             **PIvs      = DB.PIvs,
	             **PIvc      = DB.PIvc,
	             PP          = DB.PP;

	char         *BasisType     = DB.BasisType,
	             **NodeTypeG    = DB.NodeTypeG,
	             ***NodeTypeIvs = DB.NodeTypeIvs,
	             ***NodeTypeIvc = DB.NodeTypeIvc;

	// Standard datatypes
	unsigned int dim, dE, P, f, IndFType, ve, veMax, Pb, PbMin, PbMax,
	             Nve, Nf, Nbf, Eclass, NFTypes,
	             NvnP,
	             B_Nve[2], *Nfve,
	             dummy_ui, *dummyPtr_ui[2];
	double       *E_rst_vC, *rst_vC, *VeF,
	             *rst_vP, *rst_vGs, *rst_vGc, *rst_vCs, *rst_vCc, *rst_vIs, *rst_vIc, **rst_fIs, **rst_fIc,
	             *IGs, *IGc,
	             *TGs, *TGc, *TCs, *TCc,
	             *ChiRefGs_vGs, *ChiRefGc_vGc, *ChiRefCs_vCs, *ChiRefCc_vCc,
	             *ChiGs_vGs,    *ChiGc_vGc,    *ChiCs_vCs,    *ChiCc_vCc,
	             *ChiRefInvGs_vGs, *ChiRefInvGc_vGc, *ChiRefInvCs_vCs, *ChiRefInvCc_vCc,
	             *ChiInvGs_vGs,    *ChiInvGc_vGc,    *ChiInvCs_vCs,    *ChiInvCc_vCc,
	             *ChiRefGs_vP, *ChiRefGs_vGc, *ChiRefGs_vCs, *ChiRefGs_vIs, **ChiRefGs_fIs, **ChiRefGs_fIc,
	             *ChiRefGc_vP,                *ChiRefGc_vCc, *ChiRefGc_vIc, **ChiRefGc_fIs, **ChiRefGc_fIc,
	             *ChiRefCs_vIs, *ChiRefCs_vIc, **ChiRefCs_fIs, **ChiRefCs_fIc,
	             *ChiRefCc_vIs, *ChiRefCc_vIc, **ChiRefCc_fIs, **ChiRefCc_fIc,
	             *ChiGs_vP, *ChiGs_vGc, *ChiGs_vCs, *ChiGs_vIs, **ChiGs_fIs, **ChiGs_fIc,
	             *ChiGc_vP,             *ChiGc_vCc, *ChiGc_vIc, **ChiGc_fIs, **ChiGc_fIc,
	             *ChiCs_vIs, *ChiCs_vIc, **ChiCs_fIs, **ChiCs_fIc,
	             *ChiCc_vIs, *ChiCc_vIc, **ChiCc_fIs, **ChiCc_fIc,
	             **GradChiRefGs_vCs, **GradChiRefGs_vIs,
	             **GradChiRefGc_vCc, **GradChiRefGc_vIc,
	             **GradChiRefCs_vCs, **GradChiRefCc_vCc,
	             **GradChiGs_vCs, **GradChiGs_vIs,
	             **GradChiGc_vCc, **GradChiGc_vIc,
	             **GradChiCs_vCs, **GradChiCc_vCc,
	             *dummyPtr_d;

	struct BCoords {
	       	double **Is, **Ic;
	       } *BCoords[2];

	struct S_BCOORDS *BCoords_dEm1[2];
	struct S_ELEMENT *ELEMENT;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	// silence
	ChiGs_vGs = NULL; ChiGc_vGc = NULL;
	ChiCs_vCs = NULL; ChiCc_vCc = NULL;

	ELEMENT = get_ELEMENT_type(EType);

	// No need to consider the second Eclass as WEDGE basis functions will be built through a combination of lower
	// dimensional operators.
	Eclass = get_Eclass(EType);

	dE   = ELEMENT->d;
	Nve  = ELEMENT->Nve;
	Nf   = ELEMENT->Nf;
	Nfve = ELEMENT->Nfve;
	VeF  = ELEMENT->VeF;

	E_rst_vC = get_rst_vC(ELEMENT);
	rst_vC = malloc(Nve*dE * sizeof *rst_vC); // free

	select_functions(&basis,&grad_basis,&cubature,EType);

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnGc = ELEMENT->NvnGc;
	NvnCs = ELEMENT->NvnCs;
	NvnCc = ELEMENT->NvnCc;
	NvnIs = ELEMENT->NvnIs;
	NvnIc = ELEMENT->NvnIc;
	NfnIs = ELEMENT->NfnIs;
	NfnIc = ELEMENT->NfnIc;

	ICs = ELEMENT->ICs;
	ICc = ELEMENT->ICc;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vCs = ELEMENT->I_vGs_vCs;
	I_vGs_vIs = ELEMENT->I_vGs_vIs;
	I_vGc_vP  = ELEMENT->I_vGc_vP;
	I_vGc_vCc = ELEMENT->I_vGc_vCc;
	I_vGc_vIc = ELEMENT->I_vGc_vIc;
	I_vCs_vIs = ELEMENT->I_vCs_vIs;
	I_vCs_vIc = ELEMENT->I_vCs_vIc;
	I_vCc_vIs = ELEMENT->I_vCc_vIs;
	I_vCc_vIc = ELEMENT->I_vCc_vIc;

	I_vGs_fIs = ELEMENT->I_vGs_fIs;
	I_vGs_fIc = ELEMENT->I_vGs_fIc;
	I_vGc_fIs = ELEMENT->I_vGc_fIs;
	I_vGc_fIc = ELEMENT->I_vGc_fIc;
	I_vCs_fIs = ELEMENT->I_vCs_fIs;
	I_vCs_fIc = ELEMENT->I_vCs_fIc;
	I_vCc_fIs = ELEMENT->I_vCc_fIs;
	I_vCc_fIc = ELEMENT->I_vCc_fIc;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vIs = ELEMENT->D_vGs_vIs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vIc = ELEMENT->D_vGc_vIc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	// Allocate memory for arrays with multiple levels of dereferencing
	rst_fIs       = malloc(NfMax * sizeof *rst_fIs); // free
	rst_fIc       = malloc(NfMax * sizeof *rst_fIc); // free

	GradChiGs_vCs = malloc(dE * sizeof *GradChiGs_vCs); // free
	GradChiGs_vIs = malloc(dE * sizeof *GradChiGs_vIs); // free
	GradChiGc_vCc = malloc(dE * sizeof *GradChiGc_vCc); // free
	GradChiGc_vIc = malloc(dE * sizeof *GradChiGc_vIc); // free
	GradChiCs_vCs = malloc(dE * sizeof *GradChiGs_vCs); // free
	GradChiCc_vCc = malloc(dE * sizeof *GradChiGc_vCc); // free

	ChiRefGs_fIs  = malloc(NfMax * sizeof *ChiRefGs_fIs); // free
	ChiRefGs_fIc  = malloc(NfMax * sizeof *ChiRefGs_fIc); // free
	ChiRefGc_fIs  = malloc(NfMax * sizeof *ChiRefGc_fIs); // free
	ChiRefGc_fIc  = malloc(NfMax * sizeof *ChiRefGc_fIc); // free
	ChiRefCs_fIs  = malloc(NfMax * sizeof *ChiRefCs_fIs); // free
	ChiRefCs_fIc  = malloc(NfMax * sizeof *ChiRefCs_fIc); // free
	ChiRefCc_fIs  = malloc(NfMax * sizeof *ChiRefCc_fIs); // free
	ChiRefCc_fIc  = malloc(NfMax * sizeof *ChiRefCc_fIc); // free

	ChiGs_fIs     = malloc(NfMax * sizeof *ChiGs_fIs); // free
	ChiGs_fIc     = malloc(NfMax * sizeof *ChiGs_fIc); // free
	ChiGc_fIs     = malloc(NfMax * sizeof *ChiGc_fIs); // free
	ChiGc_fIc     = malloc(NfMax * sizeof *ChiGc_fIc); // free
	ChiCs_fIs     = malloc(NfMax * sizeof *ChiCs_fIs); // free
	ChiCs_fIc     = malloc(NfMax * sizeof *ChiCs_fIc); // free
	ChiCc_fIs     = malloc(NfMax * sizeof *ChiCc_fIs); // free
	ChiCc_fIc     = malloc(NfMax * sizeof *ChiCc_fIc); // free

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

	NFTypes = 1;
	BCoords_dEm1[0] = get_BCoords_dEm1(ELEMENT,0); // keep/free
	if (EType == WEDGE || EType == PYR) {
		NFTypes = 2;
		BCoords_dEm1[1] = get_BCoords_dEm1(ELEMENT,1);
	}

	for (IndFType = 0; IndFType < NFTypes; IndFType++)
		BCoords[IndFType] = malloc(sizeof *BCoords[IndFType]); // tbd


	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		B_Nve[IndFType]       = BCoords_dEm1[IndFType]->Nve;
		BCoords[IndFType]->Is = BCoords_dEm1[IndFType]->BCoords_Is;
		BCoords[IndFType]->Ic = BCoords_dEm1[IndFType]->BCoords_Ic;

		// Store Nfn* in a manner consistent with Nvn*
		for (P = 0; P <= PMax; P++) {
			NfnIs[P][IndFType] = BCoords_dEm1[IndFType]->NfnIs[P];
			NfnIc[P][IndFType] = BCoords_dEm1[IndFType]->NfnIc[P];
		}
		free(BCoords_dEm1[IndFType]->NfnIs);
		free(BCoords_dEm1[IndFType]->NfnIc);
	}

if (EType == PYR) {
printf("here\n");
for (IndFType = 0; IndFType < NFTypes; IndFType++) {
	printf("%d \n",B_Nve[IndFType]);

	for (P = 0; P <= PMax; P++) {
		printf("%d %d\n",NfnIs[P][IndFType],NfnIc[P][IndFType]);
		array_print_d(NfnIs[P][IndFType],B_Nve[IndFType],BCoords[IndFType]->Is[P],'C');
//		array_print_d(NfnIc[P][IndFType],B_Nve[IndFType],BCoords[IndFType]->Ic[P],'C');
	}
}
//exit(1);
}

	for (P = 0; P <= PMax; P++) {
		// VOLUME Operators
		cubature(&rst_vGc,&dummyPtr_d,&dummyPtr_ui[0],&NvnGc[P],&dummy_ui,0,PGc[P],         dE,NodeTypeG[Eclass]  );    free(dummyPtr_ui[0]); // free
		cubature(&rst_vCs,&dummyPtr_d,&dummyPtr_ui[0],&NvnCs[P],&dummy_ui,0,PCs[P][Eclass], dE,NodeTypeG[Eclass]  );    free(dummyPtr_ui[0]); // free
		cubature(&rst_vCc,&dummyPtr_d,&dummyPtr_ui[0],&NvnCc[P],&dummy_ui,0,PCc[P][Eclass], dE,NodeTypeG[Eclass]  );    free(dummyPtr_ui[0]); // free
		cubature(&rst_vIs,&dummyPtr_d,&dummyPtr_ui[0],&NvnIs[P],&dummy_ui,0,PIvs[P][Eclass],dE,NodeTypeIvs[P][Eclass]); free(dummyPtr_ui[0]); // free
		cubature(&rst_vIc,&dummyPtr_d,&dummyPtr_ui[0],&NvnIc[P],&dummy_ui,0,PIvc[P][Eclass],dE,NodeTypeIvc[P][Eclass]); free(dummyPtr_ui[0]); // free

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

		ChiRefGs_vGc = basis(PGs,           rst_vGc,NvnGc[P],&Nbf,dE); // free
		ChiRefGs_vCs = basis(PGs,           rst_vCs,NvnCs[P],&Nbf,dE); // free
		ChiRefGs_vIs = basis(PGs,           rst_vIs,NvnIs[P],&Nbf,dE); // free
		ChiRefGc_vP  = basis(PGc[P],        rst_vP, NvnP,    &Nbf,dE); // free
		ChiRefGc_vCc = basis(PGc[P],        rst_vCc,NvnCc[P],&Nbf,dE); // free
		ChiRefGc_vIc = basis(PGc[P],        rst_vIc,NvnIc[P],&Nbf,dE); // free
		ChiRefCs_vIs = basis(PCs[P][Eclass],rst_vIs,NvnIs[P],&Nbf,dE); // free
		ChiRefCs_vIc = basis(PCs[P][Eclass],rst_vIc,NvnIc[P],&Nbf,dE); // free
		ChiRefCc_vIs = basis(PCc[P][Eclass],rst_vIs,NvnIs[P],&Nbf,dE); // free
		ChiRefCc_vIc = basis(PCc[P][Eclass],rst_vIc,NvnIc[P],&Nbf,dE); // free

		ChiGs_vGc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vGc,TGs); // free
		ChiGs_vCs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vCs,TGs); // free
		ChiGs_vIs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vIs,TGs); // free
		ChiGc_vP  = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiRefGc_vP, TGc); // free
		ChiGc_vCc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vCc,TGc); // free
		ChiGc_vIc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vIc,TGc); // free
		ChiCs_vIs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIs,TCs); // free
		ChiCs_vIc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIc,TCs); // free
		ChiCc_vIs = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIs,TCc); // free
		ChiCc_vIc = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIc,TCc); // free

		GradChiRefGs_vCs = grad_basis(PGs,           rst_vCs,NvnCs[P],&Nbf,dE); // free
		GradChiRefGs_vIs = grad_basis(PGs,           rst_vIs,NvnIs[P],&Nbf,dE); // free
		GradChiRefGc_vCc = grad_basis(PGc[P],        rst_vCc,NvnCc[P],&Nbf,dE); // free
		GradChiRefGc_vIc = grad_basis(PGc[P],        rst_vIc,NvnIc[P],&Nbf,dE); // free
		GradChiRefCs_vCs = grad_basis(PCs[P][Eclass],rst_vCs,NvnCs[P],&Nbf,dE); // free
		GradChiRefCc_vCc = grad_basis(PCc[P][Eclass],rst_vCc,NvnCc[P],&Nbf,dE); // free

		for (dim = 0; dim < dE; dim++) {
			GradChiGs_vCs[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,GradChiRefGs_vCs[dim],TGs); // free
			GradChiGs_vIs[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnGs[0],NvnGs[0],1.0,GradChiRefGs_vIs[dim],TGs); // free
			GradChiGc_vCc[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vCc[dim],TGc); // free
			GradChiGc_vIc[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vIc[dim],TGc); // free
			GradChiCs_vCs[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiRefCs_vCs[dim],TCs); // free
			GradChiCc_vCc[dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiRefCc_vCc[dim],TCc); // free
		}

		// Returned Operators
		I_vGs_vGc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vGc,ChiInvGs_vGs); // keep
		I_vGs_vCs[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vCs,ChiInvGs_vGs); // keep
		I_vGs_vIs[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vIs,ChiInvGs_vGs); // keep
		I_vGc_vP[P]  = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiGc_vP, ChiInvGc_vGc); // keep
		I_vGc_vCc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vCc,ChiInvGc_vGc); // keep
		I_vGc_vIc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vIc,ChiInvGc_vGc); // keep
		I_vCs_vIs[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnCs[P],NvnCs[P],1.0,ChiCs_vIs,ChiInvCs_vCs); // keep
		I_vCs_vIc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnCs[P],NvnCs[P],1.0,ChiCs_vIc,ChiInvCs_vCs); // keep
		I_vCc_vIs[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnCc[P],NvnCc[P],1.0,ChiCc_vIs,ChiInvCc_vCc); // keep
		I_vCc_vIc[P] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnCc[P],NvnCc[P],1.0,ChiCc_vIc,ChiInvCc_vCc); // keep

		for (dim = 0; dim < dE; dim++) {
			D_vGs_vCs[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnGs[0],NvnGs[0],1.0,GradChiGs_vCs[dim],ChiInvGs_vGs); // keep
			D_vGs_vIs[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIs[P],NvnGs[0],NvnGs[0],1.0,GradChiGs_vIs[dim],ChiInvGs_vGs); // keep
			D_vGc_vCc[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vCc[dim],ChiInvGc_vGc); // keep
			D_vGc_vIc[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vIc[dim],ChiInvGc_vGc); // keep
			D_vCs_vCs[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiCs_vCs[dim],ChiInvCs_vCs); // keep
			D_vCc_vCc[P][dim] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiCc_vCc[dim],ChiInvCc_vCc); // keep
		}

		// FACET related operators
		// Add interpolation to P-1, P, P+1 only (ToBeDeleted)
		PbMin = P; PbMax = P;
/*
		if (P == 0) {
			PbMin = P;
			PbMax = P+1;
		} else if (P == PMax) {
			PbMin = P-1;
			PbMax = PMax;
		} else {
			PbMin = P-1;
			PbMax = P+1;
		}
*/
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			for (f = 0; f < Nf; f++) {
				IndFType = get_IndFType(Eclass,f);

				// Note: No additional index in VeF as only the conforming nodes are set up at the moment (ToBeModified)
				mm_CTN_d(Nfve[f],dE,Nve,&VeF[f*(NfrefMax*NfveMax*NveMax)],E_rst_vC,rst_vC);
//array_print_d(Nfve[IndFType],dE,rst_vC,'C');

				rst_fIs[f] = mm_Alloc_d(CblasColMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],dE,B_Nve[IndFType],1.0,BCoords[IndFType]->Is[Pb],rst_vC); // free
				rst_fIc[f] = mm_Alloc_d(CblasColMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],dE,B_Nve[IndFType],1.0,BCoords[IndFType]->Ic[Pb],rst_vC); // free

if (EType == PYR) {
if (f == 0 && P == 0)
	array_print_d(5,3,E_rst_vC,'C');

	printf("%d %d %d\n",Pb,f,IndFType);
	array_print_d(NfnIs[Pb][IndFType],dE,rst_fIs[f],'C');
//	array_print_d(NfnIc[Pb][IndFType],dE,rst_fIc[f],'C');
if (f == 4 && P == 1) {
	// SEEMS THERE IS A PROBLEM WITH THE MM_ALLOC_D ABOVE.
	array_print_d(NfnIs[Pb][IndFType],B_Nve[IndFType],BCoords[IndFType]->Is[Pb],'C');
	array_print_d(NfnIs[Pb][IndFType],5,BCoords[IndFType]->Is[Pb],'C');

	exit(1);
}
}
				ChiRefGs_fIs[f] = basis(PGs           ,rst_fIs[f],NfnIs[Pb][IndFType],&Nbf,dE); // free
				ChiRefGs_fIc[f] = basis(PGs           ,rst_fIc[f],NfnIc[Pb][IndFType],&Nbf,dE); // free
				ChiRefGc_fIs[f] = basis(PGc[P]        ,rst_fIs[f],NfnIs[Pb][IndFType],&Nbf,dE); // free
				ChiRefGc_fIc[f] = basis(PGc[P]        ,rst_fIc[f],NfnIc[Pb][IndFType],&Nbf,dE); // free
				ChiRefCs_fIs[f] = basis(PCs[P][Eclass],rst_fIs[f],NfnIs[Pb][IndFType],&Nbf,dE); // free
				ChiRefCs_fIc[f] = basis(PCs[P][Eclass],rst_fIc[f],NfnIc[Pb][IndFType],&Nbf,dE); // free
				ChiRefCc_fIs[f] = basis(PCc[P][Eclass],rst_fIs[f],NfnIs[Pb][IndFType],&Nbf,dE); // free
				ChiRefCc_fIc[f] = basis(PCc[P][Eclass],rst_fIc[f],NfnIc[Pb][IndFType],&Nbf,dE); // free

				ChiGs_fIs[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiRefGs_fIs[f],TGs); // free
				ChiGs_fIc[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiRefGs_fIc[f],TGs); // free
				ChiGc_fIs[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIs[f],TGc); // free
				ChiGc_fIc[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIc[f],TGc); // free
				ChiCs_fIs[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIs[f],TCs); // free
				ChiCs_fIc[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIc[f],TCs); // free
				ChiCc_fIs[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIs[f],TCc); // free
				ChiCc_fIc[f] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIc[f],TCc); // free

				I_vGs_fIs[P][Pb][f*7] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiGs_fIs[f],ChiInvGs_vGs); // keep
				I_vGs_fIc[P][Pb][f*8] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiGs_fIc[f],ChiInvGs_vGs); // keep
				I_vGc_fIs[P][Pb][f*9] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIs[f],ChiInvGc_vGc); // keep
				I_vGc_fIc[P][Pb][f*9] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIc[f],ChiInvGc_vGc); // keep
				I_vCs_fIs[P][Pb][f*9] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIs[f],ChiInvCs_vCs); // keep
				I_vCs_fIc[P][Pb][f*9] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIc[f],ChiInvCs_vCs); // keep
				I_vCc_fIs[P][Pb][f*9] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIs[f],ChiInvCc_vCc); // keep
				I_vCc_fIc[P][Pb][f*9] = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIc[f],ChiInvCc_vCc); // keep
			}

			for (f = 0; f < Nf; f++) {
				free(rst_fIs[f]);
				free(rst_fIc[f]);

				free(ChiRefGs_fIs[f]);
				free(ChiRefGs_fIc[f]);
				free(ChiRefGc_fIs[f]);
				free(ChiRefGc_fIc[f]);
				free(ChiRefCs_fIs[f]);
				free(ChiRefCs_fIc[f]);
				free(ChiRefCc_fIs[f]);
				free(ChiRefCc_fIc[f]);

				free(ChiGs_fIs[f]);
				free(ChiGs_fIc[f]);
				free(ChiGc_fIs[f]);
				free(ChiGc_fIc[f]);
				free(ChiCs_fIs[f]);
				free(ChiCs_fIc[f]);
				free(ChiCc_fIs[f]);
				free(ChiCc_fIc[f]);
			}
		}

		free(rst_vGc);
		free(rst_vCs);
		free(rst_vCc);
		free(rst_vIs);
		free(rst_vIc);

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
		free(ChiRefGs_vIs);
		free(ChiRefGc_vP);
		free(ChiRefGc_vCc);
		free(ChiRefGc_vIc);
		free(ChiRefCs_vIs);
		free(ChiRefCs_vIc);
		free(ChiRefCc_vIs);
		free(ChiRefCc_vIc);

		free(ChiGs_vGc);
		free(ChiGs_vCs);
		free(ChiGs_vIs);
		free(ChiGc_vP);
		free(ChiGc_vCc);
		free(ChiGc_vIc);
		free(ChiCs_vIs);
		free(ChiCs_vIc);
		free(ChiCc_vIs);
		free(ChiCc_vIc);

		array_free2_d(dE,GradChiRefGs_vCs);
		array_free2_d(dE,GradChiRefGs_vIs);
		array_free2_d(dE,GradChiRefGc_vCc);
		array_free2_d(dE,GradChiRefGc_vIc);
		array_free2_d(dE,GradChiRefCs_vCs);
		array_free2_d(dE,GradChiRefCc_vCc);

		for (dim = 0; dim < dE; dim++) {
			free(GradChiGs_vCs[dim]);
			free(GradChiGs_vIs[dim]);
			free(GradChiGc_vCc[dim]);
			free(GradChiGc_vIc[dim]);
			free(GradChiCs_vCs[dim]);
			free(GradChiCc_vCc[dim]);
		}
	}

	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		array_free2_d(PMax+1,BCoords[IndFType]->Is);
		array_free2_d(PMax+1,BCoords[IndFType]->Ic);
		free(BCoords[IndFType]);
		free(BCoords_dEm1[IndFType]);
	}

	free(E_rst_vC);
	free(rst_vC);
	free(rst_vP);
	free(rst_fIs);
	free(rst_fIc);

	free(ChiInvGs_vGs);
	free(TGs);

	free(GradChiGs_vCs);
	free(GradChiGs_vIs);
	free(GradChiGc_vCc);
	free(GradChiGc_vIc);
	free(GradChiCs_vCs);
	free(GradChiCc_vCc);

	free(ChiRefGs_fIs);
	free(ChiRefGs_fIc);
	free(ChiRefGc_fIs);
	free(ChiRefGc_fIc);
	free(ChiRefCs_fIs);
	free(ChiRefCs_fIc);
	free(ChiRefCc_fIs);
	free(ChiRefCc_fIc);

	free(ChiGs_fIs);
	free(ChiGs_fIc);
	free(ChiGc_fIs);
	free(ChiGc_fIc);
	free(ChiCs_fIs);
	free(ChiCs_fIc);
	free(ChiCc_fIs);
	free(ChiCc_fIc);
}

void setup_ELEMENT_VeF(const unsigned int EType)
{
	// Initialize DB Parameters
	unsigned int NveMax   = DB.NveMax,
	             NfveMax  = DB.NfveMax,
	             NfrefMax = DB.NfrefMax;

	// Standard datatypes
	unsigned int i, j, k, l, iMax, jMax, kMax, lMax, iStep,
	             f, Nve, *Nfve, *Nfref, Nf, *VeFcon;
	double *VeF;

	struct S_ELEMENT *ELEMENT;

	// silence
	VeF = NULL;

	ELEMENT = get_ELEMENT_type(EType);
	Nve    = ELEMENT->Nve;
	Nfve   = ELEMENT->Nfve;
	Nf     = ELEMENT->Nf;
	VeFcon = ELEMENT->VeFcon;

	// Set Nfref
	switch(EType) {
	case LINE:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = 1;
		break;
	case TRI:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = 3;
		break;
	case QUAD:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = 3;
		break;
	case TET:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = 5;
		break;
	case HEX:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = 9;
		break;
	case WEDGE:
		for (f = 0; f < 3; f++)
			ELEMENT->Nfref[f] = 9;
		for (f = 3; f < 6; f++)
			ELEMENT->Nfref[f] = 5;
		break;
	case PYR:
		for (f = 0; f < 4; f++)
			ELEMENT->Nfref[f] = 5;
		for (f = 4; f < 5; f++)
			ELEMENT->Nfref[f] = 9;
		break;
	default:
		printf("Error: Unsupported EType in setup_ELEMENT_VeF.\n"), exit(1);
		break;
	}

	Nfref   = ELEMENT->Nfref;

	unsigned int size_VeF = 0;
	double VeFref_LINE[12]  = {1.0 , 0.0 ,
	                           0.0 , 1.0 ,
	                           1.0 , 0.0 ,
	                           0.5 , 0.5 ,
	                           0.5 , 0.5 ,
	                           0.0 , 1.0 },
	       VeFref_TRI[45]   = {1.0 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 ,
	                           0.5 , 0.5 , 0.0 ,
	                           0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 ,
	                           0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 1.0 ,
	                           0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 },
	       VeFref_QUAD[144] = {1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 1.0 };

	switch(EType) {
	case LINE:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		VeF[0] = 1.0; VeF[1] = 0.0;
		VeF[2] = 0.0; VeF[3] = 1.0;

		break;
	case TRI:
		/* fall through */
	case QUAD:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*4+k]] = VeFref_LINE[j*kMax+k];
		}}}

		break;
	case TET:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*4+k]] = VeFref_TRI[j*kMax+k];
		}}}

		break;
	case HEX:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*4+k]] = VeFref_QUAD[j*kMax+k];
		}}}

		break;
	case WEDGE:
		/* fall through */
	case PYR:
		for (i = 0; i < Nf; i++)
			size_VeF += Nfve[i]*Nfref[i];
		size_VeF *= Nve;

		VeF = calloc(size_VeF , sizeof *VeF); // free

		for (i = 0, iMax = Nf; i < iMax; i++) {
			for (j = iStep = 0; j < i; j++)
				iStep += Nfref[j]*Nfve[j];
			iStep *= Nve;

			for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
			for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
				if (Nfve[i] == 3)
					VeF[iStep+j*Nve+VeFcon[i*4+k]] = VeFref_TRI[j*kMax+k];
				else if (Nfve[i] == 4)
					VeF[iStep+j*Nve+VeFcon[i*4+k]] = VeFref_QUAD[j*kMax+k];
			}
		}}

		break;
	default:
		// Already caught error above.
		break;
	}

	for (i = 0, iMax = Nf;       i < iMax; i++) {
		for (j = iStep = 0; j < i; j++)
			iStep += Nfref[j]*Nfve[j];
		iStep *= Nve;

		for (j = 0, jMax = Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];  k < kMax; k++) {
		for (l = 0, lMax = Nve;      l < lMax; l++) {
			ELEMENT->VeF[i*(NfrefMax*NfveMax*NveMax)+j*(Nfve[i]*Nve)+k*(Nve)+l] = VeF[iStep+j*(Nfve[i]*Nve)+k*(Nve)+l];
		}
	}}}

//array_print_d(Nfref[i]*Nfve[i],Nve,&ELEMENT->VeF[i*(NfrefMax*NfveMax*NveMax)],'R');

	free(VeF);
}

void setup_operators(void)
{
	// Initialize DB Parameters

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int EType;



	// LINE (Includes TP Class)
	EType = LINE;

	setup_ELEMENT_VeF(EType);
	setup_ELEMENT_plotting(EType);
	setup_ELEMENT_normals(EType);
	setup_ELEMENT_operators(EType);

	// QUAD
	EType = QUAD;

	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
	}

	// HEX
	EType = HEX;

	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
	}

	// TRI
	EType = TRI;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);

	}

	// TET
	EType = TET;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);

	}

	// PYR
	EType = PYR;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);

	}

	// WEDGE
	EType = WEDGE;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);

	}
}
