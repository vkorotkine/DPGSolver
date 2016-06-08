// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

#include "mkl.h"

/*
 *	Purpose:
 *		Set up operators to be used throughout the code.
 *
 *	Comments:
 *		Different operators are set up depending on which higher-dimensional element dependencies are present as well as
 *		whether adaptivity is desired:
 *			LINE  : Standard operators + additional SF operators if QUADs are present (standard assembly)
 *			TRI   : Standard operators + additional SF operators if WEDGEs are present (standard assembly)
 *			QUAD  : Standard operators (TP asssembly)
 *			TET   : Standard operators (standard assembly)
 *			HEX   : Standard operators (TP asssembly)
 *			WEDGE : Standard operators (TP asssembly)
 *			PYR   : Standard operators (standard assembly)
 *
 *			This is done so that the only operators computed are those which are used. The additional SF operators are
 *			needed to interpolate from VOLUME to VOLUME nodes in components of the TP operator assembly for higher
 *			dimensional differentiation and FACET operators.
 *			It must be noted, however, that the data structure is such that the most general case can be handled (3D
 *			mixed meshes using hp adaptivity and sum factorized operators) which results in notably excessive operator
 *			indexing in the much simpler cases (such as in 1D). The indexing convention is as follows:
 *				OP_v*_(v/f)*[1][2][3] : (OP)erator from order [1] to order [2] with (v)olume/(f)acet refinement index
 *				                        [3] where potential (OP)erators and '*'s are defined in the notation below.
 *			If sum factorization has not yet broken-even for a certain order, the associated operators are freed here.
 *			(ToBeModified: DON'T FORGET TO DO THIS. Remove the sum factorized output_to_paraview function as well. This
 *			can be done when sum factorization is included in the solver).
 *
 *		Standard (i.e. non sum-factorized) operators are used for all functions which are not performance critical in
 *		order to improve code readability. (ToBeModified)
 *		Ensure that operators for hp refinement are only stored when refinement is enabled (ToBeDeleted).
 *		For the PYR element rst_vGs != E_rst_vC because of the TP structures of the TP nodes and the rotational symmetry
 *		ordering of the PYR nodes in layers of 't'. This means that special consideration must be made if attempting to
 *		use the rotational symmetry of rst_vGs PYR nodes.
 *
 *	Notation:
 *
 *	General:
 *
 *		Eclass : (E)lement (class)
 *		         Options: TP, SI, PYR, WEDGE
 *
 *		N(1)n(2)(3)[4][5] : (N)umber of (1) (n)odes of (2) type on elements which are (3) of order [4] belonging to
 *		                    element class [5]
 *		                    (1): (v)olume, (f)acet
 *		                    (2): [P]lotting, [G]eometry, [C]ofactor, [I]ntegration (ToBeModified)
 *		                    (3): (s)traight, (c)urved
 *		  Example: NfnIs[1][0] == (N)umber of (f)acet (n)odes for (I)ntegration on (straight) P1 facets of Eclass 0.
 *
 *		rst_(1)(2)(3) : coordinates on the reference element (rst) of (1) nodes of (2) type
 *		                (1): (v)olume, (f)acet
 *		                (2): (C)orner
 *		                (3): (s)traight, (c)urved (optional)
 *
 *		(Grad)Chi(Ref)(1)(2)_(3)(4)(5)[6] : (Grad)ient (optional) (Ref)erence (optional) Basis functions (Chi) of type
 *		                                    (1) which are (2) evaluated at (3) nodes of (4) type which are (5) of order [6]
 *		                                    (1/4): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *		                                    (2/5): (s)traight, (c)urved
 *		                                    (3): (v)olume, (f)acet
 *
 *		_(1)(2) subscripts for operators indicate that this is a (1) operator, operating on VOLUME nodes and
 *		transferring to (2) nodes.
 *			Options (1)(2) : VV, FV, FF
 *
 *	setup_ELEMENT_plotting:
 *		connectivity  : connectivity between nodes to form P1 sub-elements (see test_imp_plotting for visualization)
 *		connect_types : VTK element type numbering for each sub-element
 *		                Note: This is not trivial as TETs and PYRs are refined into a combination of TETs and PYRs.
 *		connect_NE    : (N)umber of sub (E)lements in connectivity array
 *
 *	setup_ELEMENT_normals:
 *		Theta_[] : Angles for conversion between reference and facet coordinates (Zwanenburg(2016): Table 14)
 *		           Options: eta, zeta
 *		nr       : (n)ormal vector components for the (r)eference element
 *
 *	get_BCoords_dEm1: get (B)ary(C)entric (Coord)inates of (d)imension (E)lement (m)inus 1
 *
 *		Note: See cubature_PYR for explanation of the method to obtain the barycentric coordinates while noting that
 *		      using the P1 reference basis provides the partition of unity, linear precision and additional arbitrary
 *		      conditions needed for the construction.
 *
 *	setup_ELEMENT_operators/setup_TP_operators:
 *		w_(1)(2)(3) : Cubature (w)eights (1) nodes of (2) type which are (3)
 *		             (1): (v)olume, (f)acet
 *		             (2): (I)ntegration
 *		             (3): (s)traight, (c)urved
 *		I(1)(2) : (I)dentity matrix of type (1) which is (2)
 *		          (1): (G)eometry, (C)ofactor, (S)olution
 *		          (2): (s)traight, (c)urved
 *		T(1)(2) : (T)ransformation matrix of type (1) which is (2). (Zwanenburg(2016): eq. 2.14)
 *		          (1): (G)eometry, (C)ofactor, (S)olution
 *		          (2): (s)traight, (c)urved
 *		I_(1)(2)(3)_(4)(5)(6) : (I)nterpolation operator from (1) nodes of type (2) which are (3) to (4) nodes of type
 *		                        (5) which are (6)
 *		                        (1/4): (v)olume, (f)acet
 *		                        (2/5): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *		                        (3/6): (s)traight, (c)urved
 *		  Note: Interpolation operators to FACETs may have up to 4 levels of dereferencing *I_x_x[1][2][3]:
 *		        [0] : The operator is a pointer to  matrix
 *		        [1] : Order of the element from which you are interpolating
 *		        [2] : Order of the element to which you are interpolating
 *		        [3] : Index of operator (e.g. there is a need for many operators interpolating to different parts of the
 *		                                 VOLUME FACETs)
 *		D_(1)(2)(3)_(4)(5)(6) : (D)ifferentiation + interpolation operator from (1) nodes of type (2) which are (3) to
 *		                        (4) nodes of type (5) which are (6)
 *		                        (1/4): (v)olume, (f)acet
 *		                        (2/5): (P)lotting, (G)eometry, (C)ofactor, (I)ntegration, (S)olution
 *		                        (3/6): (s)traight, (c)urved
 *		  Note: Differentiation + interpolation operators may have up to 3 levels of dereferencing *D_x_x[1][2]:
 *		        [0] : The operator is a pointer to  matrix
 *		        [1] : Order of the element from which you are interpolating with differentiation
 *		        [2] : Index of dimension of the operator
 *
 *	setup_ELEMENT_VeF:
 *		VeF : (Ve)rtex to (F)acet operator used to project VOLUME vertex nodes to FACET vertex nodes for all supported
 *		      h-refinements.
 *
 *	sf_assemble_d/get_sf_parameters:
 *		NIn  : (N)umber of (In)puts (i.e. number of columns of the lower dimensional operator)
 *		NOut : (N)umber of (Out)puts (i.e. number of rows of the lower dimensional operator)
 *		OP   : (OP)erator
 *
 *	Returned operators : [range used depending on Adapt]
 *		- For maximum array size, see memory_constructors.c.
 *		- If no range is indicated, look at range of previous operator.
 *
 *		- (*) indicates a special case for SF operators.
 *		- [^] is a placeholder.
 *		- Further notation simplifications for table below:
 *			*    : 0:NVREFSFMAX
 *			^    : 0:NFMAX
 *			^^   : 0:NFREFMAX*NFMAX
 *
 *			rP   : (r)ange P  == 0:PMax
 *			P(P) : P as a function of P (i.e. for a range in field 1 of [0:PMax] -> [0][0], [1][1], ..., [P][P])
 *			rPb  : (r)ange Pb (always a function of P) == PbMin:PbMax (PbMin = min(P-1,0), PbMax = max(P+1,PMax))
 *			rd   : (r)ange (d)imension == 0:d-1
 *
 *		                   Adapt
 *		Returned Op (SF) | 1               2                   3               4
 *		                 |
 *		NvnGs            | [0]             [0]                 [0]             [0]
 *		NvnGc            | [P]             [rP]                [P]             [rP]
 *		NvnCs            |
 *		NvnCc            |
 *		NvnIs            |
 *		NvnIc            |
 *		NfnIs            | [P][0:1]        [rP][0:1]           [P][2]          [rP][0:1]
 *		NfnIc            |
 *		                 |
 *		ChiS_vP          | [P][PP][0]      [rP][PP][0]         [P][PP][0]      [rP][PP][0]
 *		ChiS_vIs     (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rPb][*]
 *		ChiS_vIc     (*) |
 *    	ChiInvS_vS       | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		                 |
 *		ICs              | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		ICc              |
 *		                 |
 *		I_vGs_vP         | [0][PP][0]      [0][PP][0]          [0][PP][0]      [0][PP][0]
 *		I_vGs_vGc        | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		I_vGs_vCs        |
 *		I_vGs_vIs    (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rPb][*]
 *		I_vGs_vIc    (*) |
 *		I_vGs_vS         | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		I_vGc_vP         | [0][PP][0]      [rP][PP][0]         [0][PP][0]      [rP][PP][0]
 *		I_vGc_vCc        | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		I_vGc_vIs    (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rPb][*]
 *		I_vGc_vIc    (*) |
 *		I_vGc_vS         | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		I_vCs_vIs    (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rPb][*]
 *		I_vCs_vIc    (*) |
 *		I_vCc_vIs    (*) |
 *		I_vCc_vIc    (*) |
 *		                 |
 *		D_vGs_vCs        | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		D_vGs_vIs        |
 *		D_vGc_vCc        |
 *		D_vGc_vIc        |
 *		D_vCs_vCs        |
 *		D_vCc_vCc        |
 *		                 |
 *    	ChiS_fIs         | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rPb][^^]
 *    	ChiS_fIc         |
 *		                 |
 *    	I_vGs_fIs        | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rPb][^^]
 *    	I_vGs_fIc        |
 *    	I_vGc_fIs        |
 *    	I_vGc_fIc        |
 *    	I_vCs_fIs        |
 *    	I_vCs_fIc        |
 *    	I_vCc_fIs        |
 *    	I_vCc_fIc        |
 *		                 |
 *		Is_Weak_VV       | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		Ic_Weak_VV       |
 *		Ds_Weak_VV       | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		Dc_Weak_VV       |
 *
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
	double       **w_fIs, **w_fIc, **BCoords_Is, **BCoords_Ic;
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
		printf("Error: Unsupported type (%d) in select_functions.\n",type), exit(1);
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
		rst_vC[0*Nve+2] = -1.0; rst_vC[1*Nve+2] =  1.0;
		rst_vC[0*Nve+3] =  1.0; rst_vC[1*Nve+3] =  1.0;
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
		rst_vC[0*Nve+2] = -1.0; rst_vC[1*Nve+2] =  1.0; rst_vC[2*Nve+2] = -1.0;
		rst_vC[0*Nve+3] =  1.0; rst_vC[1*Nve+3] =  1.0; rst_vC[2*Nve+3] = -1.0;
		rst_vC[0*Nve+4] = -1.0; rst_vC[1*Nve+4] = -1.0; rst_vC[2*Nve+4] =  1.0;
		rst_vC[0*Nve+5] =  1.0; rst_vC[1*Nve+5] = -1.0; rst_vC[2*Nve+5] =  1.0;
		rst_vC[0*Nve+6] = -1.0; rst_vC[1*Nve+6] =  1.0; rst_vC[2*Nve+6] =  1.0;
		rst_vC[0*Nve+7] =  1.0; rst_vC[1*Nve+7] =  1.0; rst_vC[2*Nve+7] =  1.0;
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
		rst_vC[0*Nve+2] = -1.0; rst_vC[1*Nve+2] =  1.0; rst_vC[2*Nve+2] = -1.0/5.0*sqrt(2.0);
		rst_vC[0*Nve+3] =  1.0; rst_vC[1*Nve+3] =  1.0; rst_vC[2*Nve+3] = -1.0/5.0*sqrt(2.0);
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
	double       **w_fIs, **w_fIc,
	             **BCoords_Is, **BCoords_Ic, *one;

	struct S_BCOORDS *BCoords_dEm1;
	struct S_ELEMENT *ELEMENT_F;

	BCoords_dEm1 = malloc(sizeof *BCoords_dEm1);      // keep (requires external free)
	NfnIs        = malloc(NP * sizeof *(NfnIs));      // keep (requires external free)
	NfnIc        = malloc(NP * sizeof *(NfnIc));      // keep (requires external free)
	w_fIs        = malloc(NP * sizeof *(w_fIs));      // keep (requires external free)
	w_fIc        = malloc(NP * sizeof *(w_fIc));      // keep (requires external free)
	BCoords_Is   = malloc(NP * sizeof *(BCoords_Is)); // keep (requires external free)
	BCoords_Ic   = malloc(NP * sizeof *(BCoords_Ic)); // keep (requires external free)

	one = malloc(1 * sizeof *one); // free
	one[0] = 1.0;

	EType = ELEMENT->type;
	if (EType == LINE) {
// fix this to be consistent with treatment of other elements if possible (ToBeDeleted)
		ELEMENT_F = get_ELEMENT_type(POINT);

		Nve = 1;
		for (P = 0; P <= PMax; P++) {
			NfnIs[P] = 1;
			NfnIc[P] = 1;

			w_fIs[P] = mm_Alloc_d(CBCM,CBT,CBT,1,1,1,1.0,one,one); // keep
			w_fIc[P] = mm_Alloc_d(CBCM,CBT,CBT,1,1,1,1.0,one,one); // keep

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
			cubature(&rst_fIs,&w_fIs[P],&dummyPtr_ui,&NfnIs[P],&dummy_ui,1,PIfs[P][Eclass],dE,NodeTypeIfs[P][Eclass]); free(dummyPtr_ui); // free
			cubature(&rst_fIc,&w_fIc[P],&dummyPtr_ui,&NfnIc[P],&dummy_ui,1,PIfc[P][Eclass],dE,NodeTypeIfc[P][Eclass]); free(dummyPtr_ui); // free

			ChiRefGs_fIs = basis(1,rst_fIs,NfnIs[P],&Nbf,dE); // free
			ChiRefGs_fIc = basis(1,rst_fIc,NfnIc[P],&Nbf,dE); // free

			BCoords_Is[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnIs[P],NfnGs,NfnGs,1.0,ChiRefGs_fIs,ChiRefInvGs_vGs); // keep
			BCoords_Ic[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnIc[P],NfnGs,NfnGs,1.0,ChiRefGs_fIc,ChiRefInvGs_vGs); // keep

			free(rst_fIs);
			free(rst_fIc);

			free(ChiRefGs_fIs);
			free(ChiRefGs_fIc);
		}
		free(ChiRefInvGs_vGs);
	}

	free(one);

	BCoords_dEm1->Nve   = Nve;
	BCoords_dEm1->NfnIs = NfnIs;
	BCoords_dEm1->NfnIc = NfnIc;
	BCoords_dEm1->w_fIs = w_fIs;
	BCoords_dEm1->w_fIc = w_fIc;
	BCoords_dEm1->BCoords_Is = BCoords_Is;
	BCoords_dEm1->BCoords_Ic = BCoords_Ic;

	return BCoords_dEm1;
}

static void setup_ELEMENT_VeV(const unsigned int EType)
{
	// Standard datatypes
	unsigned int i, j, jMax, VeVrefInd, Nve, *Nvve;
	double       **VeV;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	Nve  = ELEMENT->Nve;
	Nvve = ELEMENT->Nvve;
	VeV  = ELEMENT->VeV;

	double VeVref_LINE[12] = { 1.0, 0.0 ,
		                       0.0, 1.0 ,
		                       1.0, 0.0 ,
		                       0.5, 0.5 ,
		                       0.5, 0.5 ,
		                       0.0, 1.0 },
	       VeVref_TRI[45]   = {1.0 , 0.0 , 0.0 ,
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
	       VeVref_TET[120]  = {1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 1.0 },
	       VeVref_PYR[255]  = {1.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 , 0.0 ,
	                           0.5 , 0.5 , 0.0 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.5 , 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.5 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.5 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 0.0 , 1.0 };


	switch(EType) {
	case LINE:
		for (i = 0; i < NREFMAXLINE; i++)
			Nvve[i] = 2;

		VeVrefInd = 0;
		for (i = 0; i < NREFMAXLINE; i++) {
			jMax = Nve*Nvve[i];
			if (i)
				VeVrefInd += jMax;
			VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
			for (j = 0; j < jMax; j++)
				VeV[i][j] = VeVref_LINE[VeVrefInd+j];
		}
		break;
	case TRI:
		for (i = 0; i < NREFMAXTRI; i++)
			Nvve[i] = 3;

		VeVrefInd = 0;
		for (i = 0; i < NREFMAXTRI; i++) {
			jMax = Nve*Nvve[i];
			if (i)
				VeVrefInd += jMax;
			VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
			for (j = 0; j < jMax; j++)
				VeV[i][j] = VeVref_TRI[VeVrefInd+j];
		}
		break;
	case TET:
		// Original TET, 4 TETs, 2 PYRs
		for (i = 0; i < 5; i++)
			Nvve[i] = Nve;
		for (i = 5; i < NREFMAXTET; i++)
			Nvve[i] = 5;

		VeVrefInd = 0;
		for (i = 0; i < NREFMAXTET; i++) {
			jMax = Nve*Nvve[i];
			if (i)
				VeVrefInd += jMax;
			VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
			for (j = 0; j < jMax; j++)
				VeV[i][j] = VeVref_TET[VeVrefInd+j];
		}
		break;
	case PYR:
		// Original PYR, 4 PYRs, 4 TETs, 2 PYRs
		for (i = 0; i < 5; i++)
			Nvve[i] = Nve;
		for (i = 5; i < 9; i++)
			Nvve[i] = 4;
		for (i = 9; i < NREFMAXPYR; i++)
			Nvve[i] = Nve;

		VeVrefInd = 0;
		for (i = 0; i < NREFMAXPYR; i++) {
			jMax = Nve*Nvve[i];
			if (i)
				VeVrefInd += jMax;
			VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
			for (j = 0; j < jMax; j++)
				VeV[i][j] = VeVref_PYR[VeVrefInd+j];
		}
		break;
	case QUAD:
	case HEX:
	case WEDGE:
	default:
		printf("Error: Unsupported EType in setup_ELEMENT_VeV.\n"), exit(1);
		break;
	}
}

static void setup_ELEMENT_operators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS, **NfnIs, **NfnIc;
	double       **w_vIs, **w_vIc,
	             ****ChiS_vP, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS,
	             ****ICs, ****ICc,
	             ****I_vGs_vP, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vIs, ****I_vGs_vIc, ****I_vGs_vS,
	             ****I_vGc_vP,                ****I_vGc_vCc, ****I_vGc_vIs, ****I_vGc_vIc, ****I_vGc_vS,
	             ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vIs, ****I_vCc_vIc,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fIs, ****ChiS_fIc,
	             ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fIs, ****I_vCc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FF, ****Ic_Weak_FF,
	             *****Ds_Weak_VV, *****Dc_Weak_VV;

	// Initialize DB Parameters
	unsigned int Adapt       = DB.Adapt,
	             PMax        = DB.PMax,
	             NP          = DB.NP,
	             PGlobal     = DB.PGlobal,
	             PGs         = DB.PGs,
	             *PGc        = DB.PGc,
	             **PCs       = DB.PCs,
	             **PCc       = DB.PCc,
	             **PIvs      = DB.PIvs,
	             **PIvc      = DB.PIvc,
	             PP          = DB.PP,
	             EFE         = DB.EFE;

	char         *BasisType     = DB.BasisType,
	             **NodeTypeG    = DB.NodeTypeG,
	             ***NodeTypeIvs = DB.NodeTypeIvs,
	             ***NodeTypeIvc = DB.NodeTypeIvc,
	             ***NodeTypeS   = DB.NodeTypeS;

	// Standard datatypes
	unsigned int dim, dE, P, f, fh, Vf, IndFType, PSMin, PSMax, Pb, PbMin, PbMax, *fhMax,
	             Nve, Nf, Nbf, Eclass, NFTypes, Nvref, vref, vrefSF, NvrefSF, *Nfref, *ones_Nf,
	             NvnP,
	             B_Nve[2], *Nfve, *Nvve,
	             dummy_ui, *dummyPtr_ui[2];
	double       *E_rst_vC, *rst_vC, **VeF, **VeV,
	             *rst_vP, *rst_vGs, *rst_vGc, *rst_vCs, *rst_vCc, **rst_vIs, **rst_vIc, *rst_vS, *rst_fIs, *rst_fIc,
	             *diag_w_vIs, *diag_w_vIc,
	             ***w_fIs, ***w_fIc, *diag_w_fIs, *diag_w_fIc,
	             *IGs, *IGc, *IS,
	             *TGs, *TGc, *TCs, *TCc, *TS,
	             *ChiRefGs_vGs, *ChiRefGc_vGc, *ChiRefCs_vCs, *ChiRefCc_vCc, *ChiRefS_vS,
	             *ChiGs_vGs,    *ChiGc_vGc,    *ChiCs_vCs,    *ChiCc_vCc,    *ChiS_vS,
	             *ChiRefInvGs_vGs, *ChiRefInvGc_vGc, *ChiRefInvCs_vCs, *ChiRefInvCc_vCc, *ChiRefInvS_vS,
	             *ChiInvGs_vGs,    *ChiInvGc_vGc,    *ChiInvCs_vCs,    *ChiInvCc_vCc,
	             *ChiRefGs_vP, *ChiRefGs_vGc, *ChiRefGs_vCs, *ChiRefGs_vIs, *ChiRefGs_vIc, *ChiRefGs_vS,
	             *ChiRefGs_fIs, *ChiRefGs_fIc,
	             *ChiRefGc_vP,                *ChiRefGc_vCc, *ChiRefGc_vIs, *ChiRefGc_vIc, *ChiRefGc_vS,
	             *ChiRefGc_fIs, *ChiRefGc_fIc,
	             *ChiRefCs_vIs, *ChiRefCs_vIc, *ChiRefCs_fIs, *ChiRefCs_fIc,
	             *ChiRefCc_vIs, *ChiRefCc_vIc, *ChiRefCc_fIs, *ChiRefCc_fIc,
	             *ChiRefS_vP, *ChiRefS_vIs, *ChiRefS_vIc, *ChiRefS_fIs, *ChiRefS_fIc,
	             *ChiGs_vP, *ChiGs_vGc, *ChiGs_vCs, *ChiGs_vIs, *ChiGs_vIc, *ChiGs_vS,
	             *ChiGs_fIs, *ChiGs_fIc,
	             *ChiGc_vP,             *ChiGc_vCc, *ChiGc_vIs, *ChiGc_vIc, *ChiGc_vS,
	             *ChiGc_fIs, *ChiGc_fIc,
	             *ChiCs_vIs, *ChiCs_vIc, *ChiCs_fIs, *ChiCs_fIc,
	             *ChiCc_vIs, *ChiCc_vIc, *ChiCc_fIs, *ChiCc_fIc,
	             **GradChiRefGs_vCs, **GradChiRefGs_vIs,
	             **GradChiRefGc_vCc, **GradChiRefGc_vIc,
	             **GradChiRefCs_vCs, **GradChiRefCc_vCc,
	             **GradChiRefS_vIs,  **GradChiRefS_vIc,
	             **GradChiGs_vCs, **GradChiGs_vIs,
	             **GradChiGc_vCc, **GradChiGc_vIc,
	             **GradChiCs_vCs, **GradChiCc_vCc,
	             **GradChiS_vIs,  **GradChiS_vIc,
	             *dummyPtr_d;

	struct BCoords {
		double **Is, **Ic;
	} *BCoords_F[2], *BCoords_V;

	struct S_BCOORDS *BCoords_dEm1[2];
	struct S_ELEMENT *ELEMENT;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	// silence
	PSMin = PSMax = 0;
	ChiGs_vGs = NULL; ChiGc_vGc = NULL;
	ChiCs_vCs = NULL; ChiCc_vCc = NULL;
	ChiS_vS   = NULL;

	ELEMENT = get_ELEMENT_type(EType);

	// No need to consider the second Eclass as WEDGE basis functions will be built through a combination of lower
	// dimensional operators.
	Eclass = get_Eclass(EType);
	setup_ELEMENT_VeV(EType);

	dE    = ELEMENT->d;
	Nve   = ELEMENT->Nve;
	Nf    = ELEMENT->Nf;
	Nfref = ELEMENT->Nfref;
	Nfve  = ELEMENT->Nfve;
	VeF   = ELEMENT->VeF;
	VeV   = ELEMENT->VeV;
	Nvve  = ELEMENT->Nvve;

	select_functions(&basis,&grad_basis,&cubature,EType);

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnGc = ELEMENT->NvnGc;
	NvnCs = ELEMENT->NvnCs;
	NvnCc = ELEMENT->NvnCc;
	NvnIs = ELEMENT->NvnIs;
	NvnIc = ELEMENT->NvnIc;
	NvnS  = ELEMENT->NvnS;
	NfnIs = ELEMENT->NfnIs;
	NfnIc = ELEMENT->NfnIc;

	w_vIs = ELEMENT->w_vIs;
	w_vIc = ELEMENT->w_vIc;

	ChiS_vP    = ELEMENT->ChiS_vP;
	ChiS_vIs   = ELEMENT->ChiS_vIs;
	ChiS_vIc   = ELEMENT->ChiS_vIc;
	ChiInvS_vS = ELEMENT->ChiInvS_vS;

	ICs = ELEMENT->ICs;
	ICc = ELEMENT->ICc;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vCs = ELEMENT->I_vGs_vCs;
	I_vGs_vIs = ELEMENT->I_vGs_vIs;
	I_vGs_vIc = ELEMENT->I_vGs_vIc;
	I_vGs_vS  = ELEMENT->I_vGs_vS;
	I_vGc_vP  = ELEMENT->I_vGc_vP;
	I_vGc_vCc = ELEMENT->I_vGc_vCc;
	I_vGc_vIs = ELEMENT->I_vGc_vIs;
	I_vGc_vIc = ELEMENT->I_vGc_vIc;
	I_vGc_vS  = ELEMENT->I_vGc_vS;
	I_vCs_vIs = ELEMENT->I_vCs_vIs;
	I_vCs_vIc = ELEMENT->I_vCs_vIc;
	I_vCc_vIs = ELEMENT->I_vCc_vIs;
	I_vCc_vIc = ELEMENT->I_vCc_vIc;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vIs = ELEMENT->D_vGs_vIs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vIc = ELEMENT->D_vGc_vIc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	ChiS_fIs = ELEMENT->ChiS_fIs;
	ChiS_fIc = ELEMENT->ChiS_fIc;

	I_vGs_fIs = ELEMENT->I_vGs_fIs;
	I_vGs_fIc = ELEMENT->I_vGs_fIc;
	I_vGc_fIs = ELEMENT->I_vGc_fIs;
	I_vGc_fIc = ELEMENT->I_vGc_fIc;
	I_vCs_fIs = ELEMENT->I_vCs_fIs;
	I_vCs_fIc = ELEMENT->I_vCs_fIc;
	I_vCc_fIs = ELEMENT->I_vCc_fIs;
	I_vCc_fIc = ELEMENT->I_vCc_fIc;

	Is_Weak_VV = ELEMENT->Is_Weak_VV;
	Ic_Weak_VV = ELEMENT->Ic_Weak_VV;
	Is_Weak_FF = ELEMENT->Is_Weak_FF;
	Ic_Weak_FF = ELEMENT->Ic_Weak_FF;
	Ds_Weak_VV = ELEMENT->Ds_Weak_VV;
	Dc_Weak_VV = ELEMENT->Dc_Weak_VV;

	// Allocate memory for arrays with multiple levels of dereferencing
	BCoords_V     = malloc(           sizeof *BCoords_V); // free

	rst_vIs       = malloc(NVREFMAX * sizeof *rst_vIs); // free
	rst_vIc       = malloc(NVREFMAX * sizeof *rst_vIc); // free

	w_fIs         = calloc(NP , sizeof *w_fIs); // free
	w_fIc         = calloc(NP , sizeof *w_fIc); // free
	for (P = 0; P < NP; P++) {
		w_fIs[P] = calloc(NESUBCMAX , sizeof **w_fIs); // free
		w_fIc[P] = calloc(NESUBCMAX , sizeof **w_fIc); // free
	}

	GradChiGs_vCs = malloc(dE * sizeof *GradChiGs_vCs); // free
	GradChiGs_vIs = malloc(dE * sizeof *GradChiGs_vIs); // free
	GradChiGc_vCc = malloc(dE * sizeof *GradChiGc_vCc); // free
	GradChiGc_vIc = malloc(dE * sizeof *GradChiGc_vIc); // free
	GradChiCs_vCs = malloc(dE * sizeof *GradChiCs_vCs); // free
	GradChiCc_vCc = malloc(dE * sizeof *GradChiCc_vCc); // free
	GradChiS_vIs  = malloc(dE * sizeof *GradChiS_vIs);  // free
	GradChiS_vIc  = malloc(dE * sizeof *GradChiS_vIc);  // free

	// VOLUME Nodes (Order Independent)
	plotting_element_info(&rst_vP,&dummyPtr_ui[0],&dummyPtr_ui[1],&NvnP,&dummy_ui,PP,EType); // free
	free(dummyPtr_ui[0]);
	free(dummyPtr_ui[1]);

	cubature(&rst_vGs,&dummyPtr_d,&dummyPtr_ui[0],&NvnGs[0],&dummy_ui,0,PGs,dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
	// Use E_rst_vC instead of rst_vGs (!= for PYR ELEMENTs)
	free(rst_vGs);
	E_rst_vC = get_rst_vC(ELEMENT);
	rst_vC = malloc(Nve*dE * sizeof *rst_vC); // free


	// Preliminary Operators
	IGs = identity_d(NvnGs[0]); // free

	ChiRefGs_vGs = basis(PGs,E_rst_vC,NvnGs[0],&Nbf,dE); // free

	if (strstr(BasisType,"Modal") != NULL) {
		ChiGs_vGs = ChiRefGs_vGs;
	} else if (strstr(BasisType,"Nodal") != NULL) {
		ChiGs_vGs = IGs;
	}

	ChiRefInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiRefGs_vGs,IGs); // free

	ChiInvGs_vGs = inverse_d(NvnGs[0],NvnGs[0],ChiGs_vGs,IGs); // free

	TGs = mm_Alloc_d(CBRM,CBNT,CBNT,
	                 NvnGs[0],NvnGs[0],NvnGs[0],1.0,ChiRefInvGs_vGs,ChiGs_vGs); // free

	ChiRefGs_vP = basis(PGs,rst_vP,NvnP,&Nbf,dE); // free

	ChiGs_vP = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,NvnGs[0],NvnGs[0],1.0,ChiRefGs_vP,TGs); // free

	// Returned Operators
	I_vGs_vP[0][PP][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,NvnGs[0],NvnGs[0],1.0,ChiGs_vP,ChiInvGs_vGs); // keep

	free(IGs);
	free(ChiRefGs_vGs);
	free(ChiRefGs_vP);
	free(ChiGs_vP);

	// Get Barycentric coordinates for lower dimensional ELEMENTs
	NFTypes = 1;
	BCoords_dEm1[0] = get_BCoords_dEm1(ELEMENT,0); // keep/free
	if (EType == WEDGE || EType == PYR) {
		NFTypes = 2;
		BCoords_dEm1[1] = get_BCoords_dEm1(ELEMENT,1);
	}

	for (IndFType = 0; IndFType < NFTypes; IndFType++)
		BCoords_F[IndFType] = malloc(sizeof *BCoords_F[IndFType]); // tbd


	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		B_Nve[IndFType]         = BCoords_dEm1[IndFType]->Nve;
		BCoords_F[IndFType]->Is = BCoords_dEm1[IndFType]->BCoords_Is;
		BCoords_F[IndFType]->Ic = BCoords_dEm1[IndFType]->BCoords_Ic;

		// Store Nfn* in a manner consistent with Nvn*
		for (P = 0; P <= PMax; P++) {
			NfnIs[P][IndFType] = BCoords_dEm1[IndFType]->NfnIs[P];
			NfnIc[P][IndFType] = BCoords_dEm1[IndFType]->NfnIc[P];
			w_fIs[P][IndFType] = BCoords_dEm1[IndFType]->w_fIs[P];
			w_fIc[P][IndFType] = BCoords_dEm1[IndFType]->w_fIc[P];
		}
		free(BCoords_dEm1[IndFType]->NfnIs);
		free(BCoords_dEm1[IndFType]->NfnIc);
		free(BCoords_dEm1[IndFType]->w_fIs);
		free(BCoords_dEm1[IndFType]->w_fIc);
	}

	// NvrefSF used to build higher-dimensional operators.
	switch (Adapt) {
	case ADAPT_0:
		PSMin = PGlobal;
		PSMax = PGlobal;
		Nvref = NvrefSF = 1;
		break;
	case ADAPT_P:
		PSMin = 0;
		PSMax = PMax;
		Nvref = NvrefSF = 1;
		break;
	case ADAPT_H:
	default: // ADAPT_HP
		if (Adapt == ADAPT_HP) {
			PSMin = 0;
			PSMax = PMax;
		} else if (Adapt == ADAPT_H) {
			PSMin = PGlobal;
			PSMax = PGlobal;
		}
		switch (EType) {
		case LINE:
			Nvref   = NvrefSF = NREFMAXLINE;
			break;
		case TRI:
			Nvref   = NvrefSF = NREFMAXTRI;
			break;
		case TET:
			Nvref   = NREFMAXTET;
			NvrefSF = 1;
			break;
		case PYR:
			Nvref   = NREFMAXPYR;
			NvrefSF = 1;
			break;
		default:
			printf("Error: Unsupported EType in get_NvrefSF.\n"), exit(1);
			break;
		}
		break;
	}

	ones_Nf = malloc(Nf * sizeof *ones_Nf); // free
	for (f = 0; f < Nf; f++)
		ones_Nf[f] = 1;

	switch (Adapt) {
		default: // ADAPT_H or ADAPT_HP
			fhMax = Nfref;
			break;
		case ADAPT_0:
		case ADAPT_P:
			fhMax = ones_Nf;
			break;
	}

	BCoords_V->Is = calloc(NP , sizeof *(BCoords_V->Is)); // free
	BCoords_V->Ic = calloc(NP , sizeof *(BCoords_V->Ic)); // free

	for (P = PSMin; P <= PSMax; P++) {
		cubature(&rst_vGc,&dummyPtr_d,&dummyPtr_ui[0],&NvnGc[P],&dummy_ui,0,PGc[P],         dE,NodeTypeG[Eclass]  );    free(dummyPtr_ui[0]); // free
		cubature(&rst_vCs,&dummyPtr_d,&dummyPtr_ui[0],&NvnCs[P],&dummy_ui,0,PCs[P][Eclass], dE,NodeTypeG[Eclass]  );    free(dummyPtr_ui[0]); // free
		cubature(&rst_vCc,&dummyPtr_d,&dummyPtr_ui[0],&NvnCc[P],&dummy_ui,0,PCc[P][Eclass], dE,NodeTypeG[Eclass]  );    free(dummyPtr_ui[0]); // free
		cubature(&rst_vS, &dummyPtr_d,&dummyPtr_ui[0],&NvnS[P], &dummy_ui,0,P,              dE,NodeTypeS[P][Eclass]);   free(dummyPtr_ui[0]); // free

		// Preliminary Operators
		IGc          = identity_d(NvnGc[P]); // free
		ICs[P][P][0] = identity_d(NvnCs[P]); // keep
		ICc[P][P][0] = identity_d(NvnCc[P]); // keep
		IS           = identity_d(NvnS[P]);  // free

		ChiRefGc_vGc = basis(PGc[P],rst_vGc,NvnGc[P],&Nbf,dE);         // free
		ChiRefCs_vCs = basis(PCs[P][Eclass],rst_vCs,NvnCs[P],&Nbf,dE); // free
		ChiRefCc_vCc = basis(PCc[P][Eclass],rst_vCc,NvnCc[P],&Nbf,dE); // free
		ChiRefS_vS   = basis(P,             rst_vS, NvnS[P], &Nbf,dE); // free

		if (strstr(BasisType,"Modal") != NULL) {
			ChiGc_vGc = ChiRefGc_vGc;
			ChiCs_vCs = ChiRefCs_vCs;
			ChiCc_vCc = ChiRefCc_vCc;
			ChiS_vS   = ChiRefS_vS;
		} else if (strstr(BasisType,"Nodal") != NULL) {
			ChiGc_vGc = IGc;
			ChiCs_vCs = ICs[P][P][0];
			ChiCc_vCc = ICc[P][P][0];
			ChiS_vS   = IS;
		}

		ChiRefInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiRefGc_vGc,IGc);    // free
		ChiRefInvCs_vCs = inverse_d(NvnCs[P],NvnCs[P],ChiRefCs_vCs,ICs[P][P][0]); // free
		ChiRefInvCc_vCc = inverse_d(NvnCc[P],NvnCc[P],ChiRefCc_vCc,ICc[P][P][0]); // free
		ChiRefInvS_vS   = inverse_d(NvnS[P], NvnS[P], ChiRefS_vS,IS);       // free

		ChiInvGc_vGc        = inverse_d(NvnGc[P],NvnGc[P],ChiGc_vGc,IGc);    // free
		ChiInvCs_vCs        = inverse_d(NvnCs[P],NvnCs[P],ChiCs_vCs,ICs[P][P][0]); // free
		ChiInvCc_vCc        = inverse_d(NvnCc[P],NvnCc[P],ChiCc_vCc,ICc[P][P][0]); // free
		ChiInvS_vS[P][P][0] = inverse_d(NvnS[P], NvnS[P], ChiS_vS,  IS);     // keep

		TGc = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGc[P],NvnGc[P],1.0,ChiRefInvGc_vGc,ChiGc_vGc); // free
		TCs = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,ChiRefInvCs_vCs,ChiCs_vCs); // free
		TCc = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,ChiRefInvCc_vCc,ChiCc_vCc); // free
		TS  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P], NvnS[P], NvnS[P], 1.0,ChiRefInvS_vS,ChiS_vS);     // free

		switch (Adapt) {
			default: // ADAPT_P or ADAPT_HP
				if      (P == 0)    PbMin = P,   PbMax = P+1;
				else if (P == PMax) PbMin = P-1, PbMax = PMax;
				else                PbMin = P-1, PbMax = P+1;
				break;
			case ADAPT_0:
			case ADAPT_H:
				PbMin = P; PbMax = P;
				break;
		}
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			// VOLUME Operators
			cubature(&rst_vIs[0],&w_vIs[Pb],&dummyPtr_ui[0],&NvnIs[Pb],&dummy_ui,1,PIvs[Pb][Eclass],dE,NodeTypeIvs[Pb][Eclass]); free(dummyPtr_ui[0]); // free
			cubature(&rst_vIc[0],&w_vIc[Pb],&dummyPtr_ui[0],&NvnIc[Pb],&dummy_ui,1,PIvc[Pb][Eclass],dE,NodeTypeIvc[Pb][Eclass]); free(dummyPtr_ui[0]); // free

			diag_w_vIs = diag_d(w_vIs[Pb],NvnIs[Pb]); // free
			diag_w_vIc = diag_d(w_vIc[Pb],NvnIc[Pb]); // free

			ChiRefGs_vIs = basis(PGs,rst_vIs[0],NvnIs[Pb],&Nbf,dE); // free
			ChiRefGs_vIc = basis(PGs,rst_vIc[0],NvnIc[Pb],&Nbf,dE); // free

			// VOLUME Operators

			BCoords_V->Is[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIs[Pb],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vIs,ChiRefInvGs_vGs); // free
			BCoords_V->Ic[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIc[Pb],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vIc,ChiRefInvGs_vGs); // free

			for (vrefSF = 0; vrefSF < NvrefSF; vrefSF++) {
				if (vrefSF) {
					mm_CTN_d(Nvve[vrefSF],dE,Nve,VeV[vrefSF],E_rst_vC,rst_vC);
					rst_vIs[vrefSF] = mm_Alloc_d(CBCM,CBNT,CBNT,NvnIs[P],dE,Nvve[vrefSF],1.0,BCoords_V->Is[Pb],rst_vC); // free
					rst_vIc[vrefSF] = mm_Alloc_d(CBCM,CBNT,CBNT,NvnIc[P],dE,Nvve[vrefSF],1.0,BCoords_V->Ic[Pb],rst_vC); // free
					free(rst_vC);

					ChiRefGs_vIs = basis(PGs,rst_vIs[vrefSF],NvnIs[Pb],&Nbf,dE); // free
					ChiRefGs_vIc = basis(PGs,rst_vIc[vrefSF],NvnIc[Pb],&Nbf,dE); // free
				}

				ChiRefGc_vIs = basis(PGc[P],        rst_vIs[vrefSF],NvnIs[P],&Nbf,dE); // free
				ChiRefGc_vIc = basis(PGc[P],        rst_vIc[vrefSF],NvnIc[P],&Nbf,dE); // free
				ChiRefCs_vIs = basis(PCs[P][Eclass],rst_vIs[vrefSF],NvnIs[P],&Nbf,dE); // free
				ChiRefCs_vIc = basis(PCs[P][Eclass],rst_vIc[vrefSF],NvnIc[P],&Nbf,dE); // free
				ChiRefCc_vIs = basis(PCc[P][Eclass],rst_vIs[vrefSF],NvnIs[P],&Nbf,dE); // free
				ChiRefCc_vIc = basis(PCc[P][Eclass],rst_vIc[vrefSF],NvnIc[P],&Nbf,dE); // free
				ChiRefS_vIs  = basis(P,             rst_vIs[vrefSF],NvnIs[P],&Nbf,dE); // free
				ChiRefS_vIc  = basis(P,             rst_vIc[vrefSF],NvnIc[P],&Nbf,dE); // free

				ChiGs_vIs               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vIs,TGs); // free
				ChiGs_vIc               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vIc,TGs); // free
				ChiGc_vIs               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vIs,TGc); // free
				ChiGc_vIc               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vIc,TGc); // free
				ChiCs_vIs               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIs,TCs); // free
				ChiCs_vIc               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIc,TCs); // free
				ChiCc_vIs               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIs,TCc); // free
				ChiCc_vIc               = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIc,TCc); // free

				// Returned Operators
				ChiS_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIs,TS);    // keep
				ChiS_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIc,TS);    // keep

				I_vGs_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vIs,ChiInvGs_vGs); // keep
				I_vGs_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vIc,ChiInvGs_vGs); // keep
				I_vGc_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vIs,ChiInvGc_vGc); // keep
				I_vGc_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vIc,ChiInvGc_vGc); // keep
				I_vCs_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnCs[P],NvnCs[P],1.0,ChiCs_vIs,ChiInvCs_vCs); // keep
				I_vCs_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnCs[P],NvnCs[P],1.0,ChiCs_vIc,ChiInvCs_vCs); // keep
				I_vCc_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnCc[P],NvnCc[P],1.0,ChiCc_vIs,ChiInvCc_vCc); // keep
				I_vCc_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnCc[P],NvnCc[P],1.0,ChiCc_vIc,ChiInvCc_vCc); // keep

				free(ChiRefGs_vIs);
				free(ChiRefGs_vIc);
				free(ChiRefGc_vIs);
				free(ChiRefGc_vIc);
				free(ChiRefCs_vIs);
				free(ChiRefCs_vIc);
				free(ChiRefCc_vIs);
				free(ChiRefCc_vIc);

				free(ChiGs_vIs);
				free(ChiGs_vIc);
				free(ChiGc_vIs);
				free(ChiGc_vIc);
				free(ChiCs_vIs);
				free(ChiCs_vIc);
				free(ChiCc_vIs);
				free(ChiCc_vIc);

				free(ChiRefS_vIs);
				free(ChiRefS_vIc);
			}

			// Operators (P -> P only)
			if (P == Pb) {
				ChiRefGs_vGc = basis(PGs,           rst_vGc,   NvnGc[P],&Nbf,dE); // free
				ChiRefGs_vCs = basis(PGs,           rst_vCs,   NvnCs[P],&Nbf,dE); // free
				ChiRefGs_vS  = basis(PGs,           rst_vS,    NvnS[P], &Nbf,dE); // free
				ChiRefGc_vP  = basis(PGc[P],        rst_vP,    NvnP,    &Nbf,dE); // free
				ChiRefGc_vCc = basis(PGc[P],        rst_vCc,   NvnCc[P],&Nbf,dE); // free
				ChiRefGc_vS  = basis(PGc[P],        rst_vS,    NvnS[P], &Nbf,dE); // free
				ChiRefS_vP   = basis(P,             rst_vP,    NvnP,    &Nbf,dE); // free
				ChiRefS_vIs  = basis(P,             rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				ChiRefS_vIc  = basis(P,             rst_vIc[0],NvnIc[P],&Nbf,dE); // free

				ChiGs_vGc         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vGc,TGs); // free
				ChiGs_vCs         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[0],NvnGs[0],1.0,ChiRefGs_vCs,TGs); // free
				ChiGs_vS          = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P], NvnGs[0],NvnGs[0],1.0,ChiRefGs_vS, TGs); // free
				ChiGc_vP          = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiRefGc_vP, TGc); // free
				ChiGc_vCc         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vCc,TGc); // free
				ChiGc_vS          = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P], NvnGc[P],NvnGc[P],1.0,ChiRefGc_vS, TGc); // free
				ChiS_vP[P][PP][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnS[P], NvnS[P], 1.0,ChiRefS_vP  ,TS);  // keep

				GradChiRefGs_vCs = grad_basis(PGs,           rst_vCs,   NvnCs[P],&Nbf,dE); // free
				GradChiRefGs_vIs = grad_basis(PGs,           rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				GradChiRefGc_vCc = grad_basis(PGc[P],        rst_vCc,   NvnCc[P],&Nbf,dE); // free
				GradChiRefGc_vIc = grad_basis(PGc[P],        rst_vIc[0],NvnIc[P],&Nbf,dE); // free
				GradChiRefCs_vCs = grad_basis(PCs[P][Eclass],rst_vCs,   NvnCs[P],&Nbf,dE); // free
				GradChiRefCc_vCc = grad_basis(PCc[P][Eclass],rst_vCc,   NvnCc[P],&Nbf,dE); // free
				GradChiRefS_vIs  = grad_basis(P,             rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				GradChiRefS_vIc  = grad_basis(P,             rst_vIc[0],NvnIc[P],&Nbf,dE); // free

				for (dim = 0; dim < dE; dim++) {
					GradChiGs_vCs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[0],NvnGs[0],1.0,GradChiRefGs_vCs[dim],TGs); // free
					GradChiGs_vIs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[0],NvnGs[0],1.0,GradChiRefGs_vIs[dim],TGs); // free
					GradChiGc_vCc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vCc[dim],TGc); // free
					GradChiGc_vIc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vIc[dim],TGc); // free
					GradChiCs_vCs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiRefCs_vCs[dim],TCs); // free
					GradChiCc_vCc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiRefCc_vCc[dim],TCc); // free
					GradChiS_vIs[dim]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnS[P], NvnS[P], 1.0,GradChiRefS_vIs[dim],TS);   // free
					GradChiS_vIc[dim]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnS[P], NvnS[P], 1.0,GradChiRefS_vIc[dim],TS);   // free
				}

				// Returned Operators
				// VOLUME related operators
				I_vGs_vGc[P][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vGc,ChiInvGs_vGs); // keep
				I_vGs_vCs[P][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[0],NvnGs[0],1.0,ChiGs_vCs,ChiInvGs_vGs); // keep
				I_vGs_vS[P][Pb][0]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P] ,NvnGs[0],NvnGs[0],1.0,ChiGs_vS, ChiInvGs_vGs); // keep
				I_vGc_vP[P][PP][0]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiGc_vP, ChiInvGc_vGc); // keep
				I_vGc_vCc[P][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vCc,ChiInvGc_vGc); // keep
				I_vGc_vS[P][Pb][0]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P], NvnGc[P],NvnGc[P],1.0,ChiGc_vS, ChiInvGc_vGc); // keep

				if (EFE) {
					Is_Weak_VV[P][Pb][0] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIs[P],NvnIs[P],1.0,ChiS_vIs[P][P][0],diag_w_vIs); // keep
					Ic_Weak_VV[P][Pb][0] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIc[P],NvnIc[P],1.0,ChiS_vIc[P][P][0],diag_w_vIc); // keep
				} else {
					;
				}

				for (dim = 0; dim < dE; dim++) {
					D_vGs_vCs[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[0],NvnGs[0],1.0,GradChiGs_vCs[dim],ChiInvGs_vGs); // keep
					D_vGs_vIs[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[0],NvnGs[0],1.0,GradChiGs_vIs[dim],ChiInvGs_vGs); // keep
					D_vGc_vCc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vCc[dim],ChiInvGc_vGc); // keep
					D_vGc_vIc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vIc[dim],ChiInvGc_vGc); // keep
					D_vCs_vCs[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiCs_vCs[dim],ChiInvCs_vCs); // keep
					D_vCc_vCc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiCc_vCc[dim],ChiInvCc_vCc); // keep

					if (EFE) {
						Ds_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIs[P],NvnIs[P],1.0,GradChiS_vIs[dim],diag_w_vIs); // keep
						Dc_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIc[P],NvnIc[P],1.0,GradChiS_vIc[dim],diag_w_vIc); // keep
					} else {
						;
					}
				}

				free(ChiRefGs_vGc);
				free(ChiRefGs_vCs);
				free(ChiRefGs_vS);
				free(ChiRefGc_vP);
				free(ChiRefGc_vCc);
				free(ChiRefGc_vS);
				free(ChiRefS_vP);
				free(ChiRefS_vIs);
				free(ChiRefS_vIc);

				free(ChiGs_vGc);
				free(ChiGs_vCs);
				free(ChiGs_vS);
				free(ChiGc_vP);
				free(ChiGc_vCc);
				free(ChiGc_vS);

				array_free2_d(dE,GradChiRefGs_vCs);
				array_free2_d(dE,GradChiRefGs_vIs);
				array_free2_d(dE,GradChiRefGc_vCc);
				array_free2_d(dE,GradChiRefGc_vIc);
				array_free2_d(dE,GradChiRefCs_vCs);
				array_free2_d(dE,GradChiRefCc_vCc);
				array_free2_d(dE,GradChiRefS_vIs);
				array_free2_d(dE,GradChiRefS_vIc);

				for (dim = 0; dim < dE; dim++) {
					free(GradChiGs_vCs[dim]);
					free(GradChiGs_vIs[dim]);
					free(GradChiGc_vCc[dim]);
					free(GradChiGc_vIc[dim]);
					free(GradChiCs_vCs[dim]);
					free(GradChiCc_vCc[dim]);
					free(GradChiS_vIs[dim]);
					free(GradChiS_vIc[dim]);
				}
			}

			// FACET related operators
//Is_Weak_FV[P][Pb] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIs[Pb],NvnIs[Pb],ChiS_vIs[P][Pb],diag_w_vIs); // keep
//Ic_Weak_FV[P][Pb] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIc[Pb],NvnIc[Pb],ChiS_vIc[P][Pb],diag_w_vIc); // keep

			for (f = 0; f < Nf; f++) {
				IndFType = get_IndFType(Eclass,f);

				for (fh = 0; fh < fhMax[f]; fh++) {
					Vf = f*NFREFMAX+fh;

					mm_CTN_d(Nfve[f],dE,Nve,VeF[Vf],E_rst_vC,rst_vC);
//array_print_d(Nfve[IndFType],dE,rst_vC,'C');

					rst_fIs = mm_Alloc_d(CBCM,CBNT,CBNT,NfnIs[Pb][IndFType],dE,B_Nve[IndFType],1.0,BCoords_F[IndFType]->Is[Pb],rst_vC); // free
					rst_fIc = mm_Alloc_d(CBCM,CBNT,CBNT,NfnIc[Pb][IndFType],dE,B_Nve[IndFType],1.0,BCoords_F[IndFType]->Ic[Pb],rst_vC); // free

					diag_w_fIs = diag_d(w_fIs[Pb][IndFType],NfnIs[Pb][IndFType]); // free
					diag_w_fIc = diag_d(w_fIc[Pb][IndFType],NfnIc[Pb][IndFType]); // free

					ChiRefS_fIs  = basis(P,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefS_fIc  = basis(P,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					// Returned Operators
					ChiS_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnS[P],NvnS[P],1.0,ChiRefS_fIs,TS); // keep
					ChiS_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnS[P],NvnS[P],1.0,ChiRefS_fIc,TS); // keep

					Is_Weak_FF[P][Pb][Vf] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NfnIs[Pb][IndFType],NfnIs[Pb][IndFType],-1.0,ChiS_fIs[P][Pb][Vf],diag_w_fIs); // keep
					Ic_Weak_FF[P][Pb][Vf] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NfnIc[Pb][IndFType],NfnIc[Pb][IndFType],-1.0,ChiS_fIc[P][Pb][Vf],diag_w_fIc); // keep

					if (fh == 0) {
						ChiRefGs_fIs = basis(PGs           ,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
						ChiRefGs_fIc = basis(PGs           ,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
						ChiRefGc_fIs = basis(PGc[P]        ,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
						ChiRefGc_fIc = basis(PGc[P]        ,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
						ChiRefCs_fIs = basis(PCs[P][Eclass],rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
						ChiRefCs_fIc = basis(PCs[P][Eclass],rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
						ChiRefCc_fIs = basis(PCc[P][Eclass],rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
						ChiRefCc_fIc = basis(PCc[P][Eclass],rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

						ChiGs_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiRefGs_fIs,TGs); // free
						ChiGs_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiRefGs_fIc,TGs); // free
						ChiGc_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIs,TGc); // free
						ChiGc_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIc,TGc); // free
						ChiCs_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIs,TCs); // free
						ChiCs_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIc,TCs); // free
						ChiCc_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIs,TCc); // free
						ChiCc_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIc,TCc); // free

						// Returned Operators
						I_vGs_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiGs_fIs,ChiInvGs_vGs); // keep
						I_vGs_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[0],NvnGs[0],1.0,ChiGs_fIc,ChiInvGs_vGs); // keep
						I_vGc_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIs,ChiInvGc_vGc); // keep
						I_vGc_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIc,ChiInvGc_vGc); // keep
						I_vCs_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIs,ChiInvCs_vCs); // keep
						I_vCs_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIc,ChiInvCs_vCs); // keep
						I_vCc_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIs,ChiInvCc_vCc); // keep
						I_vCc_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIc,ChiInvCc_vCc); // keep

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

					free(rst_fIs);
					free(rst_fIc);

					free(diag_w_fIs);
					free(diag_w_fIc);

					free(ChiRefS_fIs);
					free(ChiRefS_fIc);
				}
			}
			free(diag_w_vIs);
			free(diag_w_vIc);
		}

		free(rst_vGc);
		free(rst_vCs);
		free(rst_vCc);
		free(rst_vS);

		for (vref = 0; vref < Nvref; vref++) {
			free(rst_vIs[vref]);
			free(rst_vIc[vref]);
		}
		free(rst_vIs);
		free(rst_vIc);

		free(IGc);
		free(IS);

		free(ChiRefGc_vGc);
		free(ChiRefCs_vCs);
		free(ChiRefCc_vCc);
		free(ChiRefS_vS);

		free(ChiRefInvGc_vGc);
		free(ChiRefInvCs_vCs);
		free(ChiRefInvCc_vCc);
		free(ChiRefInvS_vS);

		free(ChiInvGc_vGc);
		free(ChiInvCs_vCs);
		free(ChiInvCc_vCc);

		free(TGc);
		free(TCs);
		free(TCc);
		free(TS);
	}

	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		array_free2_d(PMax+1,BCoords_F[IndFType]->Is);
		array_free2_d(PMax+1,BCoords_F[IndFType]->Ic);
		free(BCoords_F[IndFType]);
		free(BCoords_dEm1[IndFType]);
	}

	array_free3_d(NP,NESUBCMAX,w_fIs);
	array_free3_d(NP,NESUBCMAX,w_fIc);

	array_free2_d(NP,BCoords_V->Is);
	array_free2_d(NP,BCoords_V->Ic);
	free(BCoords_V);

	free(E_rst_vC);
	free(rst_vC);
	free(rst_vP);

	free(ChiRefInvGs_vGs);
	free(ChiInvGs_vGs);
	free(TGs);

	free(GradChiGs_vCs);
	free(GradChiGs_vIs);
	free(GradChiGc_vCc);
	free(GradChiGc_vIc);
	free(GradChiCs_vCs);
	free(GradChiCc_vCc);
	free(GradChiS_vIs);
	free(GradChiS_vIc);

	free(ones_Nf);
}

void setup_ELEMENT_VeF(const unsigned int EType)
{
	/*
	 *	Comments:
	 *		While VeV defines VeVref_LINE and VeVref_TRI, the redundant definitions are included here as the ELEMENTs
	 *		are set up in order and VeV may not yet be defined.
	 */

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
			ELEMENT->Nfref[f] = NREFMAXPOINT;
		break;
	case TRI:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = NREFMAXLINE;
		break;
	case QUAD:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = NREFMAXLINE;
		break;
	case TET:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = NREFMAXTRI;
		break;
	case HEX:
		for (f = 0; f < Nf; f++)
			ELEMENT->Nfref[f] = NREFMAXQUAD;
		break;
	case WEDGE:
		for (f = 0; f < 3; f++)
			ELEMENT->Nfref[f] = NREFMAXQUAD;
		for (f = 3; f < 5; f++)
			ELEMENT->Nfref[f] = NREFMAXTRI;
		break;
	case PYR:
		for (f = 0; f < 4; f++)
			ELEMENT->Nfref[f] = NREFMAXTRI;
		for (f = 4; f < 5; f++)
			ELEMENT->Nfref[f] = NREFMAXQUAD;
		break;
	default:
		printf("Error: Unsupported EType in setup_ELEMENT_VeF.\n"), exit(1);
		break;
	}

	Nfref   = ELEMENT->Nfref;

	unsigned int size_VeF = 0;
	double VeVref_LINE[12]  = {1.0 , 0.0 ,
	                           0.0 , 1.0 ,
	                           1.0 , 0.0 ,
	                           0.5 , 0.5 ,
	                           0.5 , 0.5 ,
	                           0.0 , 1.0 },
	       VeVref_TRI[45]   = {1.0 , 0.0 , 0.0 ,
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
	       VeVref_QUAD[144] = {1.0 , 0.0 , 0.0 , 0.0 ,
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
	case QUAD:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*4+k]] = VeVref_LINE[j*kMax+k];
		}}}

		break;
	case TET:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*4+k]] = VeVref_TRI[j*kMax+k];
		}}}

		break;
	case HEX:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*4+k]] = VeVref_QUAD[j*kMax+k];
		}}}

		break;
	case WEDGE:
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
					VeF[iStep+j*Nve+VeFcon[i*4+k]] = VeVref_TRI[j*kMax+k];
				else if (Nfve[i] == 4)
					VeF[iStep+j*Nve+VeFcon[i*4+k]] = VeVref_QUAD[j*kMax+k];
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
			ELEMENT->VeF[i*NFREFMAX+j] = malloc(Nfve[i]*Nve * sizeof **(ELEMENT->VeF)); // keep
			for (k = 0, kMax = Nfve[i];  k < kMax; k++) {
			for (l = 0, lMax = Nve;      l < lMax; l++) {
				ELEMENT->VeF[i*NFREFMAX+j][k*Nve+l] = VeF[iStep+j*(Nfve[i]*Nve)+k*(Nve)+l];
			}}
		}
	}

/*
printf("\n\nEType: %d\n\n\n",EType);
i = 0;
for (j = 0; j < Nfref[i]; j++)
	array_print_d(Nfve[i],Nve,ELEMENT->VeF[i*NFREFMAX+j],'R');
*/

	free(VeF);
}

static double *sf_assemble_d(const unsigned int NIn[3], const unsigned int NOut[3], const unsigned int d,
                             double *BOP[3])
{
	/*
	 *	Purpose:
	 *		Assemble ST(andard) OP(erators) for TP elements using the lower dimensional operators.
	 *
	 *	Comments:
	 *		Standard operators are assembled using sparse BLAS multiplications of lower dimensional operators.
	 *		Sparse matrices are stored in compressed sparse row (CSR) format.
	 *
	 *	Notation:
	 *		BOP : (B)lock (OP)erator in each TP direction.
	 *
	 *	References:
	 *		Intel MKL Sparse BLAS CSR Matrix Storage Format: https://software.intel.com/en-us/node/599835
	 */

	char         transa, matdescra[6];
	unsigned int dim, i, j, k, Bcol, iMax, jMax, kMax, BcolMax, Gcol,
	             Indd, IndG, IndGrow;
	double *OP_ST;

	unsigned int NNZ_BOP[3] = {NIn[0]*NOut[0], NIn[1]*NOut[1], NIn[2]*NOut[2]};
	MKL_INT      BRows[3] = {NIn[1]*NIn[2], NOut[0]*NIn[2], NOut[0]*NOut[1]},
	             dims_OP_ST[2] = {NOut[0]*NOut[1]*NOut[2], NIn[0]*NIn[1]*NIn[2]},
	             dims_DOPr[2]  = {NOut[0]*NIn[1]*NIn[2],   NIn[0]*NIn[1]*NIn[2]},
	             dims_DOPs[2]  = {NOut[0]*NOut[1]*NIn[2],  NOut[0]*NIn[1]*NIn[2]},
	             dims_DOPt[2]  = {NOut[0]*NOut[1]*NOut[2], NOut[0]*NOut[1]*NIn[2]};

	MKL_INT      OPr_rowIndex[BRows[0]*NOut[0]+1], OPs_rowIndex[BRows[1]*NOut[1]+1], OPt_rowIndex[BRows[2]*NOut[2]+1],
	             OPr_cols[BRows[0]*NNZ_BOP[0]],    OPs_cols[BRows[1]*NNZ_BOP[1]],    OPt_cols[BRows[2]*NNZ_BOP[2]];
	double       OPr_vals[BRows[0]*NNZ_BOP[0]],    OPs_vals[BRows[1]*NNZ_BOP[1]],    OPt_vals[BRows[2]*NNZ_BOP[2]],
	             OPr_ST[dims_DOPr[0]*dims_DOPr[1]], OPInter_ST[dims_DOPs[0]*dims_DOPr[1]],
	             alpha, beta, one_d[1] = {1.0};

	OP_ST = malloc(dims_OP_ST[0]*dims_OP_ST[1] * sizeof *OP_ST); // keep (requires external free)

	if (d == 1)
		printf("Error: d must be greater than 1 in sf_assemble_d.\n"), exit(1);

	for (dim = 0; dim < 3; dim++) {
		if (BOP[dim] == NULL)
			BOP[dim] = one_d;
	}

	// r
	Indd = 0;

	IndG = 0; IndGrow = 0;
	for (k = 0, kMax = NIn[2];  k < kMax; k++) {
	for (j = 0, jMax = NIn[1];  j < jMax; j++) {
	for (i = 0, iMax = NOut[0]; i < iMax; i++) {
		OPr_rowIndex[IndGrow] = IndG;
		IndGrow++;

		for (Bcol = 0, BcolMax = NIn[Indd]; Bcol < BcolMax; Bcol++) {
			Gcol = Bcol + j*BcolMax + k*BcolMax*jMax;

			OPr_cols[IndG] = Gcol;
			OPr_vals[IndG] = BOP[Indd][i*BcolMax+Bcol];

			IndG++;
		}
	}}}
	OPr_rowIndex[IndGrow] = BRows[Indd]*NNZ_BOP[Indd];

	// s
	Indd = 1;

	IndG = 0; IndGrow = 0;
	for (k = 0, kMax = NIn[2];  k < kMax; k++) {
	for (j = 0, jMax = NOut[1]; j < jMax; j++) {
	for (i = 0, iMax = NOut[0]; i < iMax; i++) {
		OPs_rowIndex[IndGrow] = IndG;
		IndGrow++;

		for (Bcol = 0, BcolMax = NIn[Indd]; Bcol < BcolMax; Bcol++) {
			Gcol = i + Bcol*iMax + k*iMax*BcolMax;

			OPs_cols[IndG] = Gcol;
			OPs_vals[IndG] = BOP[Indd][j*BcolMax+Bcol];

			IndG++;
		}
	}}}
	OPs_rowIndex[IndGrow] = BRows[Indd]*NNZ_BOP[Indd];

/*
array_print_i(1,BRows[0]*NOut[0]+1,OPr_rowIndex,'R');
array_print_i(1,BRows[0]*NNZ_BOP[0],OPr_cols,'R');
array_print_d(1,BRows[0]*NNZ_BOP[0],OPr_vals,'R');
*/

	MKL_INT job[8] = {1,0,0,2,BRows[0]*NNZ_BOP[0],1}, info = 0;
	mkl_ddnscsr(job,&dims_DOPr[0],&dims_DOPr[1],OPr_ST,&dims_DOPr[1],OPr_vals,OPr_cols,OPr_rowIndex,&info);
	if (info != 0)
		printf("Problem converting from sparse to dense in mkl_ddnscsr.\n"), exit(1);

	alpha = 1.0;
	beta  = 0.0;
	transa = 'N';
	matdescra[0] = 'G';
	matdescra[1] = 'L'; // not used
	matdescra[2] = 'N'; // not used
	matdescra[3] = 'C';

	if (d == 2) {
		mkl_dcsrmm(&transa,&dims_DOPs[0],&dims_DOPr[1],&dims_DOPs[1],&alpha,matdescra,OPs_vals,OPs_cols,
		           OPs_rowIndex,&OPs_rowIndex[1],OPr_ST,&dims_DOPr[1],&beta,OP_ST,&dims_DOPr[1]);
	} else {
		mkl_dcsrmm(&transa,&dims_DOPs[0],&dims_DOPr[1],&dims_DOPs[1],&alpha,matdescra,OPs_vals,OPs_cols,
		           OPs_rowIndex,&OPs_rowIndex[1],OPr_ST,&dims_DOPr[1],&beta,OPInter_ST,&dims_DOPr[1]);

		Indd = 2;
		IndG = 0; IndGrow = 0;
		for (k = 0, kMax = NOut[2]; k < kMax; k++) {
		for (j = 0, jMax = NOut[1]; j < jMax; j++) {
		for (i = 0, iMax = NOut[0]; i < iMax; i++) {
			OPt_rowIndex[IndGrow] = IndG;
			IndGrow++;

			for (Bcol = 0, BcolMax = NIn[Indd]; Bcol < BcolMax; Bcol++) {
				Gcol = i + j*iMax + Bcol*iMax*jMax;

				OPt_cols[IndG] = Gcol;
				OPt_vals[IndG] = BOP[Indd][k*BcolMax+Bcol];

				IndG++;
			}
		}}}
		OPt_rowIndex[IndGrow] = BRows[Indd]*NNZ_BOP[Indd];

		mkl_dcsrmm(&transa,&dims_DOPt[0],&dims_DOPr[1],&dims_DOPt[1],&alpha,matdescra,OPt_vals,OPt_cols,
		           OPt_rowIndex,&OPt_rowIndex[1],OPInter_ST,&dims_DOPr[1],&beta,OP_ST,&dims_DOPr[1]);
	}

//array_print_d(dims_OP_ST[0],dims_OP_ST[1],OP_ST,'R');
//exit(1);

	return OP_ST;
}

void get_sf_parameters(const unsigned int NIn0, const unsigned int NOut0, double *OP0,
                       const unsigned int NIn1, const unsigned int NOut1, double *OP1,
                       unsigned int NIn_SF[3], unsigned int NOut_SF[3], double *OP_SF[3],
                       const unsigned int d, const unsigned int dim1, const unsigned int Eclass)
{
	/*
	 *	Purpose:
	 *		Set up (s)um (f)actorization parameters.
	 *
	 *	Comments:
	 *		Clearly the input parameters are not very elegant, however, this results in the least overhead when calling
	 *		this function (as compared to using a structure for example).
	 */

	unsigned int dim;

	if (Eclass == C_TP) {
		for (dim = 0; dim < 3; dim++) {
			if (dim == dim1) {
				NIn_SF[dim]  = NIn1;
				NOut_SF[dim] = NOut1;
				OP_SF[dim]   = OP1;
			} else if (dim < d) {
				NIn_SF[dim]  = NIn0;
				NOut_SF[dim] = NOut0;
				OP_SF[dim]   = OP0;
			} else {
				NIn_SF[dim]  = 1;
				NOut_SF[dim] = 1;
				OP_SF[dim]   = NULL;
			}
		}
	} else if (Eclass == C_WEDGE) {
		for (dim = 0; dim < 3; dim++) {
			if (dim == 0) {
				NIn_SF[dim]  = NIn0;
				NOut_SF[dim] = NOut0;
				OP_SF[dim]   = OP0;
			} else if (dim == 2) {
				NIn_SF[dim]  = NIn1;
				NOut_SF[dim] = NOut1;
				OP_SF[dim]   = OP1;
			} else {
				NIn_SF[dim]  = 1;
				NOut_SF[dim] = 1;
				OP_SF[dim]   = NULL;
			}
		}
	} else {
		printf("Error: Unsupported Eclass in get_sf_parameters.\n"), exit(1);
	}
}

static void get_sf_parametersF(const unsigned int NIn0, const unsigned int NOut0, double **OP0,
                               const unsigned int NIn1, const unsigned int NOut1, double **OP1,
                               unsigned int NIn_SF[3], unsigned int NOut_SF[3], double *OP_SF[3],
                               const unsigned int d, const unsigned int Vf, const unsigned int Eclass)
{
	/*
	 *	Purpose:
	 *		Set up (s)um (f)actorization parameters for (F)ACET operators.
	 *
	 *	Comments:
	 *		While get_sf_parameters is sufficient for QUADs, HEX ELEMENTs may require three different FACET operators if
	 *		h-refinement is employed.
	 */

	unsigned int f, fh, dimF, dimV1, dimV2;

	f  = Vf/NFREFMAX;
	fh = Vf % NFREFMAX;

	dimF = f/2;

	switch (d) {
	default: // 3D
		switch (Eclass) {
		default: // C_TP
			// FACET term (standard treatment)
			NIn_SF[dimF]  = NIn1;
			NOut_SF[dimF] = NOut1;
			OP_SF[dimF]   = OP1[(f%2)*NFREFMAX];

			// VOLUME term NIn/NOut (standard treatment)
			dimV1 = (dimF+1)%3;
			dimV2 = (dimV1+1)%3;

			NIn_SF[dimV1]  = NIn0;
			NOut_SF[dimV1] = NOut0;
			NIn_SF[dimV2]  = NIn0;
			NOut_SF[dimV2] = NOut0;

			// VOLUME term OP_SF
			switch (fh) {
			case 0: OP_SF[dimV1] = OP0[0]; OP_SF[dimV2] = OP0[0]; break; // Conforming
			case 1: OP_SF[dimV1] = OP0[1]; OP_SF[dimV2] = OP0[0]; break;
			case 2: OP_SF[dimV1] = OP0[2]; OP_SF[dimV2] = OP0[0]; break;
			case 3: OP_SF[dimV1] = OP0[0]; OP_SF[dimV2] = OP0[1]; break;
			case 4: OP_SF[dimV1] = OP0[0]; OP_SF[dimV2] = OP0[2]; break;
			case 5: OP_SF[dimV1] = OP0[1]; OP_SF[dimV2] = OP0[1]; break;
			case 6: OP_SF[dimV1] = OP0[2]; OP_SF[dimV2] = OP0[1]; break;
			case 7: OP_SF[dimV1] = OP0[1]; OP_SF[dimV2] = OP0[2]; break;
			case 8: OP_SF[dimV1] = OP0[2]; OP_SF[dimV2] = OP0[2]; break;
			default:
				printf("Error: Unsupported fh in get_sf_parametersF.\n"), exit(1);
				break;
			}
			break;
		case C_WEDGE:
			// OP_SF[1] (standard)
			NIn_SF[1]  = 1;
			NOut_SF[1] = 1;
			OP_SF[1]   = NULL;

			// dim = 0 and dim = 2
			NIn_SF[0]  = NIn0;
			NOut_SF[0] = NOut0;
			NIn_SF[2]  = NIn1;
			NOut_SF[2] = NOut1;

			if (f < 3) { // QUAD FACETs
				switch (fh) {
				case 0: OP_SF[0] = OP0[0]; OP_SF[2] = OP1[0]; break;
				case 1: OP_SF[0] = OP0[1]; OP_SF[2] = OP1[0]; break;
				case 2: OP_SF[0] = OP0[2]; OP_SF[2] = OP1[0]; break;
				case 3: OP_SF[0] = OP0[0]; OP_SF[2] = OP1[1]; break;
				case 4: OP_SF[0] = OP0[0]; OP_SF[2] = OP1[2]; break;
				case 5: OP_SF[0] = OP0[1]; OP_SF[2] = OP1[1]; break;
				case 6: OP_SF[0] = OP0[2]; OP_SF[2] = OP1[1]; break;
				case 7: OP_SF[0] = OP0[1]; OP_SF[2] = OP1[2]; break;
				case 8: OP_SF[0] = OP0[2]; OP_SF[2] = OP1[2]; break;
				default:
					printf("Error: Unsupported fh for f < 3 in get_sf_parametersF (C_WEDGE).\n"), exit(1);
					break;
				}
			} else if (f < 5) { // TRI FACETs
				if (fh > 4)
					printf("Error: Unsupported fh for f < 5 in get_sf_parametersF (C_WEDGE).\n"), exit(1);

				OP_SF[0] = OP0[fh];
				OP_SF[2] = OP1[((f+1)%2)*NFREFMAX];
			} else {
				printf("Error: Unsupported f in get_sf_parametersF (C_WEDGE).\n"), exit(1);
			}
			break;
		}
		break;
	case 2:
		get_sf_parameters(NIn0,NOut0,OP0[fh],NIn1,NOut1,OP1[(f%2)*NFREFMAX],NIn_SF,NOut_SF,OP_SF,d,dimF,C_TP);
		break;
	}
}

static void setup_TP_operators(const unsigned int EType)
{
	/*
	 *	Purpose:
	 *		Compute operators for elements which are tensor-products of lower dimensional elements.
	 *
	 *	Comments:
	 *		Add comment about general idea of how this is done. (ToBeModified)
	 */

	// Returned operators
	unsigned int *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS, **NfnIs, **NfnIc;
	double       **w_vIs, **w_vIc,
	             ****ChiS_vP, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS,
	             ****I_vGs_vP, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vS,
	             ****I_vGc_vP, ****I_vGc_vCc,                ****I_vGc_vS,
	             ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vIs, ****I_vCc_vIc,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fIs, ****ChiS_fIc,
	             ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fIs, ****I_vCc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FF, ****Ic_Weak_FF,
	             *****Ds_Weak_VV, *****Dc_Weak_VV;

	// Initialize DB Parameters
	unsigned int Adapt = DB.Adapt,
	             PP    = DB.PP,
	             PMax  = DB.PMax,
	             EFE   = DB.EFE;

	// Standard datatypes
	unsigned int dim, P, f, fh, Pb, PbMin, PbMax, *fhMax,
	             Eclass, dE, Nf, Vf, *Nfref, *ones_Nf,
	             NIn[3], NOut[3];
	double       *OP[3];

	struct S_ELEMENT *ELEMENT, *ELEMENTclass[2];

	// silence
	fhMax = NULL;
	ELEMENTclass[1] = NULL;

	Eclass = get_Eclass(EType);
	ELEMENT = get_ELEMENT_type(EType);

	ELEMENTclass[0] = ELEMENT->ELEMENTclass[0];
	if (EType == WEDGE)
		ELEMENTclass[1] = ELEMENT->ELEMENTclass[1];

	dE    = ELEMENT->d;
	Nf    = ELEMENT->Nf;
	Nfref = ELEMENT->Nfref;

	ones_Nf = malloc(Nf * sizeof *ones_Nf); // free
	for (f = 0; f < Nf; f++)
		ones_Nf[f] = 1;

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnGc = ELEMENT->NvnGc;
	NvnCs = ELEMENT->NvnCs;
	NvnCc = ELEMENT->NvnCc;
	NvnIs = ELEMENT->NvnIs;
	NvnIc = ELEMENT->NvnIc;
	NvnS  = ELEMENT->NvnS;
	NfnIs = ELEMENT->NfnIs;
	NfnIc = ELEMENT->NfnIc;

	w_vIs = ELEMENT->w_vIs;
	w_vIc = ELEMENT->w_vIc;

	ChiS_vP    = ELEMENT->ChiS_vP;
	ChiS_vIs   = ELEMENT->ChiS_vIs;
	ChiS_vIc   = ELEMENT->ChiS_vIc;
	ChiInvS_vS = ELEMENT->ChiInvS_vS;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vCs = ELEMENT->I_vGs_vCs;
	I_vGs_vS  = ELEMENT->I_vGs_vS;
	I_vGc_vP  = ELEMENT->I_vGc_vP;
	I_vGc_vCc = ELEMENT->I_vGc_vCc;
	I_vGc_vS  = ELEMENT->I_vGc_vS;
	I_vCs_vIs = ELEMENT->I_vCs_vIs;
	I_vCs_vIc = ELEMENT->I_vCs_vIc;
	I_vCc_vIs = ELEMENT->I_vCc_vIs;
	I_vCc_vIc = ELEMENT->I_vCc_vIc;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vIs = ELEMENT->D_vGs_vIs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vIc = ELEMENT->D_vGc_vIc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	ChiS_fIs   = ELEMENT->ChiS_fIs;
	ChiS_fIc   = ELEMENT->ChiS_fIc;

	I_vGs_fIs = ELEMENT->I_vGs_fIs;
	I_vGs_fIc = ELEMENT->I_vGs_fIc;
	I_vGc_fIs = ELEMENT->I_vGc_fIs;
	I_vGc_fIc = ELEMENT->I_vGc_fIc;
	I_vCs_fIs = ELEMENT->I_vCs_fIs;
	I_vCs_fIc = ELEMENT->I_vCs_fIc;
	I_vCc_fIs = ELEMENT->I_vCc_fIs;
	I_vCc_fIc = ELEMENT->I_vCc_fIc;

	Is_Weak_VV = ELEMENT->Is_Weak_VV;
	Ic_Weak_VV = ELEMENT->Ic_Weak_VV;
	Is_Weak_FF = ELEMENT->Is_Weak_FF;
	Ic_Weak_FF = ELEMENT->Ic_Weak_FF;
	Ds_Weak_VV = ELEMENT->Ds_Weak_VV;
	Dc_Weak_VV = ELEMENT->Dc_Weak_VV;

	switch (Adapt) {
		default: // ADAPT_H or ADAPT_HP
			fhMax = Nfref;
			break;
		case ADAPT_0:
		case ADAPT_P:
			fhMax = ones_Nf;
			break;
	}

	if (Eclass == C_TP) {
		NvnGs[0] = pow(ELEMENTclass[0]->NvnGs[0],dE);

		get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnP,ELEMENTclass[0]->I_vGs_vP[0][PP][0],
		                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
		I_vGs_vP[0][PP][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

		for (P = 0; P <= PMax; P++) {
			NvnGc[P] = pow(ELEMENTclass[0]->NvnGc[P],dE);
			NvnCs[P] = pow(ELEMENTclass[0]->NvnCs[P],dE);
			NvnCc[P] = pow(ELEMENTclass[0]->NvnCc[P],dE);
			NvnS[P]  = pow(ELEMENTclass[0]->NvnS[P],dE);

			switch (Adapt) {
				default: // ADAPT_P or ADAPT_HP
					if      (P == 0)    PbMin = P,   PbMax = P+1;
					else if (P == PMax) PbMin = P-1, PbMax = PMax;
					else                PbMin = P-1, PbMax = P+1;
					break;
				case ADAPT_0:
				case ADAPT_H:
					PbMin = P; PbMax = P;
					break;
			}
			for (Pb = PbMin; Pb <= PbMax; Pb++) {
				NvnIs[Pb] = pow(ELEMENTclass[0]->NvnIs[Pb],dE);
				NvnIc[Pb] = pow(ELEMENTclass[0]->NvnIc[Pb],dE);

				NfnIs[Pb][0] = pow(ELEMENTclass[0]->NvnIs[Pb],dE-1)*(ELEMENTclass[0]->NfnIs[Pb][0]);
				NfnIc[Pb][0] = pow(ELEMENTclass[0]->NvnIc[Pb],dE-1)*(ELEMENTclass[0]->NfnIc[Pb][0]);

				get_sf_parameters(1,ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->w_vIs[Pb],0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				w_vIs[Pb] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->w_vIc[Pb],0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				w_vIc[Pb] = sf_assemble_d(NIn,NOut,dE,OP);

				// Note: Most VOLUME operators need not interpolate between different orders
				if (P == Pb) {
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->ChiInvS_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiInvS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnP,ELEMENTclass[0]->ChiS_vP[P][PP][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiS_vP[P][PP][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->ChiS_vIs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->ChiS_vIc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnGc[Pb],ELEMENTclass[0]->I_vGs_vGc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vGc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vCs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->I_vGs_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnP,ELEMENTclass[0]->I_vGc_vP[P][PP][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vP[P][PP][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vCc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->I_vGc_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCs_vIs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCs_vIc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCc_vIs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCc_vIc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (EFE) {
						get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
										  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
						Is_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
										  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
						Ic_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					} else {
						;
					}

					for (dim = 0; dim < dE; dim++) {
						get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[P][Pb][0],
						                  ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->D_vGs_vCs[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGs_vCs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vGs_vIs[P][Pb][0],
						                  ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->D_vGs_vIs[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGs_vIs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
						                  ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->D_vGc_vCc[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],
						                  ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->D_vGc_vIc[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGc_vIc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->ICs[P][Pb][0],
						                  ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->D_vCs_vCs[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vCs_vCs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->ICc[P][Pb][0],
						                  ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->D_vCc_vCc[P][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vCc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						if (EFE) {
							get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
						                      ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->Ds_Weak_VV[P][Pb][0][0],
						                      NIn,NOut,OP,dE,dim,Eclass);
							Ds_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
						                      ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->Dc_Weak_VV[P][Pb][0][0],
						                      NIn,NOut,OP,dE,dim,Eclass);
							Dc_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						} else {
							;
						}
					}
				}

				for (f = 0; f < Nf; f++) {
				for (fh = 0; fh < fhMax[f]; fh++) {
					Vf = f*NFREFMAX+fh;
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->ChiS_vIs[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->ChiS_fIs[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->ChiS_vIc[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->ChiS_fIc[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					get_sf_parametersF(ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb],
					                   ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_FF[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					Is_Weak_FF[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb],
					                   ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_FF[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					Ic_Weak_FF[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (fh == 0) {
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->I_vGs_vIs[P][Pb],
										   ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->I_vGs_fIs[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->I_vGs_vIc[P][Pb],
										   ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->I_vGs_fIc[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->I_vGc_vIs[P][Pb],
										   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->I_vGc_fIs[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->I_vGc_vIc[P][Pb],
										   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->I_vGc_fIc[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->I_vCs_vIs[P][Pb],
										   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->I_vCs_fIs[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->I_vCs_vIc[P][Pb],
										   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->I_vCs_fIc[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->I_vCc_vIs[P][Pb],
										   ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->I_vCc_fIs[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->I_vCc_vIc[P][Pb],
										   ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->I_vCc_fIc[P][Pb],
										   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}
				}}
			}
		}
	} else if (Eclass == C_WEDGE) {
		unsigned int NOut0, NOut1;
		double       *OP0, *OP1, **OPF0, **OPF1;

		NvnGs[0] = (ELEMENTclass[0]->NvnGs[0])*(ELEMENTclass[1]->NvnGs[0]);

		get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnP,ELEMENTclass[0]->I_vGs_vP[0][PP][0],
		                  ELEMENTclass[1]->NvnGs[0],ELEMENTclass[1]->NvnP,ELEMENTclass[1]->I_vGs_vP[0][PP][0],
		                  NIn,NOut,OP,dE,3,Eclass);
		I_vGs_vP[0][PP][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

		for (P = 0; P <= PMax; P++) {
			NvnGc[P] = (ELEMENTclass[0]->NvnGc[P])*(ELEMENTclass[1]->NvnGc[P]);
			NvnCs[P] = (ELEMENTclass[0]->NvnCs[P])*(ELEMENTclass[1]->NvnCs[P]);
			NvnCc[P] = (ELEMENTclass[0]->NvnCc[P])*(ELEMENTclass[1]->NvnCc[P]);
			NvnS[P]  = (ELEMENTclass[0]->NvnS[P])*(ELEMENTclass[1]->NvnS[P]);

			switch (Adapt) {
				default: // ADAPT_P or ADAPT_HP
					if      (P == 0)    PbMin = P,   PbMax = P+1;
					else if (P == PMax) PbMin = P-1, PbMax = PMax;
					else                PbMin = P-1, PbMax = P+1;
					break;
				case ADAPT_0:
				case ADAPT_H:
					PbMin = P; PbMax = P;
					break;
			}
			for (Pb = 0; Pb <= PbMax; Pb++) {
				NvnIs[P] = (ELEMENTclass[0]->NvnIs[P])*(ELEMENTclass[1]->NvnIs[P]);
				NvnIc[P] = (ELEMENTclass[0]->NvnIc[P])*(ELEMENTclass[1]->NvnIc[P]);

				NfnIs[P][0] = (ELEMENTclass[0]->NfnIs[P][0])*(ELEMENTclass[1]->NvnIs[P]);
				NfnIs[P][1] = (ELEMENTclass[0]->NvnIs[P])*   (ELEMENTclass[1]->NfnIs[P][0]);
				NfnIc[P][0] = (ELEMENTclass[0]->NfnIc[P][0])*(ELEMENTclass[1]->NvnIc[P]);
				NfnIc[P][1] = (ELEMENTclass[0]->NvnIc[P])*   (ELEMENTclass[1]->NfnIc[P][0]);

				get_sf_parameters(1,ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->w_vIs[Pb],
				                  1,ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->w_vIs[Pb],NIn,NOut,OP,dE,3,Eclass);
				w_vIs[Pb] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->w_vIc[Pb],
				                  1,ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->w_vIc[Pb],NIn,NOut,OP,dE,3,Eclass);
				w_vIc[Pb] = sf_assemble_d(NIn,NOut,dE,OP);

				// Note: Most VOLUME operators need not interpolate between different orders
				if (P == Pb) {
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->ChiInvS_vS[P][Pb][0],
									  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[Pb],ELEMENTclass[1]->ChiInvS_vS[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					ChiInvS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnP,ELEMENTclass[0]->ChiS_vP[P][PP][0],
									  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnP,ELEMENTclass[1]->ChiS_vP[P][PP][0],
									  NIn,NOut,OP,dE,3,Eclass);
					ChiS_vP[P][PP][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->ChiS_vIs[P][Pb][0],
									  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->ChiS_vIs[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->ChiS_vIc[P][Pb][0],
									  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->ChiS_vIc[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnGc[Pb],ELEMENTclass[0]->I_vGs_vGc[P][Pb][0],
									  ELEMENTclass[1]->NvnGs[0],ELEMENTclass[1]->NvnGc[Pb],ELEMENTclass[1]->I_vGs_vGc[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vGc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[P][Pb][0],
									  ELEMENTclass[1]->NvnGs[0],ELEMENTclass[1]->NvnCs[Pb],ELEMENTclass[1]->I_vGs_vCs[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vCs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->I_vGs_vS[P][Pb][0],
									  ELEMENTclass[1]->NvnGs[0],ELEMENTclass[1]->NvnS[Pb],ELEMENTclass[1]->I_vGs_vS[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnP,ELEMENTclass[0]->I_vGc_vP[P][PP][0],
									  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnP,ELEMENTclass[1]->I_vGc_vP[P][PP][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vP[P][PP][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
									  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnCc[Pb],ELEMENTclass[1]->I_vGc_vCc[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vCc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->I_vGc_vS[P][Pb][0],
									  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnS[Pb],ELEMENTclass[1]->I_vGc_vS[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCs_vIs[P][Pb][0],
									  ELEMENTclass[1]->NvnCs[P],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->I_vCs_vIs[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCs_vIc[P][Pb][0],
									  ELEMENTclass[1]->NvnCs[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->I_vCs_vIc[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vCs_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vCc_vIs[P][Pb][0],
									  ELEMENTclass[1]->NvnCc[P],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->I_vCc_vIs[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vCc_vIc[P][Pb][0],
									  ELEMENTclass[1]->NvnCc[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->I_vCc_vIc[P][Pb][0],
									  NIn,NOut,OP,dE,3,Eclass);
					I_vCc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					for (dim = 0; dim < dE; dim++) {
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGs_vCs[P][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGs_vCs[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGs_vCs[P][Pb][0],      OP1 = ELEMENTclass[1]->D_vGs_vCs[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnCs[Pb],OP0,
										  ELEMENTclass[1]->NvnGs[0],ELEMENTclass[1]->NvnCs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGs_vCs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGs_vIs[P][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGs_vIs[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGs_vIs[P][Pb][0],      OP1 = ELEMENTclass[1]->D_vGs_vIs[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGs[0],ELEMENTclass[0]->NvnIs[Pb],OP0,
										  ELEMENTclass[1]->NvnGs[0],ELEMENTclass[1]->NvnIs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGs_vIs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGc_vCc[P][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGc_vCc[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],      OP1 = ELEMENTclass[1]->D_vGc_vCc[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],OP0,
										  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnCc[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGc_vIc[P][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGc_vIc[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],      OP1 = ELEMENTclass[1]->D_vGc_vIc[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],OP0,
										  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnIc[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGc_vIc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vCs_vCs[P][Pb][0][dim], OP1 = ELEMENTclass[1]->ICs[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->ICs[P][Pb][0],            OP1 = ELEMENTclass[1]->D_vCs_vCs[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnCs[Pb],OP0,
										  ELEMENTclass[1]->NvnCs[P],ELEMENTclass[1]->NvnCs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vCs_vCs[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vCc_vCc[P][Pb][0][dim], OP1 = ELEMENTclass[1]->ICc[P][Pb][0];
						else         OP0 = ELEMENTclass[0]->ICc[P][Pb][0],            OP1 = ELEMENTclass[1]->D_vCc_vCc[P][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnCc[Pb],OP0,
										  ELEMENTclass[1]->NvnCc[P],ELEMENTclass[1]->NvnCc[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vCc_vCc[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						if (EFE) {
							if (dim < 2) OP0 = ELEMENTclass[0]->Ds_Weak_VV[P][Pb][0][dim], OP1 = ELEMENTclass[1]->Is_Weak_VV[P][Pb][0];
							else         OP0 = ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],      OP1 = ELEMENTclass[1]->Ds_Weak_VV[P][Pb][0][0];
							get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],OP0,
											  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
							Ds_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							if (dim < 2) OP0 = ELEMENTclass[0]->Dc_Weak_VV[P][Pb][0][dim], OP1 = ELEMENTclass[1]->Ic_Weak_VV[P][Pb][0];
							else         OP0 = ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],      OP1 = ELEMENTclass[1]->Dc_Weak_VV[P][Pb][0][0];
							get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],OP0,
											  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIc[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
							Dc_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						} else {
							;
						}
					}
				}

				for (f = 0; f < Nf; f++) {
				for (fh = 0; fh < fhMax[f]; fh++) {
					Vf = f*NFREFMAX+fh;

					if (f < 3) { OPF0  = ELEMENTclass[0]->ChiS_fIs[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_vIs[P][Pb];
								 NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],    NOut1 = ELEMENTclass[1]->NvnIs[Pb];
					} else {     OPF0  = ELEMENTclass[0]->ChiS_vIs[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_fIs[P][Pb];
								 NOut0 = ELEMENTclass[0]->NvnIs[Pb],       NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],NOut0,OPF0,
									   ELEMENTclass[1]->NvnS[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					if (f < 3) { OPF0  = ELEMENTclass[0]->ChiS_fIc[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_vIc[P][Pb];
								 NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],    NOut1 = ELEMENTclass[1]->NvnIc[Pb];
					} else {     OPF0  = ELEMENTclass[0]->ChiS_vIc[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_fIc[P][Pb];
								 NOut0 = ELEMENTclass[0]->NvnIc[Pb],       NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],NOut0,OPF0,
									  ELEMENTclass[1]->NvnS[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					if (fh == 0) {
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGs_vIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[0],NOut0,OPF0,
										  ELEMENTclass[1]->NvnGs[0],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGs_vIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[0],NOut0,OPF0,
										  ELEMENTclass[1]->NvnGs[0],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGc_fIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_vIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGc_vIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_fIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],NOut0,OPF0,
										  ELEMENTclass[1]->NvnGc[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGc_fIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_vIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGc_vIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_fIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],NOut0,OPF0,
										  ELEMENTclass[1]->NvnGc[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCs_fIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_vIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCs_vIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_fIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],NOut0,OPF0,
										  ELEMENTclass[1]->NvnCs[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCs_fIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_vIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCs_vIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_fIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],NOut0,OPF0,
										  ELEMENTclass[1]->NvnCs[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCc_fIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_vIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCc_vIs[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_fIs[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],NOut0,OPF0,
										  ELEMENTclass[1]->NvnCc[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCc_fIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_vIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCc_vIc[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_fIc[P][Pb];
									 NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],NOut0,OPF0,
										  ELEMENTclass[1]->NvnCc[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}
				}}
			}
		}
	} else {
		printf("Error: Unsupported Eclass in setup_TP_operators.\n"), exit(1);
	}
	free(ones_Nf);
}

static void setup_ELEMENT_FACET_ordering(const unsigned int FType)
{
	// Returned operators
	unsigned int ***nOrd_fIs, ***nOrd_fIc;

	// Initialize DB Parameters
	unsigned int d              = DB.d,
	             PMax           = DB.PMax,
	             **PIfs         = DB.PIfs,
	             **PIfc         = DB.PIfc;
	char         ***NodeTypeIfs = DB.NodeTypeIfs,
	             ***NodeTypeIfc = DB.NodeTypeIfc;

	// Standard datatypes
	unsigned int P, dE,
	             NfnIs, NfnIc, NsIs, NsIc,
	             IndOrd, NOrd, Eclass,
	             *symms_fIs, *symms_fIc;
	double       *rst_fIs, *rst_fIc, *w;

	struct S_ELEMENT *ELEMENT;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	select_functions(&basis,&grad_basis,&cubature,FType);

	// silence
	NfnIs = NfnIc = NsIs = NsIc = 0;

	ELEMENT = get_ELEMENT_type(FType);
	Eclass  = get_Eclass(FType);

	dE = d-1;

	if      (FType == POINT) NOrd = 1;
	else if (FType == LINE)  NOrd = 2;
	else if (FType == QUAD)  NOrd = 8;
	else if (FType == TRI)   NOrd = 6;
	else    printf("Error: Unsupported FType in setup_ELEMENT_FACET_ordering.\n"), exit(1);

	nOrd_fIs = ELEMENT->nOrd_fIs;
	nOrd_fIc = ELEMENT->nOrd_fIc;

	for (P = 0; P <= PMax; P++) {
		cubature(&rst_fIs,&w,&symms_fIs,&NfnIs,&NsIs,0,PIfs[P][Eclass],dE,NodeTypeIfs[P][Eclass]); // free
		cubature(&rst_fIc,&w,&symms_fIc,&NfnIc,&NsIc,0,PIfc[P][Eclass],dE,NodeTypeIfc[P][Eclass]); // free

		for (IndOrd = 0; IndOrd < NOrd; IndOrd++) {
			nOrd_fIs[P][IndOrd] = malloc(NfnIs * sizeof ***nOrd_fIs); //  keep
			nOrd_fIc[P][IndOrd] = malloc(NfnIc * sizeof ***nOrd_fIc); //  keep

			get_facet_ordering(d,IndOrd,FType,NfnIs,NsIs,symms_fIs,rst_fIs,nOrd_fIs[P][IndOrd]);
			get_facet_ordering(d,IndOrd,FType,NfnIc,NsIc,symms_fIc,rst_fIc,nOrd_fIc[P][IndOrd]);
		}

		free(rst_fIs);
		free(rst_fIc);
		free(symms_fIs);
		free(symms_fIc);
	}
}

void setup_operators(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int EType;

	// POINT
	EType = POINT;
	if (d == 1)
		setup_ELEMENT_FACET_ordering(EType);

	// LINE (Includes TP Class)
	EType = LINE;
	setup_ELEMENT_VeF(EType);
	setup_ELEMENT_plotting(EType);
	setup_ELEMENT_normals(EType);
	setup_ELEMENT_operators(EType);
	if (d == 2)
		setup_ELEMENT_FACET_ordering(EType);

	// QUAD
	EType = QUAD;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_TP_operators(EType);
		if (d == 3)
			setup_ELEMENT_FACET_ordering(EType);
	}

	// HEX
	EType = HEX;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_TP_operators(EType);
	}

	// TRI
	EType = TRI;
	if (is_ELEMENT_present(EType)) {
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);
		if (d == 3)
			setup_ELEMENT_FACET_ordering(EType);
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
		setup_TP_operators(EType);
	}

	// Free unused operators (including unused lower dimensional operators (if applicable) and sum factorized operators)
}
