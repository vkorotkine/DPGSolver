// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_operators.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
 
#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_OpCSR.h"

#include "bases.h"
#include "cubature.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "sum_factorization.h"
#include "array_free.h"
#include "plotting_element_info.h"
#include "solver_functions.h"
#include "adaptation.h"
#include "memory_destructors.h"

/*
 *	Purpose:
 *		Set up operators to be used throughout the code.
 *
 *	Comments:
 *		Different operators are set up depending on which higher-dimensional element dependencies are present as well as
 *		whether adaptivity is desired:
 *			LINE  : Standard operators + additional SF operators if QUADs are present (standard assembly)
 *			TRI   : Standard operators + additional SF operators if WEDGEs are present (standard assembly)
 *			QUAD  : Standard operators + sparse versions of heavily used operators (TP asssembly)
 *			TET   : Standard operators (standard assembly)
 *			HEX   : Standard operators + sparse versions of heavily used operators (TP asssembly)
 *			WEDGE : Standard operators + sparse versions of heavily used operators (TP asssembly)
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
 *		order to improve code readability.
 *		Ensure that operators for hp refinement are only stored when refinement is enabled (ToBeDeleted).
 *
 *		For the PYR element rst_vGs != E_rst_vC because of the TP structures of the TP nodes and the rotational symmetry
 *		ordering of the PYR nodes in layers of 't'. This means that special consideration must be made if attempting to
 *		use the rotational symmetry of rst_vGs PYR nodes.
 *
 *		The most intuitive method for obtaining projections between different h-refined spaces is to use the
 *		(Q)uadrature (M)irror (F)ilter coefficient matrices from the multiresolution framework. As the code discards the
 *		detail coefficients when projecting from fine to coarser spaces, the wavelet basis functions are not needed.
 *		Furthermore, the information embedded in the detail coefficients is simply an elegant representation of lack of
 *		smoothness.
 *		When imagining the use of a modal basis in the discussion of Kopriva(1996) in the context of h-adaptation
 *		related L2 projections, it is easily seen that the QMF coefficient matrices perform exactly such projections.
 *		See Alpert(2002) for their definitions.
 *
 *		Consider removing the "-ve" sign included in the I()_Weak_FF operator to avoid confusion. (ToBeModified)
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
 *		Ihat_(1)(2)(3)_(4)(5)(6) : Analogous to I (above) but where interpolation is of coefficients instead of values.
 *		  Note: values and coefficients are the same for the nodal scheme, but not for modal.
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
 *		*_sp : Standard (*) operator stored in (sp)arse format ((C)ompressed (S)parse (R)ow), only used for TP ELEMENTs.
 *
 *	setup_L2_projection_operators:
 *		L2hat_vS_vS : L2 projection operator from coefficients of the VOLUMEs on level l+1 to the VOLUME on level l.
 *
 *	setup_ELEMENT_VeV:
 *	setup_ELEMENT_VeF:
 *		VeV : (Ve)rtex to (V)OLUME operator used to project VOLUME vertex nodes to VOLUME vertex nodes for all supported
 *		      h-refinements.
 *		VeF : (Ve)rtex to (F)ACET operator used to project VOLUME vertex nodes to FACET vertex nodes for all supported
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
 *			**   : 0:NREFMAX
 *			^    : 0:NFMAX
 *			^^   : 0:NFREFMAX*NFMAX
 *
 *			rP   : (r)ange P  == 0:PMax
 *			P(P) : P as a function of P (i.e. for a range in field 1 of [0:PMax] -> [0][0], [1][1], ..., [P][P])
 *			rPb  : (r)ange Pb (always a function of P) == PbMin:PbMax (PbMin = min(P-1,0), PbMax = max(P+1,PMax))
 *			rd   : (r)ange (d)imension == 0:d-1
 *
 *		                   Adapt
 *		Returned Op (SF) | ADAPT_0         ADAPT_P             ADAPT_H         ADAPT_HP
 *		                 |
 *		NvnGs            | [1]             [1]                 [1]             [1]
 *		NvnGc            | [P]             [rP]                [P]             [rP]
 *		NvnCs            |
 *		NvnCc            |
 *		NvnIs            |
 *		NvnIc            |
 *		NfnIs            | [P][0:1]        [rP][0:1]           [P][2]          [rP][0:1]
 *		NfnIc            |
 *		                 |
 *		w_vIs            | [P]             [rP]                [P]             [rP]
 *		w_vIc            |
 *
 *		ChiS_vP          | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		ChiS_vS      (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		ChiS_vIs     (*) |
 *		ChiS_vIc     (*) |
 *		ChiInvS_vS       | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		                 |
 *		ICs              | [P][P][0]       [rP][P(P)][0]       [P][P][0]       [rP][P(P)][0]
 *		ICc              |
 *		                 |
 *		I_vGs_vP         | [1][P][0]       [1][rP][0]          [1][P][0]       [1][rP][0]
 *		I_vGs_vGs        | NOT_USED        NOT_USED            [1][1][**]      [1][1][**]
 *		I_vGs_vGc        | [1][P][0]       [1][rP][0]          [1][P][0]       [1][rP][0]
 *		I_vGs_vCs        |
 *		I_vGs_vS     (*) | [1][P][0]       [1][rPb][0]         [1][P][*]       [1][rP][*]
 *		I_vGs_vIs    (*) |
 *		I_vGs_vIc    (*) |
 *		I_vGc_vP         | [1][P][0]       [1][P(P)][0]        [1][P][0]       [1][P(P)][0]
 *		I_vGc_vCc        |
 *		I_vGc_vS     (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		I_vGc_vIs    (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		I_vGc_vIc    (*) |
 *		I_vCs_vS     (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		I_vCs_vIs    (*) |
 *		I_vCs_vIc    (*) |
 *		I_vCc_vS     (*) |
 *		I_vCc_vIs    (*) |
 *		I_vCc_vIc    (*) |
 *		                 |
 *		Ihat_vS_vS       | NOT_USED        [rP][rPb][0]        [P][P][*]       [rP][rP][0] & [P][P][*]
 *		                 |
 *		D_vGs_vCs        | [1][P][0][rd]   [1][rP][0][rd]      [1][P][0][rd]   [1][rP][0][rd]
 *		D_vGs_vIs        |
 *		D_vGc_vCc        | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		D_vGc_vIc        |
 *		D_vCs_vCs        |
 *		D_vCc_vCc        |
 *		                 |
 *		ChiS_fIs         | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rP][^^]
 *		ChiS_fIc         |
 *		                 |
 *		I_vGs_fS         | [1][P][^]       [P][rP][^]          [1][P][^^]      [1][rP][^^]
 *		I_vGs_fIs        |
 *		I_vGs_fIc        |
 *		I_vGc_fS         | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rP][^^]
 *		I_vGc_fIs        |
 *		I_vGc_fIc        |
 *		I_vCs_fIs        |
 *		I_vCs_fIc        |
 *		I_vCc_fIs        |
 *		I_vCc_fIc        |
 *		                 |
 *		Is_Weak_VV   (*) | [P][P][0]       [rP][rPb][0]        [P][P][*]       [rP][rP][*]
 *		Ic_Weak_VV   (*) |
 *		Ds_Weak_VV       | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		Dc_Weak_VV       |
 *		                 |
 *		Is_Weak_FF       | [P][P][^]       [rP][rPb][^]       [P][P][^^]       [rP][rP][^^]
 *		Ic_Weak_FF       |
 *		                 |
 *		ChiS_fIs_sp      | [P][P][^]       [rP][rPb][^]        [P][P][^^]      [rP][rP][^^]
 *		ChiS_fIc_sp      |
 *		Ds_Weak_VV_sp    | [P][P][0][rd]   [rP][P(P)][0][rd]   [P][P][0][rd]   [rP][P(P)][0][rd]
 *		Dc_Weak_VV_sp    |
 *		Is_Weak_FF_sp    | [P][P][^]       [rP][rPb][^]       [P][P][^^]       [rP][rP][^^]
 *		Ic_Weak_FF_sp    |
 *
 *
 *	References:
 *		Zwanenburg(2016)_Equivalence between the Energy Stable Flux Reconstruction and Filtered Discontinuous Galerkin Schemes
 *		Karniadakis(1999)_Spectral hp Element Methods for CFD
 *		Hesthaven(2000)_Stable Spectral Methods on Tetrahedral Elements
 *		Stiller(2008)_Factorization Techniques for Nodal Spectral Elements in Curved Domains
 *		Alpert(2002)-Adaptive_Solution_of_Partial_Differential_Equations_in_Multiwavelet_Bases
 *		Kopriva(1996)-A_Conservative_Staggered-Grid_Chebyshev_Multidomain_Method_for_Compressible_Flows_II._A_Semi-Structured_Method
 *
 *		ToBeDeleted: Potentially add in MATH578 report instead of Alpert(2002) for the QMF discussion.
 */

struct S_BCOORDS {
	unsigned int Nve,
	             *NfnS, *NfnIs, *NfnIc;
	double       **w_fIs, **w_fIc, **BCoords_S, **BCoords_Is, **BCoords_Ic;
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
		printf("Error: WEDGE elements use a combination of TRI and LINE basis functions/nodes.\n"), EXIT_MSG;
		break;
	case PYR:
		*basis      = basis_PYR;
		*grad_basis = grad_basis_PYR;
		*cubature   = cubature_PYR;
		break;
	default:
		printf("Error: Unsupported type = %d.\n",type), EXIT_MSG;
		break;
	}
}

static void setup_ELEMENT_plotting(const unsigned int EType)
{
	// Standard datatypes
	unsigned int P, PSMin, PSMax, NvnP, NE, u1 = 1,
	             *connectivity, *types;
	double       *rst_vP;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {
		plotting_element_info(&rst_vP,&connectivity,&types,&NvnP,&NE,max(P,u1),EType); // free
		ELEMENT->connectivity[P]  = connectivity;
		ELEMENT->connect_types[P] = types;
		ELEMENT->connect_NE[P]    = NE;
		ELEMENT->NvnP[P] = NvnP;

		free(rst_vP);
	}
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
		printf("Error: Unsupported EType.\n"), EXIT_MSG;
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
		printf("Error: Unsupported EType.\n"), EXIT_MSG;
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
	             *NfnS, *NfnIs, *NfnIc;
	double       **w_fIs, **w_fIc,
	             **BCoords_S, **BCoords_Is, **BCoords_Ic, *one;

	struct S_BCOORDS *BCoords_dEm1;
	struct S_ELEMENT *ELEMENT_F;

	BCoords_dEm1 = malloc(sizeof *BCoords_dEm1);      // keep (requires external free)
	NfnS         = malloc(NP * sizeof *(NfnS));       // keep (requires external free)
	NfnIs        = malloc(NP * sizeof *(NfnIs));      // keep (requires external free)
	NfnIc        = malloc(NP * sizeof *(NfnIc));      // keep (requires external free)
	w_fIs        = malloc(NP * sizeof *(w_fIs));      // keep (requires external free)
	w_fIc        = malloc(NP * sizeof *(w_fIc));      // keep (requires external free)
	BCoords_S    = malloc(NP * sizeof *(BCoords_S));  // keep (requires external free)
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
			NfnS[P]  = 1;
			NfnIs[P] = 1;
			NfnIc[P] = 1;

			w_fIs[P] = mm_Alloc_d(CBCM,CBT,CBT,1,1,1,1.0,one,one); // keep
			w_fIc[P] = mm_Alloc_d(CBCM,CBT,CBT,1,1,1,1.0,one,one); // keep

			BCoords_S[P]  = malloc(1 * sizeof **(BCoords_S));
			BCoords_Is[P] = malloc(1 * sizeof **(BCoords_Is));
			BCoords_Ic[P] = malloc(1 * sizeof **(BCoords_Ic));

			BCoords_S[P][0]  = 1.0;
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
			printf("Error: Unsupported EType/IndFType combination.\n"), EXIT_MSG;
		}

		// Initialize DB Parameters
		unsigned int **PIfs         = DB.PIfs,
		             **PIfc         = DB.PIfc;
		char         **NodeTypeG    = DB.NodeTypeG,
		             ***NodeTypeS   = DB.NodeTypeS,
		             ***NodeTypeIfs = DB.NodeTypeIfs,
		             ***NodeTypeIfc = DB.NodeTypeIfc;

		// Standard datatypes
		unsigned int dE, Nbf,
		             NfnGs, EType, Eclass,
		             dummy_ui, *dummyPtr_ui;
		double       *rst_vGs, *rst_fS, *rst_fIs, *rst_fIc,
		             *dummyPtr_d[2],
		             *IGs,
		             *ChiRefGs_vGs,
		             *ChiRefInvGs_vGs,
		             *ChiRefGs_fS, *ChiRefGs_fIs, *ChiRefGs_fIc;

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
			cubature(&rst_fS,&dummyPtr_d[0],&dummyPtr_ui,&NfnS[P], &dummy_ui,0,P,              dE,NodeTypeS[P][Eclass]);   free(dummyPtr_ui); // free
			cubature(&rst_fIs,&w_fIs[P],    &dummyPtr_ui,&NfnIs[P],&dummy_ui,1,PIfs[P][Eclass],dE,NodeTypeIfs[P][Eclass]); free(dummyPtr_ui); // free
			cubature(&rst_fIc,&w_fIc[P],    &dummyPtr_ui,&NfnIc[P],&dummy_ui,1,PIfc[P][Eclass],dE,NodeTypeIfc[P][Eclass]); free(dummyPtr_ui); // free

			ChiRefGs_fS  = basis(1,rst_fS, NfnS[P], &Nbf,dE); // free
			ChiRefGs_fIs = basis(1,rst_fIs,NfnIs[P],&Nbf,dE); // free
			ChiRefGs_fIc = basis(1,rst_fIc,NfnIc[P],&Nbf,dE); // free

			BCoords_S[P]  = mm_Alloc_d(CBCM,CBT,CBT,NfnS[P], NfnGs,NfnGs,1.0,ChiRefGs_fS, ChiRefInvGs_vGs); // keep
			BCoords_Is[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnIs[P],NfnGs,NfnGs,1.0,ChiRefGs_fIs,ChiRefInvGs_vGs); // keep
			BCoords_Ic[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnIc[P],NfnGs,NfnGs,1.0,ChiRefGs_fIc,ChiRefInvGs_vGs); // keep

			free(rst_fS);
			free(rst_fIs);
			free(rst_fIc);

			free(ChiRefGs_fS);
			free(ChiRefGs_fIs);
			free(ChiRefGs_fIc);
		}
		free(ChiRefInvGs_vGs);
	}

	free(one);

	BCoords_dEm1->Nve   = Nve;
	BCoords_dEm1->NfnS  = NfnS;
	BCoords_dEm1->NfnIs = NfnIs;
	BCoords_dEm1->NfnIc = NfnIc;
	BCoords_dEm1->w_fIs = w_fIs;
	BCoords_dEm1->w_fIc = w_fIc;
	BCoords_dEm1->BCoords_S  = BCoords_S;
	BCoords_dEm1->BCoords_Is = BCoords_Is;
	BCoords_dEm1->BCoords_Ic = BCoords_Ic;

	return BCoords_dEm1;
}

static void setup_ELEMENT_VeV(const unsigned int EType)
{
	/*
	 *	Comments:
	 *		For VeVref_TET, note that there are 3 options for the internal TET positioning. Here only one option is
	 *		supported but the others can be added if, for example, it is desired to refine TETs by splitting along the
	 *		longest physical diagonal. Similarly, if it is found that TETs are advantageous over PYRs, the top two PYRs
	 *		in the refined PYR can be split into four TETs.
	 */

	// Standard datatypes
	unsigned int i, j, jMax, VeVrefInd, Nve, *Nvve, Nvref, EcType;
	double       **VeV;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	Nve   = ELEMENT->Nve;
	Nvve  = ELEMENT->Nvve;
	VeV   = ELEMENT->VeV;
	Nvref = ELEMENT->Nvref;

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
	       VeVref_TET[144]  = {1.0 , 0.0 , 0.0 , 0.0 , // Conforming TET
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 , // Corner TETs
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 , // 2
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 , // 3
	                           0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.5 , // 4
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 , // Internal TETs
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 , // 6
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 , // 7
	                           0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 0.5 , 0.5 , 0.0 , // 8
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 },
	       VeVref_PYR[255]  = {1.0 , 0.0 , 0.0 , 0.0 , 0.0 , // Conforming
	                           0.0 , 1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 , 0.0 , // First 4 PYRs
	                           0.5 , 0.5 , 0.0 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 , 0.0 , // 2
	                           0.0 , 1.0 , 0.0 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 , 0.0 , // 3
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 , // 4
	                           0.0 , 0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 , // 4 TETs
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 , // 6
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 , // 7
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.5 , 0.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 , // 8
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.5 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 , // Last 2 PYRs
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, 0.0 ,
	                           0.5 , 0.0 , 0.0 , 0.0 , 0.5 , // 10
	                           0.0 , 0.5 , 0.0 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 0.0 , 1.0 };


	switch(EType) {
	case LINE:
		for (i = 0; i < Nvref; i++)
			Nvve[i] = Nve;

		VeVrefInd = 0;
		for (i = 0; i < Nvref; i++) {
			jMax = Nve*Nvve[i];
			if (i)
				VeVrefInd += jMax;
			VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
			for (j = 0; j < jMax; j++)
				VeV[i][j] = VeVref_LINE[VeVrefInd+j];
		}
		break;
	case TRI:
		for (i = 0; i < Nvref; i++)
			Nvve[i] = Nve;

		VeVrefInd = 0;
		for (i = 0; i < Nvref; i++) {
			jMax = Nve*Nvve[i];
			if (i)
				VeVrefInd += jMax;
			VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
			for (j = 0; j < jMax; j++)
				VeV[i][j] = VeVref_TRI[VeVrefInd+j];
		}
		break;
	case TET:
		// Original TET, 8 TETs (4 corner TETs + 2 PYRs divided into 2 TETs each)
		for (i = 0; i < Nvref; i++)
			Nvve[i] = Nve;

		VeVrefInd = 0;
		for (i = 0; i < Nvref; i++) {
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
		for (i = 0; i < Nvref; i++) {
			EcType = get_VOLUMEc_type(EType,i);
			if      (EcType == PYR) Nvve[i] = Nve;
			else if (EcType == TET) Nvve[i] = 4;
			else
				printf("Error: Unsupported EcType (PYR).\n"), EXIT_MSG;
		}

		VeVrefInd = 0;
		for (i = 0; i < Nvref; i++) {
			if (i)
				VeVrefInd += jMax;
			jMax = Nve*Nvve[i];
			VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
			for (j = 0; j < jMax; j++)
				VeV[i][j] = VeVref_PYR[VeVrefInd+j];
		}
		break;
	case QUAD:
	case HEX:
	case WEDGE:
	default:
		printf("Error: Unsupported EType.\n"), EXIT_MSG;
		break;
	}
}

static double get_L2_scaling(const unsigned int EType, const unsigned int vref)
{
	/*
	 *	Purpose:
	 *		Returns the VOLUME ratio of the refined VOLUME as compared to the initial VOLUME.
	 */

	if (vref == 0)
		return 1.0;

	switch (EType) {
	case LINE:
		return 0.5;
		break;
	case TRI:
		if (vref < 5)
			return 0.25;
		else
			printf("Error: Unsupported vref (%d) (anisotropic) for TRIs.\n",vref), EXIT_MSG;
		break;
	case TET:
		if (vref < 9)
			return 0.125;
		else
			printf("Error: Unsupported vref (%d) (anisotropic) for TETs.\n",vref), EXIT_MSG;
		break;
	case PYR:
		if (vref < 5 || vref > 8)
			return 0.125;
		else
			return 0.125; // 1/8 of reference TET
		break;
	case QUAD:
	case HEX:
	case WEDGE:
	default:
		printf("Error: Unsupported EType (%d).\n",EType), EXIT_MSG;
		break;
	}
}

static void setup_ELEMENT_operators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS, **NfnS, **NfnIs, **NfnIc;
	double       **w_vIs, **w_vIc,
	             ****ChiS_vP, ****ChiS_vS, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS,
	             ****ICs, ****ICc,
	             ****I_vGs_vP, ****I_vGs_vGs, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vIs, ****I_vGs_vIc, ****I_vGs_vS,
	             ****I_vGc_vP,                ****I_vGc_vCc, ****I_vGc_vIs, ****I_vGc_vIc, ****I_vGc_vS,
	             ****I_vCs_vS, ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vS, ****I_vCc_vIs, ****I_vCc_vIc,
	             ****Ihat_vS_vS,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fS, ****ChiS_fIs, ****ChiS_fIc,
				 *****GradChiS_fIs, *****GradChiS_fIc,
	             ****I_vGs_fS, ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fS, ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vCs_fS, ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fS, ****I_vCc_fIs, ****I_vCc_fIc,
				 *****D_vGs_fIs, *****D_vGs_fIc,
				 *****D_vGc_fIs, *****D_vGc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FF, ****Ic_Weak_FF,
	             *****Ds_Weak_VV, *****Dc_Weak_VV;

	struct S_OpCSR ****ChiS_fIs_sp, ****ChiS_fIc_sp,
	               ****Is_Weak_FF_sp, ****Ic_Weak_FF_sp;

	// Initialize DB Parameters
	unsigned int PMax         = DB.PMax,
	             PGlobal      = DB.PGlobal,
	             NP           = DB.NP,
	             PGs          = DB.PGs,
	             *PGc         = DB.PGc,
	             **PCs        = DB.PCs,
	             **PCc        = DB.PCc,
	             **PIvs       = DB.PIvs,
	             **PIvc       = DB.PIvc,
	             Collocated   = DB.Collocated,
	             EFE          = DB.EFE,
	             *VFPartUnity = DB.VFPartUnity;

	char         *BasisType     = DB.BasisType,
	             **NodeTypeG    = DB.NodeTypeG,
	             ***NodeTypeIvs = DB.NodeTypeIvs,
	             ***NodeTypeIvc = DB.NodeTypeIvc,
	             ***NodeTypeS   = DB.NodeTypeS;

	// Standard datatypes
	unsigned int i, iMax, dim, dE, P, f, vh, fh, Vf, IndFType, PSMin, PSMax, Pb, PbMin, PbMax, fhMax,
	             Nve, Nf, Nbf, Eclass, NFTypes, Nvref, vref, vrefSF, NvrefSF, *Nfref, *ones_Nf, NEhref,
	             Indh, NvnP, u1 = 1,
	             B_Nve[2], *Nfve, *Nvve, *EType_h,
	             dummy_ui, *dummyPtr_ui[2];
	double       *E_rst_vC, *rst_vC, **VeF, **VeV,
	             *rst_vP, *rst_vGs, *rst_vGc, *rst_vCs, *rst_vCc, **rst_vIs, **rst_vIc, **rst_vS,
	             *rst_fS, *rst_fIs, *rst_fIc,
	             *wInv_vIs, *wInv_vIc,
	             *diag_w_vIs, *diag_w_vIc, *diag_wInv_vIs, *diag_wInv_vIc,
	             ***w_fIs, ***w_fIc, *diag_w_fIs, *diag_w_fIc,
	             *IGs, *IGc, *IS,
	             *TGs, *TGc, *TCs, *TCc, *TS,
	             *ChiRefGc_vGc, *ChiRefCs_vCs, *ChiRefCc_vCc, *ChiRefS_vS,
	             *ChiGc_vGc,    *ChiCs_vCs,    *ChiCc_vCc,
	             *ChiRefInvGs_vGs, *ChiRefInvGc_vGc, *ChiRefInvCs_vCs, *ChiRefInvCc_vCc, *ChiRefInvS_vS,
	             *ChiInvGs_vGs,    *ChiInvGc_vGc,    *ChiInvCs_vCs,    *ChiInvCc_vCc,
	             *ChiRefGs_vP, *ChiRefGs_vGs, *ChiRefGs_vGc, *ChiRefGs_vCs, *ChiRefGs_vIs, *ChiRefGs_vIc, *ChiRefGs_vS,
	             *ChiRefGs_fS, *ChiRefGs_fIs, *ChiRefGs_fIc,
	             *ChiRefGc_vP,                *ChiRefGc_vCc, *ChiRefGc_vIs, *ChiRefGc_vIc, *ChiRefGc_vS,
	             *ChiRefGc_fS, *ChiRefGc_fIs, *ChiRefGc_fIc,
	             *ChiRefCs_vS, *ChiRefCs_vIs, *ChiRefCs_vIc, *ChiRefCs_fS, *ChiRefCs_fIs, *ChiRefCs_fIc,
	             *ChiRefCc_vS, *ChiRefCc_vIs, *ChiRefCc_vIc, *ChiRefCc_fS, *ChiRefCc_fIs, *ChiRefCc_fIc,
	             *ChiRefS_vP, *ChiRefS_vIs, *ChiRefS_vIc, *ChiRefS_fS, *ChiRefS_fIs, *ChiRefS_fIc,
	             *ChiGs_vP, *ChiGs_vGs, *ChiGs_vGc, *ChiGs_vCs, *ChiGs_vIs, *ChiGs_vIc, *ChiGs_vS,
	             *ChiGs_fS, *ChiGs_fIs, *ChiGs_fIc,
	             *ChiGc_vP,             *ChiGc_vCc, *ChiGc_vIs, *ChiGc_vIc, *ChiGc_vS,
	             *ChiGc_fS, *ChiGc_fIs, *ChiGc_fIc,
	             *ChiCs_vS, *ChiCs_vIs, *ChiCs_vIc, *ChiCs_fS, *ChiCs_fIs, *ChiCs_fIc,
	             *ChiCc_vS, *ChiCc_vIs, *ChiCc_vIc, *ChiCc_fS, *ChiCc_fIs, *ChiCc_fIc,
	             **GradChiRefGs_vCs, **GradChiRefGs_vIs,
	             **GradChiRefGc_vCc, **GradChiRefGc_vIc,
	             **GradChiRefCs_vCs, **GradChiRefCc_vCc,
	             **GradChiRefS_vIs,  **GradChiRefS_vIc,
				 **GradChiRefGs_fIs, **GradChiRefGs_fIc,
				 **GradChiRefGc_fIs, **GradChiRefGc_fIc,
				 **GradChiRefS_fIs,  **GradChiRefS_fIc,
	             **GradChiGs_vCs, **GradChiGs_vIs,
	             **GradChiGc_vCc, **GradChiGc_vIc,
	             **GradChiCs_vCs, **GradChiCc_vCc,
	             **GradChiS_vIs,  **GradChiS_vIc,
				 **GradChiGs_fIs, **GradChiGs_fIc,
				 **GradChiGc_fIs, **GradChiGc_fIc,
	             *dummyPtr_d;

	struct BCoords {
		double **Is, **Ic, **S;
	} *BCoords_F[2], **BCoords_V;

	struct S_BCOORDS *BCoords_dEm1[2];
	struct S_ELEMENT *ELEMENT, *ELEMENT_h;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	// silence
	ChiGs_vGs = NULL; ChiGc_vGc = NULL;
	ChiCs_vCs = NULL; ChiCc_vCc = NULL;
	ChiRefGs_vIs = ChiRefGs_vIc = NULL;

	ELEMENT   = get_ELEMENT_type(EType);
	ELEMENT_h = ELEMENT;

	// No need to consider the second Eclass as WEDGE basis functions will be built through a combination of lower
	// dimensional operators.
	Eclass = get_Eclass(EType);
//	if (EType != TET && EType != PYR)
		setup_ELEMENT_VeV(EType);

	dE      = ELEMENT->d;
	Nve     = ELEMENT->Nve;
	Nf      = ELEMENT->Nf;
	Nfref   = ELEMENT->Nfref;
	Nfve    = ELEMENT->Nfve;
	VeF     = ELEMENT->VeF;
	VeV     = ELEMENT->VeV;
	Nvve    = ELEMENT->Nvve;
	Nvref   = ELEMENT->Nvref;
	NvrefSF = ELEMENT->NvrefSF;
	NEhref  = ELEMENT->NEhref;
	EType_h = ELEMENT->type_h;

	select_functions(&basis,&grad_basis,&cubature,EType);

	// Stored operators
	NvnGs = ELEMENT->NvnGs;
	NvnGc = ELEMENT->NvnGc;
	NvnCs = ELEMENT->NvnCs;
	NvnCc = ELEMENT->NvnCc;
	NvnIs = ELEMENT->NvnIs;
	NvnIc = ELEMENT->NvnIc;
	NvnS  = ELEMENT->NvnS;
	NfnS  = ELEMENT->NfnS;
	NfnIs = ELEMENT->NfnIs;
	NfnIc = ELEMENT->NfnIc;

	w_vIs = ELEMENT->w_vIs;
	w_vIc = ELEMENT->w_vIc;

	ChiS_vP    = ELEMENT->ChiS_vP;
	ChiS_vS    = ELEMENT->ChiS_vS;
	ChiS_vIs   = ELEMENT->ChiS_vIs;
	ChiS_vIc   = ELEMENT->ChiS_vIc;
	ChiInvS_vS = ELEMENT->ChiInvS_vS;

	ICs = ELEMENT->ICs;
	ICc = ELEMENT->ICc;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGs = ELEMENT->I_vGs_vGs;
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
	I_vCs_vS  = ELEMENT->I_vCs_vS;
	I_vCs_vIs = ELEMENT->I_vCs_vIs;
	I_vCs_vIc = ELEMENT->I_vCs_vIc;
	I_vCc_vS  = ELEMENT->I_vCc_vS;
	I_vCc_vIs = ELEMENT->I_vCc_vIs;
	I_vCc_vIc = ELEMENT->I_vCc_vIc;

	Ihat_vS_vS  = ELEMENT->Ihat_vS_vS;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vIs = ELEMENT->D_vGs_vIs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vIc = ELEMENT->D_vGc_vIc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	ChiS_fS     = ELEMENT->ChiS_fS;
	ChiS_fIs    = ELEMENT->ChiS_fIs;
	ChiS_fIc    = ELEMENT->ChiS_fIc;
	ChiS_fIs_sp = ELEMENT->ChiS_fIs_sp;
	ChiS_fIc_sp = ELEMENT->ChiS_fIc_sp;

	GradChiS_fIs = ELEMENT->GradChiS_fIs;
	GradChiS_fIc = ELEMENT->GradChiS_fIc;

	I_vGs_fS  = ELEMENT->I_vGs_fS;
	I_vGs_fIs = ELEMENT->I_vGs_fIs;
	I_vGs_fIc = ELEMENT->I_vGs_fIc;
	I_vGc_fS  = ELEMENT->I_vGc_fS;
	I_vGc_fIs = ELEMENT->I_vGc_fIs;
	I_vGc_fIc = ELEMENT->I_vGc_fIc;
	I_vCs_fS  = ELEMENT->I_vCs_fS;
	I_vCs_fIs = ELEMENT->I_vCs_fIs;
	I_vCs_fIc = ELEMENT->I_vCs_fIc;
	I_vCc_fS  = ELEMENT->I_vCc_fS;
	I_vCc_fIs = ELEMENT->I_vCc_fIs;
	I_vCc_fIc = ELEMENT->I_vCc_fIc;

	D_vGs_fIs = ELEMENT->D_vGs_fIs;
	D_vGs_fIc = ELEMENT->D_vGs_fIc;
	D_vGc_fIs = ELEMENT->D_vGc_fIs;
	D_vGc_fIc = ELEMENT->D_vGc_fIc;

	Is_Weak_VV    = ELEMENT->Is_Weak_VV;
	Ic_Weak_VV    = ELEMENT->Ic_Weak_VV;
	Is_Weak_FF    = ELEMENT->Is_Weak_FF;
	Ic_Weak_FF    = ELEMENT->Ic_Weak_FF;
	Is_Weak_FF_sp = ELEMENT->Is_Weak_FF_sp;
	Ic_Weak_FF_sp = ELEMENT->Ic_Weak_FF_sp;

	Ds_Weak_VV = ELEMENT->Ds_Weak_VV;
	Dc_Weak_VV = ELEMENT->Dc_Weak_VV;

	// Allocate memory for arrays with multiple levels of dereferencing
	rst_vS        = malloc(NVREFMAX * sizeof *rst_vS);  // free
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

	GradChiGs_fIs = malloc(dE * sizeof *GradChiGs_fIs); // free
	GradChiGs_fIc = malloc(dE * sizeof *GradChiGs_fIc); // free
	GradChiGc_fIs = malloc(dE * sizeof *GradChiGc_fIs); // free
	GradChiGc_fIc = malloc(dE * sizeof *GradChiGc_fIc); // free

	// VOLUME Nodes (Order Independent)
	cubature(&rst_vGs,&dummyPtr_d,&dummyPtr_ui[0],&NvnGs[1],&dummy_ui,0,PGs,dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
	// Use E_rst_vC instead of rst_vGs (!= for PYR ELEMENTs)
	free(rst_vGs);
	E_rst_vC = get_rst_vC(ELEMENT); // free
	rst_vC = malloc((Nve+1)*dE * sizeof *rst_vC); // free (+1 for h-refined TET -> PYR) ToBeModified (remove +1)


	// Preliminary Operators
	IGs = identity_d(NvnGs[1]); // free

	ChiRefGs_vGs = basis(PGs,E_rst_vC,NvnGs[1],&Nbf,dE); // free

	if (strstr(BasisType,"Modal"))
		ChiGs_vGs = ChiRefGs_vGs;
	else if (strstr(BasisType,"Nodal"))
		ChiGs_vGs = IGs;

	ChiRefInvGs_vGs = inverse_d(NvnGs[1],NvnGs[1],ChiRefGs_vGs,IGs); // free
	ChiInvGs_vGs    = inverse_d(NvnGs[1],NvnGs[1],ChiGs_vGs,IGs);    // free
	TGs             = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGs[1],NvnGs[1],NvnGs[1],1.0,ChiRefInvGs_vGs,ChiGs_vGs); // free

	free(IGs);
	free(ChiRefGs_vGs);

	// Get Barycentric coordinates for lower dimensional ELEMENTs
	NFTypes = 1;
	BCoords_dEm1[0] = get_BCoords_dEm1(ELEMENT,0); // keep/free
	if (EType == WEDGE || EType == PYR) {
		NFTypes = 2;
		BCoords_dEm1[1] = get_BCoords_dEm1(ELEMENT,1);
	}

	for (IndFType = 0; IndFType < NFTypes; IndFType++)
		BCoords_F[IndFType] = malloc(sizeof *BCoords_F[IndFType]); // free


	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		B_Nve[IndFType]         = BCoords_dEm1[IndFType]->Nve;
		BCoords_F[IndFType]->S  = BCoords_dEm1[IndFType]->BCoords_S;
		BCoords_F[IndFType]->Is = BCoords_dEm1[IndFType]->BCoords_Is;
		BCoords_F[IndFType]->Ic = BCoords_dEm1[IndFType]->BCoords_Ic;

		// Store Nfn* in a manner consistent with Nvn*
		for (P = 0; P <= PMax; P++) {
			NfnS[P][IndFType]  = BCoords_dEm1[IndFType]->NfnS[P];
			NfnIs[P][IndFType] = BCoords_dEm1[IndFType]->NfnIs[P];
			NfnIc[P][IndFType] = BCoords_dEm1[IndFType]->NfnIc[P];
			w_fIs[P][IndFType] = BCoords_dEm1[IndFType]->w_fIs[P];
			w_fIc[P][IndFType] = BCoords_dEm1[IndFType]->w_fIc[P];
		}
		free(BCoords_dEm1[IndFType]->NfnS);
		free(BCoords_dEm1[IndFType]->NfnIs);
		free(BCoords_dEm1[IndFType]->NfnIc);
		free(BCoords_dEm1[IndFType]->w_fIs);
		free(BCoords_dEm1[IndFType]->w_fIc);
	}

	ones_Nf = malloc(Nf * sizeof *ones_Nf); // free
	for (f = 0; f < Nf; f++)
		ones_Nf[f] = 1;

	BCoords_V = malloc(NEhref * sizeof *BCoords_V); // free
	for (i = 0, iMax = NEhref; i < iMax; i++) {
		BCoords_V[i]     = malloc(     sizeof *BCoords_V[i]);       // free
		BCoords_V[i]->S  = calloc(NP , sizeof *(BCoords_V[i]->S));  // free
		BCoords_V[i]->Is = calloc(NP , sizeof *(BCoords_V[i]->Is)); // free
		BCoords_V[i]->Ic = calloc(NP , sizeof *(BCoords_V[i]->Ic)); // free
	}

	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {
		get_Pb_range(P,&PbMin,&PbMax);

		// Build preliminary operators needed for P adaptation.
		IS         = NULL;
		ChiRefS_vS = NULL;

		for (Pb = PbMax; Pb >= P; Pb--) {
			cubature(&rst_vS[0],&dummyPtr_d,&dummyPtr_ui[0],&NvnS[Pb],&dummy_ui,0,Pb,dE,NodeTypeS[Pb][Eclass]); free(dummyPtr_ui[0]); // free

			IS         = identity_d(NvnS[Pb]);  // free
			ChiRefS_vS = basis(Pb,rst_vS[0],NvnS[Pb],&Nbf,dE); // free
			if (strstr(BasisType,"Modal")) {
				ChiS_vS[P][Pb][0] = ChiRefS_vS;
			} else if (strstr(BasisType,"Nodal")) {
				ChiS_vS[P][Pb][0] = IS;
			}

			if (ChiInvS_vS[Pb][Pb][0] == NULL)
				ChiInvS_vS[Pb][Pb][0] = inverse_d(NvnS[Pb],NvnS[Pb],ChiS_vS[P][Pb][0],IS); // keep

			free(rst_vS[0]);
			if (Pb != P) {
				free(IS);
				free(ChiRefS_vS);
			} else {
				// Avoid negative Pb for P = 0
				break;
			}
		}
		cubature(&rst_vGc,&dummyPtr_d,&dummyPtr_ui[0],&NvnGc[P],&dummy_ui,0,PGc[P],        dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
		cubature(&rst_vCs,&dummyPtr_d,&dummyPtr_ui[0],&NvnCs[P],&dummy_ui,0,PCs[P][Eclass],dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free
		cubature(&rst_vCc,&dummyPtr_d,&dummyPtr_ui[0],&NvnCc[P],&dummy_ui,0,PCc[P][Eclass],dE,NodeTypeG[Eclass]); free(dummyPtr_ui[0]); // free

		// Preliminary Operators
		IGc          = identity_d(NvnGc[P]); // free
		ICs[P][P][0] = identity_d(NvnCs[P]); // keep
		ICc[P][P][0] = identity_d(NvnCc[P]); // keep

		ChiRefGc_vGc = basis(PGc[P],rst_vGc,NvnGc[P],&Nbf,dE);         // free
		ChiRefCs_vCs = basis(PCs[P][Eclass],rst_vCs,NvnCs[P],&Nbf,dE); // free
		ChiRefCc_vCc = basis(PCc[P][Eclass],rst_vCc,NvnCc[P],&Nbf,dE); // free

		if (strstr(BasisType,"Modal")) {
			ChiGc_vGc = ChiRefGc_vGc;
			ChiCs_vCs = ChiRefCs_vCs;
			ChiCc_vCc = ChiRefCc_vCc;
		} else if (strstr(BasisType,"Nodal")) {
			ChiGc_vGc = IGc;
			ChiCs_vCs = ICs[P][P][0];
			ChiCc_vCc = ICc[P][P][0];
		}

		ChiRefInvGc_vGc = inverse_d(NvnGc[P],NvnGc[P],ChiRefGc_vGc,IGc);          // free
		ChiRefInvCs_vCs = inverse_d(NvnCs[P],NvnCs[P],ChiRefCs_vCs,ICs[P][P][0]); // free
		ChiRefInvCc_vCc = inverse_d(NvnCc[P],NvnCc[P],ChiRefCc_vCc,ICc[P][P][0]); // free
		ChiRefInvS_vS   = inverse_d(NvnS[P], NvnS[P], ChiRefS_vS,IS);             // free

		ChiInvGc_vGc        = inverse_d(NvnGc[P],NvnGc[P],ChiGc_vGc,IGc);          // free
		ChiInvCs_vCs        = inverse_d(NvnCs[P],NvnCs[P],ChiCs_vCs,ICs[P][P][0]); // free
		ChiInvCc_vCc        = inverse_d(NvnCc[P],NvnCc[P],ChiCc_vCc,ICc[P][P][0]); // free

		TGc = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGc[P],NvnGc[P],1.0,ChiRefInvGc_vGc,ChiGc_vGc);      // free
		TCs = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,ChiRefInvCs_vCs,ChiCs_vCs);      // free
		TCc = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,ChiRefInvCc_vCc,ChiCc_vCc);      // free
		TS  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P], NvnS[P], NvnS[P], 1.0,ChiRefInvS_vS,ChiS_vS[P][P][0]); // free

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

		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			if (w_vIs[Pb])
				free(w_vIs[Pb]);
			cubature(&rst_vIs[0],&w_vIs[Pb],&dummyPtr_ui[0],&NvnIs[Pb],&dummy_ui,1,PIvs[Pb][Eclass],dE,NodeTypeIvs[Pb][Eclass]); free(dummyPtr_ui[0]); // free
			if (w_vIc[Pb])
				free(w_vIc[Pb]);
			cubature(&rst_vIc[0],&w_vIc[Pb],&dummyPtr_ui[0],&NvnIc[Pb],&dummy_ui,1,PIvc[Pb][Eclass],dE,NodeTypeIvc[Pb][Eclass]); free(dummyPtr_ui[0]); // free
			free(rst_vIs[0]);
			free(rst_vIc[0]);

			wInv_vIs = malloc(NvnIs[P] * sizeof *wInv_vIs); // free
			wInv_vIc = malloc(NvnIc[P] * sizeof *wInv_vIc); // free

			for (i = 0, iMax = NvnIs[P]; i < iMax; i++)
				wInv_vIs[i] = 1./w_vIs[P][i];

			for (i = 0, iMax = NvnIc[P]; i < iMax; i++)
				wInv_vIc[i] = 1./w_vIc[P][i];

			diag_w_vIs    = diag_d(w_vIs[Pb],NvnIs[Pb]); // free
			diag_w_vIc    = diag_d(w_vIc[Pb],NvnIc[Pb]); // free
			diag_wInv_vIs = diag_d(wInv_vIs, NvnIs[P]);  // free
			diag_wInv_vIc = diag_d(wInv_vIc, NvnIc[P]);  // free

			free(wInv_vIs);
			free(wInv_vIc);

			// VOLUME Operators
			for (i = iMax = NEhref; i--; ) {
				ELEMENT = get_ELEMENT_type(EType_h[i]);

				Eclass = get_Eclass(EType_h[i]);
				select_functions(&basis,&grad_basis,&cubature,EType_h[i]);

				NvnGs = ELEMENT->NvnGs;


				cubature(&rst_vS[0], &dummyPtr_d,&dummyPtr_ui[0],&NvnS[Pb], &dummy_ui,0,Pb,              dE,NodeTypeS[Pb][Eclass]);   free(dummyPtr_ui[0]); // free
				cubature(&rst_vIs[0],&dummyPtr_d,&dummyPtr_ui[0],&NvnIs[Pb],&dummy_ui,0,PIvs[Pb][Eclass],dE,NodeTypeIvs[Pb][Eclass]); free(dummyPtr_ui[0]); // free
				cubature(&rst_vIc[0],&dummyPtr_d,&dummyPtr_ui[0],&NvnIc[Pb],&dummy_ui,0,PIvc[Pb][Eclass],dE,NodeTypeIvc[Pb][Eclass]); free(dummyPtr_ui[0]); // free

				ChiRefGs_vS  = basis(PGs,rst_vS[0], NvnS[Pb], &Nbf,dE); // free
				ChiRefGs_vIs = basis(PGs,rst_vIs[0],NvnIs[Pb],&Nbf,dE); // free
				ChiRefGs_vIc = basis(PGs,rst_vIc[0],NvnIc[Pb],&Nbf,dE); // free

				free(E_rst_vC);
				free(ChiRefInvGs_vGs);

				E_rst_vC        = get_rst_vC(ELEMENT);                           // free
				IGs             = identity_d(NvnGs[1]);                          // free
				ChiRefGs_vGs    = basis(PGs,E_rst_vC,NvnGs[1],&Nbf,dE);          // free
				ChiRefInvGs_vGs = inverse_d(NvnGs[1],NvnGs[1],ChiRefGs_vGs,IGs); // free
				free(IGs);
				free(ChiRefGs_vGs);

				BCoords_V[i]->S[Pb]  = mm_Alloc_d(CBCM,CBT,CBT,NvnS[Pb], NvnGs[1],NvnGs[1],1.0,ChiRefGs_vS, ChiRefInvGs_vGs); // free
				BCoords_V[i]->Is[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIs[Pb],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vIs,ChiRefInvGs_vGs); // free
				BCoords_V[i]->Ic[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIc[Pb],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vIc,ChiRefInvGs_vGs); // free

				if (i) {
					free(rst_vS[0]);
					free(rst_vIs[0]);
					free(rst_vIc[0]);

					free(ChiRefGs_vIs);
					free(ChiRefGs_vIc);
				}

				free(ChiRefGs_vS);
			}

			// Note: Nvref >= NvrefSF
			for (vh = 0; vh < Nvref; vh++) {
				vrefSF = vh; // Used for operators defined only for setup_TP_operators
				mm_CTN_d(Nvve[vh],dE,Nve,VeV[vh],E_rst_vC,rst_vC);

				if (vh) {
					Indh = get_IndEhref(EType,vh);
					if (EType == PYR)
						ELEMENT_h = get_ELEMENT_type(EType_h[Indh]);
					else
						ELEMENT_h = ELEMENT;

					rst_vS[vh]  = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], dE,Nvve[vh],1.0,BCoords_V[Indh]->S[Pb], rst_vC); // free
					rst_vIs[vh] = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],dE,Nvve[vh],1.0,BCoords_V[Indh]->Is[Pb],rst_vC); // free
					rst_vIc[vh] = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],dE,Nvve[vh],1.0,BCoords_V[Indh]->Ic[Pb],rst_vC); // free

					ChiRefGs_vIs = basis(PGs,rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
					ChiRefGs_vIc = basis(PGs,rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				}

// Note: Nvve[vh] == ELEMENT_h->NvnGs[1]. Replace. (ToBeDeleted)
				ChiRefGs_vGs = basis(PGs,           rst_vC,     Nvve[vh], &Nbf,dE); // free
				ChiRefGs_vS  = basis(PGs,           rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefGc_vIs = basis(PGc[P],        rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefGc_vIc = basis(PGc[P],        rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				ChiRefGc_vS  = basis(PGc[P],        rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefCs_vS  = basis(PCs[P][Eclass],rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefCs_vIs = basis(PCs[P][Eclass],rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefCs_vIc = basis(PCs[P][Eclass],rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				ChiRefCc_vS  = basis(PCc[P][Eclass],rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free
				ChiRefCc_vIs = basis(PCc[P][Eclass],rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefCc_vIc = basis(PCc[P][Eclass],rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				ChiRefS_vIs  = basis(P,             rst_vIs[vh],ELEMENT_h->NvnIs[Pb],&Nbf,dE); // free
				ChiRefS_vIc  = basis(P,             rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE); // free
				ChiRefS_vS   = basis(P,             rst_vS[vh], ELEMENT_h->NvnS[Pb], &Nbf,dE); // free

				ChiGs_vGs = mm_Alloc_d(CBRM,CBNT,CBNT,Nvve[vh], NvnGs[1],NvnGs[1],1.0,ChiRefGs_vGs,TGs); // free
				ChiGs_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnGs[1],NvnGs[1],1.0,ChiRefGs_vS, TGs); // free
				ChiGs_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vIs,TGs); // free
				ChiGs_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vIc,TGs); // free
				ChiGc_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnGc[P],NvnGc[P],1.0,ChiRefGc_vS, TGc); // free
				ChiGc_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vIs,TGc); // free
				ChiGc_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vIc,TGc); // free
				ChiCs_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCs[P],NvnCs[P],1.0,ChiRefCs_vS, TCs); // free
				ChiCs_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIs,TCs); // free
				ChiCs_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCs[P],NvnCs[P],1.0,ChiRefCs_vIc,TCs); // free
				ChiCc_vS  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCc[P],NvnCc[P],1.0,ChiRefCc_vS, TCc); // free
				ChiCc_vIs = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIs,TCc); // free
				ChiCc_vIc = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCc[P],NvnCc[P],1.0,ChiRefCc_vIc,TCc); // free

				// Returned SF Operators
				ChiS_vIc[P][Pb][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIc,TS); // keep
				if (vrefSF < NvrefSF) {
					ChiS_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnS[P],NvnS[P],1.0,ChiRefS_vS,TS);  // keep
					ChiS_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIs,TS); // keep

					if (P == Pb) {
						I_vGs_vS[1][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb] ,NvnGs[1],NvnGs[1],1.0,ChiGs_vS, ChiInvGs_vGs); // keep
						I_vGs_vIs[1][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGs[1],NvnGs[1],1.0,ChiGs_vIs,ChiInvGs_vGs); // keep
						I_vGs_vIc[1][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGs[1],NvnGs[1],1.0,ChiGs_vIc,ChiInvGs_vGs); // keep
					}
					I_vGc_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnGc[P],NvnGc[P],1.0,ChiGc_vS, ChiInvGc_vGc); // keep
					I_vGc_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnGc[P],NvnGc[P],1.0,ChiGc_vIs,ChiInvGc_vGc); // keep
					I_vGc_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnGc[P],NvnGc[P],1.0,ChiGc_vIc,ChiInvGc_vGc); // keep
					I_vCs_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCs[P],NvnCs[P],1.0,ChiCs_vS, ChiInvCs_vCs); // keep
					I_vCs_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCs[P],NvnCs[P],1.0,ChiCs_vIs,ChiInvCs_vCs); // keep
					I_vCs_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCs[P],NvnCs[P],1.0,ChiCs_vIc,ChiInvCs_vCs); // keep
					I_vCc_vS[P][Pb][vrefSF]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnCc[P],NvnCc[P],1.0,ChiCc_vS, ChiInvCc_vCc); // keep
					I_vCc_vIs[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIs[Pb],NvnCc[P],NvnCc[P],1.0,ChiCc_vIs,ChiInvCc_vCc); // keep
					I_vCc_vIc[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnCc[P],NvnCc[P],1.0,ChiCc_vIc,ChiInvCc_vCc); // keep

					if (EFE) {
						// Incorrect operator for PYR -> TET (vrefSF != 0) but this operator is not used.
						Is_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIs[Pb],NvnIs[Pb], 1.0,ChiS_vIs[P][Pb][vrefSF],diag_w_vIs); // keep
						Ic_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIc[Pb],NvnIc[Pb], 1.0,ChiS_vIc[P][Pb][vrefSF],diag_w_vIc); // keep
					} else {
						printf("Error: Unsupported EFE.\n"), EXIT_MSG;
					}

					if (Collocated) {
						// Incorrect operator for PYR -> TET (vrefSF != 0) but this operator is not used.
						dummyPtr_d = Is_Weak_VV[P][Pb][vrefSF];
						Is_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnIs[Pb],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Ic_Weak_VV[P][Pb][vrefSF];
						Ic_Weak_VV[P][Pb][vrefSF] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnIc[Pb],NvnS[P],1.0,diag_wInv_vIc,dummyPtr_d); // keep
						free(dummyPtr_d);
					}
				}

				// Returned Adaptation Operators
				if (vh == 0 || P == Pb) {
					i = get_IndEhref(EType,vh);
					dummyPtr_d = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb], NvnS[P],NvnS[P],1.0,ChiRefS_vS,TS);  // free
					Ihat_vS_vS[P][Pb][vh]  = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnS[Pb],NvnS[P],ELEMENT_h->NvnS[Pb],
					                                    1.0,ELEMENT_h->ChiInvS_vS[Pb][Pb][0],dummyPtr_d); // keep
					free(dummyPtr_d);
					if (P == PGlobal && Pb == PGlobal) {
						I_vGs_vGs[1][1][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,Nvve[vh],NvnGs[1],NvnGs[1],1.0,ChiGs_vGs,ChiInvGs_vGs); // keep
					}
				}

				free(ChiRefGs_vGs);
				free(ChiRefGs_vS);
				free(ChiRefGs_vIs);
				free(ChiRefGs_vIc);
				free(ChiRefGc_vS);
				free(ChiRefGc_vIs);
				free(ChiRefGc_vIc);
				free(ChiRefCs_vS);
				free(ChiRefCs_vIs);
				free(ChiRefCs_vIc);
				free(ChiRefCc_vS);
				free(ChiRefCc_vIs);
				free(ChiRefCc_vIc);
				free(ChiRefS_vS);

				free(ChiGs_vGs);
				free(ChiGs_vS);
				free(ChiGs_vIs);
				free(ChiGs_vIc);
				free(ChiGc_vS);
				free(ChiGc_vIs);
				free(ChiGc_vIc);
				free(ChiCs_vS);
				free(ChiCs_vIs);
				free(ChiCs_vIc);
				free(ChiCc_vS);
				free(ChiCc_vIs);
				free(ChiCc_vIc);

				free(ChiRefS_vIs);
				free(ChiRefS_vIc);
			}

			plotting_element_info(&rst_vP,&dummyPtr_ui[0],&dummyPtr_ui[1],&NvnP,&dummy_ui,max(Pb,u1),EType); // free
			free(dummyPtr_ui[0]);
			free(dummyPtr_ui[1]);

			ChiRefGc_vP  = basis(PGc[P],        rst_vP,    NvnP,    &Nbf,dE); // free
			ChiRefS_vP   = basis(P,             rst_vP,    NvnP,    &Nbf,dE); // free

			ChiGc_vP            = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiRefGc_vP, TGc); // free

			// Returned Operators
			ChiS_vP[P][Pb][0]   = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnS[P], NvnS[P], 1.0,ChiRefS_vP,TS);           // keep
			I_vGc_vP[P][Pb][0]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGc[P],NvnGc[P],1.0,ChiGc_vP,  ChiInvGc_vGc); // keep

			free(ChiRefGc_vP);
			free(ChiRefS_vP);
			free(ChiGc_vP);

			// Operators (1/P -> P only)
			if (P == Pb) {
				ChiRefGs_vP  = basis(PGs,           rst_vP,    NvnP,    &Nbf,dE); // free
				ChiRefGs_vGc = basis(PGs,           rst_vGc,   NvnGc[P],&Nbf,dE); // free
				ChiRefGs_vCs = basis(PGs,           rst_vCs,   NvnCs[P],&Nbf,dE); // free
				ChiRefGc_vCc = basis(PGc[P],        rst_vCc,   NvnCc[P],&Nbf,dE); // free
				ChiRefS_vIs  = basis(P,             rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				ChiRefS_vIc  = basis(P,             rst_vIc[0],NvnIc[P],&Nbf,dE); // free

				ChiGs_vP          = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGs[1],NvnGs[1],1.0,ChiRefGs_vP, TGs); // free
				ChiGs_vGc         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vGc,TGs); // free
				ChiGs_vCs         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,ChiRefGs_vCs,TGs); // free
				ChiGc_vCc         = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiRefGc_vCc,TGc); // free

				GradChiRefGs_vCs = grad_basis(PGs,           rst_vCs,   NvnCs[P],&Nbf,dE); // free
				GradChiRefGs_vIs = grad_basis(PGs,           rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				GradChiRefGc_vCc = grad_basis(PGc[P],        rst_vCc,   NvnCc[P],&Nbf,dE); // free
				GradChiRefGc_vIc = grad_basis(PGc[P],        rst_vIc[0],NvnIc[P],&Nbf,dE); // free
				GradChiRefCs_vCs = grad_basis(PCs[P][Eclass],rst_vCs,   NvnCs[P],&Nbf,dE); // free
				GradChiRefCc_vCc = grad_basis(PCc[P][Eclass],rst_vCc,   NvnCc[P],&Nbf,dE); // free
				GradChiRefS_vIs  = grad_basis(P,             rst_vIs[0],NvnIs[P],&Nbf,dE); // free
				GradChiRefS_vIc  = grad_basis(P,             rst_vIc[0],NvnIc[P],&Nbf,dE); // free

				for (dim = 0; dim < dE; dim++) {
					GradChiGs_vCs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_vCs[dim],TGs); // free
					GradChiGs_vIs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_vIs[dim],TGs); // free
					GradChiGc_vCc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vCc[dim],TGc); // free
					GradChiGc_vIc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_vIc[dim],TGc); // free
					GradChiCs_vCs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiRefCs_vCs[dim],TCs); // free
					GradChiCc_vCc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiRefCc_vCc[dim],TCc); // free
					GradChiS_vIs[dim]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnS[P], NvnS[P], 1.0,GradChiRefS_vIs[dim],TS);   // free
					GradChiS_vIc[dim]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnS[P], NvnS[P], 1.0,GradChiRefS_vIc[dim],TS);   // free
				}

				// Returned Operators
				// VOLUME related operators
				I_vGs_vP[1][Pb][0]  = mm_Alloc_d(CBRM,CBNT,CBNT,NvnP,    NvnGs[1],NvnGs[1],1.0,ChiGs_vP, ChiInvGs_vGs); // keep
				I_vGs_vGc[1][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnGc[P],NvnGs[1],NvnGs[1],1.0,ChiGs_vGc,ChiInvGs_vGs); // keep
				I_vGs_vCs[1][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,ChiGs_vCs,ChiInvGs_vGs); // keep
				I_vGc_vCc[P][Pb][0] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,ChiGc_vCc,ChiInvGc_vGc); // keep

				for (dim = 0; dim < dE; dim++) {
					D_vGs_vCs[1][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnGs[1],NvnGs[1],1.0,GradChiGs_vCs[dim],ChiInvGs_vGs); // keep
					D_vGs_vIs[1][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnGs[1],NvnGs[1],1.0,GradChiGs_vIs[dim],ChiInvGs_vGs); // keep
					D_vGc_vCc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vCc[dim],ChiInvGc_vGc); // keep
					D_vGc_vIc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnGc[P],NvnGc[P],1.0,GradChiGc_vIc[dim],ChiInvGc_vGc); // keep
					D_vCs_vCs[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCs[P],NvnCs[P],NvnCs[P],1.0,GradChiCs_vCs[dim],ChiInvCs_vCs); // keep
					D_vCc_vCc[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnCc[P],NvnCc[P],NvnCc[P],1.0,GradChiCc_vCc[dim],ChiInvCc_vCc); // keep

					if (EFE) {
						Ds_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIs[P],NvnIs[P],1.0,GradChiS_vIs[dim],diag_w_vIs); // keep
						Dc_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NvnIc[P],NvnIc[P],1.0,GradChiS_vIc[dim],diag_w_vIc); // keep
					} else {
						printf("Error: Unsupported EFE.\n"), EXIT_MSG;
					}
					if (Collocated) {
						dummyPtr_d = Ds_Weak_VV[P][Pb][0][dim];
						Ds_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NvnIs[Pb],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Dc_Weak_VV[P][Pb][0][dim];
						Dc_Weak_VV[P][Pb][0][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NvnIc[Pb],NvnS[P],1.0,diag_wInv_vIc,dummyPtr_d); // keep
						free(dummyPtr_d);
					}
				}

				free(ChiRefGs_vP);
				free(ChiRefGs_vGc);
				free(ChiRefGs_vCs);
				free(ChiRefGc_vCc);
				free(ChiRefS_vIs);
				free(ChiRefS_vIc);

				free(ChiGs_vP);
				free(ChiGs_vGc);
				free(ChiGs_vCs);
				free(ChiGc_vCc);

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
			for (f = 0; f < Nf; f++) {
				IndFType = get_IndFType(Eclass,f);

				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;

					mm_CTN_d(Nfve[f],dE,Nve,VeF[Vf],E_rst_vC,rst_vC);

					rst_fS  = mm_Alloc_d(CBCM,CBNT,CBNT,NfnS[Pb][IndFType], dE,B_Nve[IndFType],1.0,BCoords_F[IndFType]->S[Pb], rst_vC); // free
					rst_fIs = mm_Alloc_d(CBCM,CBNT,CBNT,NfnIs[Pb][IndFType],dE,B_Nve[IndFType],1.0,BCoords_F[IndFType]->Is[Pb],rst_vC); // free
					rst_fIc = mm_Alloc_d(CBCM,CBNT,CBNT,NfnIc[Pb][IndFType],dE,B_Nve[IndFType],1.0,BCoords_F[IndFType]->Ic[Pb],rst_vC); // free

					diag_w_fIs = diag_d(w_fIs[Pb][IndFType],NfnIs[Pb][IndFType]); // free
					diag_w_fIc = diag_d(w_fIc[Pb][IndFType],NfnIc[Pb][IndFType]); // free

					ChiRefS_fS  = basis(P,rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefS_fIs = basis(P,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefS_fIc = basis(P,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					GradChiRefS_fIs = grad_basis(P,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					GradChiRefS_fIc = grad_basis(P,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					// Returned Operators
					ChiS_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnS[P],NvnS[P],1.0,ChiRefS_fS, TS); // keep
					ChiS_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnS[P],NvnS[P],1.0,ChiRefS_fIs,TS); // keep
					ChiS_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnS[P],NvnS[P],1.0,ChiRefS_fIc,TS); // keep

					for (dim = 0; dim < dE; dim++) {
						GradChiS_fIs[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnS[P],NvnS[P],1.0,GradChiRefS_fIs[dim],TS); // keep
						GradChiS_fIc[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnS[P],NvnS[P],1.0,GradChiRefS_fIc[dim],TS); // keep
					}

					Is_Weak_FF[P][Pb][Vf] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NfnIs[Pb][IndFType],NfnIs[Pb][IndFType],-1.0,ChiS_fIs[P][Pb][Vf],diag_w_fIs); // keep
					Ic_Weak_FF[P][Pb][Vf] = mm_Alloc_d(CBRM,CBT,CBNT,NvnS[P],NfnIc[Pb][IndFType],NfnIc[Pb][IndFType],-1.0,ChiS_fIc[P][Pb][Vf],diag_w_fIc); // keep

					if (Collocated) {
						dummyPtr_d = Is_Weak_FF[P][Pb][Vf];
						Is_Weak_FF[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIs[P],NfnIs[Pb][IndFType],NvnS[P],1.0,diag_wInv_vIs,dummyPtr_d); // keep
						free(dummyPtr_d);
						dummyPtr_d = Ic_Weak_FF[P][Pb][Vf];
						Ic_Weak_FF[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnIc[P],NfnIc[Pb][IndFType],NvnS[P],1.0,diag_wInv_vIc,dummyPtr_d); // keep
						free(dummyPtr_d);
					}

					if (VFPartUnity[Eclass]) {
						convert_to_CSR_d(NfnIs[Pb][IndFType],NvnS[P],ChiS_fIs[P][Pb][Vf],&ChiS_fIs_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NfnIc[Pb][IndFType],NvnS[P],ChiS_fIc[P][Pb][Vf],&ChiS_fIc_sp[P][Pb][Vf]); // keep

						convert_to_CSR_d(NvnS[P],NfnIs[Pb][IndFType],Is_Weak_FF[P][Pb][Vf],&Is_Weak_FF_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NvnS[P],NfnIc[Pb][IndFType],Ic_Weak_FF[P][Pb][Vf],&Ic_Weak_FF_sp[P][Pb][Vf]); // keep
					}

					ChiRefGs_fS  = basis(PGs           ,rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefGs_fIs = basis(PGs           ,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefGs_fIc = basis(PGs           ,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
					ChiRefGc_fS  = basis(PGc[P]        ,rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefGc_fIs = basis(PGc[P]        ,rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefGc_fIc = basis(PGc[P]        ,rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
					ChiRefCs_fS  = basis(PCs[P][Eclass],rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefCs_fIs = basis(PCs[P][Eclass],rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefCs_fIc = basis(PCs[P][Eclass],rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free
					ChiRefCc_fS  = basis(PCc[P][Eclass],rst_fS, NfnS[Pb][IndFType], &Nbf,dE); // free
					ChiRefCc_fIs = basis(PCc[P][Eclass],rst_fIs,NfnIs[Pb][IndFType],&Nbf,dE); // free
					ChiRefCc_fIc = basis(PCc[P][Eclass],rst_fIc,NfnIc[Pb][IndFType],&Nbf,dE); // free

					GradChiRefGs_fIs = grad_basis(PGs,   rst_fIs,NfnIs[P][IndFType],&Nbf,dE); // free
					GradChiRefGs_fIc = grad_basis(PGs,   rst_fIc,NfnIc[P][IndFType],&Nbf,dE); // free
					GradChiRefGc_fIs = grad_basis(PGc[P],rst_fIs,NfnIs[P][IndFType],&Nbf,dE); // free
					GradChiRefGc_fIc = grad_basis(PGc[P],rst_fIc,NfnIc[P][IndFType],&Nbf,dE); // free

					ChiGs_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGs[1],NvnGs[1],1.0,ChiRefGs_fS, TGs); // free
					ChiGs_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiRefGs_fIs,TGs); // free
					ChiGs_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiRefGs_fIc,TGs); // free
					ChiGc_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGc[P],NvnGc[P],1.0,ChiRefGc_fS, TGc); // free
					ChiGc_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIs,TGc); // free
					ChiGc_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiRefGc_fIc,TGc); // free
					ChiCs_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCs[P],NvnCs[P],1.0,ChiRefCs_fS, TCs); // free
					ChiCs_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIs,TCs); // free
					ChiCs_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiRefCs_fIc,TCs); // free
					ChiCc_fS  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCc[P],NvnCc[P],1.0,ChiRefCc_fS, TCc); // free
					ChiCc_fIs = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIs,TCc); // free
					ChiCc_fIc = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiRefCc_fIc,TCc); // free

					for (dim = 0; dim < dE; dim++) {
						GradChiGs_fIs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_fIs[dim],TGs); // free
						GradChiGs_fIc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiRefGs_fIc[dim],TGs); // free
						GradChiGc_fIs[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_fIs[dim],TGc); // free
						GradChiGc_fIc[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiRefGc_fIc[dim],TGc); // free
					}

					// Returned Operators
					if (P == Pb) {
						// Outside of fh == 0 if condition as operator is used for testing h-adaptation connectivity.
						I_vGs_fS[1][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGs[1],NvnGs[1],1.0,ChiGs_fS, ChiInvGs_vGs); // keep
					}

					if (fh == 0) {
						if (P == Pb) {
							I_vGs_fIs[1][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiGs_fIs,ChiInvGs_vGs); // keep
							I_vGs_fIc[1][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,ChiGs_fIc,ChiInvGs_vGs); // keep

							for (dim = 0; dim < dE; dim++) {
								D_vGs_fIs[1][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiGs_fIs[dim],ChiInvGs_vGs); // keep
								D_vGs_fIc[1][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGs[1],NvnGs[1],1.0,GradChiGs_fIc[dim],ChiInvGs_vGs); // keep
							}
						}
						I_vGc_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnGc[P],NvnGc[P],1.0,ChiGc_fS, ChiInvGc_vGc); // keep
						I_vGc_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIs,ChiInvGc_vGc); // keep
						I_vGc_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,ChiGc_fIc,ChiInvGc_vGc); // keep
						I_vCs_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCs[P],NvnCs[P],1.0,ChiCs_fS, ChiInvCs_vCs); // keep
						I_vCs_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIs,ChiInvCs_vCs); // keep
						I_vCs_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCs[P],NvnCs[P],1.0,ChiCs_fIc,ChiInvCs_vCs); // keep
						I_vCc_fS[P][Pb][Vf]  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnS[Pb][IndFType], NvnCc[P],NvnCc[P],1.0,ChiCc_fS, ChiInvCc_vCc); // keep
						I_vCc_fIs[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIs,ChiInvCc_vCc); // keep
						I_vCc_fIc[P][Pb][Vf] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnCc[P],NvnCc[P],1.0,ChiCc_fIc,ChiInvCc_vCc); // keep

						for (dim = 0; dim < dE; dim++) {
							D_vGc_fIs[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiGc_fIs[dim],ChiInvGc_vGc); // keep
							D_vGc_fIc[P][Pb][Vf][dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc[Pb][IndFType],NvnGc[P],NvnGc[P],1.0,GradChiGc_fIc[dim],ChiInvGc_vGc); // keep
						}
					}

					free(ChiRefGs_fS);
					free(ChiRefGs_fIs);
					free(ChiRefGs_fIc);
					free(ChiRefGc_fS);
					free(ChiRefGc_fIs);
					free(ChiRefGc_fIc);
					free(ChiRefCs_fS);
					free(ChiRefCs_fIs);
					free(ChiRefCs_fIc);
					free(ChiRefCc_fS);
					free(ChiRefCc_fIs);
					free(ChiRefCc_fIc);

					free(ChiGs_fS);
					free(ChiGs_fIs);
					free(ChiGs_fIc);
					free(ChiGc_fS);
					free(ChiGc_fIs);
					free(ChiGc_fIc);
					free(ChiCs_fS);
					free(ChiCs_fIs);
					free(ChiCs_fIc);
					free(ChiCc_fS);
					free(ChiCc_fIs);
					free(ChiCc_fIc);

					array_free2_d(dE,GradChiRefGs_fIs);
					array_free2_d(dE,GradChiRefGs_fIc);
					array_free2_d(dE,GradChiRefGc_fIs);
					array_free2_d(dE,GradChiRefGc_fIc);

					for (dim = 0; dim < dE; dim++) {
						free(GradChiGs_fIs[dim]);
						free(GradChiGs_fIc[dim]);
						free(GradChiGc_fIs[dim]);
						free(GradChiGc_fIc[dim]);
					}

					free(rst_fS);
					free(rst_fIs);
					free(rst_fIc);

					free(diag_w_fIs);
					free(diag_w_fIc);

					free(ChiRefS_fS);
					free(ChiRefS_fIs);
					free(ChiRefS_fIc);

					array_free2_d(dE,GradChiRefS_fIs);
					array_free2_d(dE,GradChiRefS_fIc);
				}
			}
			for (vref = 0; vref < Nvref; vref++) {
				free(rst_vS[vref]);
				free(rst_vIs[vref]);
				free(rst_vIc[vref]);
			}
			for (i = 0; i < NEhref; i++) {
				free(BCoords_V[i]->S[Pb]);
				free(BCoords_V[i]->Is[Pb]);
				free(BCoords_V[i]->Ic[Pb]);
			}

			free(rst_vP);

			free(diag_w_vIs);
			free(diag_w_vIc);
			free(diag_wInv_vIs);
			free(diag_wInv_vIc);
		}
		free(rst_vGc);
		free(rst_vCs);
		free(rst_vCc);

		free(ChiInvGc_vGc);
		free(ChiInvCs_vCs);
		free(ChiInvCc_vCc);

		free(TGc);
		free(TCs);
		free(TCc);
		free(TS);
	}
	free(rst_vS);
	free(rst_vIs);
	free(rst_vIc);

	free(ChiRefInvGs_vGs);

	for (IndFType = 0; IndFType < NFTypes; IndFType++) {
		array_free2_d(PMax+1,BCoords_F[IndFType]->S);
		array_free2_d(PMax+1,BCoords_F[IndFType]->Is);
		array_free2_d(PMax+1,BCoords_F[IndFType]->Ic);
		free(BCoords_F[IndFType]);
		free(BCoords_dEm1[IndFType]);
	}

	array_free3_d(NP,NESUBCMAX,w_fIs);
	array_free3_d(NP,NESUBCMAX,w_fIc);

	for (i = 0; i < NEhref; i++) {
		free(BCoords_V[i]->S);
		free(BCoords_V[i]->Is);
		free(BCoords_V[i]->Ic);
		free(BCoords_V[i]);
	}
	free(BCoords_V);

	free(E_rst_vC);
	free(rst_vC);

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

	free(GradChiGs_fIs);
	free(GradChiGs_fIc);
	free(GradChiGc_fIs);
	free(GradChiGc_fIc);

	free(ones_Nf);
}

static void setup_ELEMENT_VeF(const unsigned int EType)
{
	/*
	 *	Comments:
	 *		While VeV defines VeVref_LINE and VeVref_TRI, the redundant definitions are included here as the ELEMENTs
	 *		are set up in order and VeV may not yet be defined.
	 */

	// Standard datatypes
	unsigned int i, j, k, l, iMax, jMax, kMax, lMax, iStep,
	             Nve, *Nfve, *Nfref, Nf, *VeFcon;
	double *VeF;

	struct S_ELEMENT *ELEMENT;

	// silence
	VeF = NULL;

	ELEMENT = get_ELEMENT_type(EType);
	Nve    = ELEMENT->Nve;
	Nfve   = ELEMENT->Nfve;
	Nf     = ELEMENT->Nf;
	VeFcon = ELEMENT->VeFcon;

	Nfref   = ELEMENT->Nfref;

	unsigned int size_VeF = 0;
	double VeVref_LINE[12]  = {1.0 , 0.0 , // 0
	                           0.0 , 1.0 ,
	                           1.0 , 0.0 , // 1
	                           0.5 , 0.5 ,
	                           0.5 , 0.5 , // 2
	                           0.0 , 1.0 },
	       VeVref_TRI[45]   = {1.0 , 0.0 , 0.0 , // 0
	                           0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , // 1
	                           0.5 , 0.5 , 0.0 ,
	                           0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , // 2
	                           0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , // 3
	                           0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 1.0 ,
	                           0.0 , 0.5 , 0.5 , // 4
	                           0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.5 , 0.0 },
	       VeVref_QUAD[144] = {1.0 , 0.0 , 0.0 , 0.0 , // 0
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 , // 1
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25,
	                           0.5 , 0.5 , 0.0 , 0.0 , // 2
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.25, 0.25, 0.25, 0.25,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 , // 3
	                           0.25, 0.25, 0.25, 0.25,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.25, 0.25, 0.25, 0.25, // 4
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 , // 5
	                           0.5 , 0.5 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.5 , 0.5 , 0.0 , 0.0 , // 6
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.0 , 0.0 , 0.5 , 0.5 ,
	                           0.0 , 0.0 , 0.0 , 1.0 ,
	                           1.0 , 0.0 , 0.0 , 0.0 , // 7
	                           0.0 , 1.0 , 0.0 , 0.0 ,
	                           0.5 , 0.0 , 0.5 , 0.0 ,
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.5 , 0.0 , 0.5 , 0.0 , // 8
	                           0.0 , 0.5 , 0.0 , 0.5 ,
	                           0.0 , 0.0 , 1.0 , 0.0 ,
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
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref_LINE[j*kMax+k];
		}}}

		break;
	case TET:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref_TRI[j*kMax+k];
		}}}

		break;
	case HEX:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref_QUAD[j*kMax+k];
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
					VeF[iStep+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref_TRI[j*kMax+k];
				else if (Nfve[i] == 4)
					VeF[iStep+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref_QUAD[j*kMax+k];
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
	free(VeF);
}

static void setup_TP_operators(const unsigned int EType)
{
	/*
	 *	Purpose:
	 *		Compute operators for elements which are tensor-products of lower dimensional elements.
	 *
	 *	Comments:
	 *		This functions works by taking tensor-products of the operators of lower-dimensional ELEMENT operators. This
	 *		allows for a reduction in the computational complexity of computing these operators and is very intuitive
	 *		once the concept is understood.
	 *
	 *	References:
	 *		Thesis (ToBeModified)
	 */

	// Returned operators
	unsigned int *NvnGs, *NvnGc, *NvnCs, *NvnCc, *NvnIs, *NvnIc, *NvnS, **NfnS, **NfnIs, **NfnIc;
	double       **w_vIs, **w_vIc,
	             ****ChiS_vP, ****ChiS_vS, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS,
	             ****I_vGs_vP, ****I_vGs_vGs, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vS, ****I_vGs_vIs,
	             ****I_vGc_vP,                ****I_vGc_vCc,                ****I_vGc_vS, ****I_vGc_vIc,
	             ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vIs, ****I_vCc_vIc,
	             ****Ihat_vS_vS,
	             ****L2hat_vS_vS,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fS, ****ChiS_fIs, ****ChiS_fIc,
	             ****I_vGs_fS, ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fS, ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vCs_fS, ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fS, ****I_vCc_fIs, ****I_vCc_fIc,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FF, ****Ic_Weak_FF,
	             *****Ds_Weak_VV, *****Dc_Weak_VV;

	struct S_OpCSR ****ChiS_fIs_sp, ****ChiS_fIc_sp,
	               *****Ds_Weak_VV_sp, *****Dc_Weak_VV_sp,
	               ****Is_Weak_FF_sp, ****Ic_Weak_FF_sp;

	// Initialize DB Parameters
	unsigned int EFE          = DB.EFE,
	             PGlobal      = DB.PGlobal,
	             Collocated   = DB.Collocated,
	             *VFPartUnity = DB.VFPartUnity;

	// Standard datatypes
	unsigned int dim, P, vh, f, fh, Pb, PSMin, PSMax, PbMin, PbMax, fhMax, Nvref, IndClass,
	             Eclass, dE, Nf, Vf, *Nfref, *ones_Nf,
	             NIn[3], NOut[3];
	double       *OP[3];

	struct S_ELEMENT *ELEMENT, *ELEMENTclass[2];

	// silence
	ELEMENTclass[1] = NULL;

	Eclass = get_Eclass(EType);
	ELEMENT = get_ELEMENT_type(EType);

	ELEMENTclass[0] = ELEMENT->ELEMENTclass[0];
	if (EType == WEDGE)
		ELEMENTclass[1] = ELEMENT->ELEMENTclass[1];

	dE    = ELEMENT->d;
	Nf    = ELEMENT->Nf;
	Nfref = ELEMENT->Nfref;
	Nvref = ELEMENT->Nvref;

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
	NfnS  = ELEMENT->NfnS;
	NfnIs = ELEMENT->NfnIs;
	NfnIc = ELEMENT->NfnIc;

	w_vIs = ELEMENT->w_vIs;
	w_vIc = ELEMENT->w_vIc;

	ChiS_vP    = ELEMENT->ChiS_vP;
	ChiS_vS    = ELEMENT->ChiS_vS;
	ChiS_vIs   = ELEMENT->ChiS_vIs;
	ChiS_vIc   = ELEMENT->ChiS_vIc;
	ChiInvS_vS = ELEMENT->ChiInvS_vS;

	I_vGs_vP  = ELEMENT->I_vGs_vP;
	I_vGs_vGs = ELEMENT->I_vGs_vGs;
	I_vGs_vGc = ELEMENT->I_vGs_vGc;
	I_vGs_vCs = ELEMENT->I_vGs_vCs;
	I_vGs_vS  = ELEMENT->I_vGs_vS;
	I_vGs_vIs = ELEMENT->I_vGs_vIs;
	I_vGc_vP  = ELEMENT->I_vGc_vP;
	I_vGc_vCc = ELEMENT->I_vGc_vCc;
	I_vGc_vS  = ELEMENT->I_vGc_vS;
	I_vGc_vIc = ELEMENT->I_vGc_vIc;
	I_vCs_vIs = ELEMENT->I_vCs_vIs;
	I_vCs_vIc = ELEMENT->I_vCs_vIc;
	I_vCc_vIs = ELEMENT->I_vCc_vIs;
	I_vCc_vIc = ELEMENT->I_vCc_vIc;

	Ihat_vS_vS  = ELEMENT->Ihat_vS_vS;
	L2hat_vS_vS = ELEMENT->L2hat_vS_vS;

	D_vGs_vCs = ELEMENT->D_vGs_vCs;
	D_vGs_vIs = ELEMENT->D_vGs_vIs;
	D_vGc_vCc = ELEMENT->D_vGc_vCc;
	D_vGc_vIc = ELEMENT->D_vGc_vIc;
	D_vCs_vCs = ELEMENT->D_vCs_vCs;
	D_vCc_vCc = ELEMENT->D_vCc_vCc;

	ChiS_fS     = ELEMENT->ChiS_fS;
	ChiS_fIs    = ELEMENT->ChiS_fIs;
	ChiS_fIc    = ELEMENT->ChiS_fIc;
	ChiS_fIs_sp = ELEMENT->ChiS_fIs_sp;
	ChiS_fIc_sp = ELEMENT->ChiS_fIc_sp;

	I_vGs_fS  = ELEMENT->I_vGs_fS;
	I_vGs_fIs = ELEMENT->I_vGs_fIs;
	I_vGs_fIc = ELEMENT->I_vGs_fIc;
	I_vGc_fS  = ELEMENT->I_vGc_fS;
	I_vGc_fIs = ELEMENT->I_vGc_fIs;
	I_vGc_fIc = ELEMENT->I_vGc_fIc;
	I_vCs_fS  = ELEMENT->I_vCs_fS;
	I_vCs_fIs = ELEMENT->I_vCs_fIs;
	I_vCs_fIc = ELEMENT->I_vCs_fIc;
	I_vCc_fS  = ELEMENT->I_vCc_fS;
	I_vCc_fIs = ELEMENT->I_vCc_fIs;
	I_vCc_fIc = ELEMENT->I_vCc_fIc;

	Is_Weak_VV    = ELEMENT->Is_Weak_VV;
	Ic_Weak_VV    = ELEMENT->Ic_Weak_VV;
	Is_Weak_FF    = ELEMENT->Is_Weak_FF;
	Ic_Weak_FF    = ELEMENT->Ic_Weak_FF;
	Is_Weak_FF_sp = ELEMENT->Is_Weak_FF_sp;
	Ic_Weak_FF_sp = ELEMENT->Ic_Weak_FF_sp;

	Ds_Weak_VV    = ELEMENT->Ds_Weak_VV;
	Dc_Weak_VV    = ELEMENT->Dc_Weak_VV;
	Ds_Weak_VV_sp = ELEMENT->Ds_Weak_VV_sp;
	Dc_Weak_VV_sp = ELEMENT->Dc_Weak_VV_sp;

	get_PS_range(&PSMin,&PSMax);
	if (Eclass == C_TP) {
		NvnGs[1] = pow(ELEMENTclass[0]->NvnGs[1],dE);

		for (P = PSMin; P <= PSMax; P++) {
			NvnGc[P] = pow(ELEMENTclass[0]->NvnGc[P],dE);
			NvnCs[P] = pow(ELEMENTclass[0]->NvnCs[P],dE);
			NvnCc[P] = pow(ELEMENTclass[0]->NvnCc[P],dE);
			NvnS[P]  = pow(ELEMENTclass[0]->NvnS[P],dE);

			get_Pb_range(P,&PbMin,&PbMax);
			for (Pb = PbMin; Pb <= PbMax; Pb++) {
				NvnIs[Pb] = pow(ELEMENTclass[0]->NvnIs[Pb],dE);
				NvnIc[Pb] = pow(ELEMENTclass[0]->NvnIc[Pb],dE);

				NfnS[Pb][0]  = pow(ELEMENTclass[0]->NvnS[Pb], dE-1)*(ELEMENTclass[0]->NfnS[Pb][0]);
				NfnIs[Pb][0] = pow(ELEMENTclass[0]->NvnIs[Pb],dE-1)*(ELEMENTclass[0]->NfnIs[Pb][0]);
				NfnIc[Pb][0] = pow(ELEMENTclass[0]->NvnIc[Pb],dE-1)*(ELEMENTclass[0]->NfnIc[Pb][0]);

				get_sf_parameters(1,ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->w_vIs[Pb],0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				if (w_vIs[Pb] == NULL)
					w_vIs[Pb] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->w_vIc[Pb],0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				if (w_vIc[Pb] == NULL)
					w_vIc[Pb] = sf_assemble_d(NIn,NOut,dE,OP);

				for (vh = 0; vh < Nvref; vh++) {
					if (vh == 0 || P == Pb) {
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->Ihat_vS_vS[P][Pb],
						                   0,0,NULL,NIn,NOut,OP,dE,vh,Eclass);
						Ihat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->L2hat_vS_vS[P][Pb],
						                   0,0,NULL,NIn,NOut,OP,dE,vh,Eclass);
						L2hat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (P == PGlobal && Pb == PGlobal) {
							get_sf_parametersV(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->I_vGs_vGs[1][1],
							                   0,0,NULL,NIn,NOut,OP,dE,vh,Eclass);
							I_vGs_vGs[1][1][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
					}
				}

				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->ChiS_vP[P][Pb][0],
				                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				ChiS_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->ChiS_vIc[P][Pb][0],
				                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				ChiS_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->I_vGc_vP[P][Pb][0],
				                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
				I_vGc_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

				// Note: Most VOLUME operators need not interpolate between different orders
				if (P == Pb) {
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->ChiS_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnP[P],ELEMENTclass[0]->I_vGs_vP[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vP[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->ChiInvS_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiInvS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->ChiS_vIs[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGc[Pb],ELEMENTclass[0]->I_vGs_vGc[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vGc[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vCs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGs_vS[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vS[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vIs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vCc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGc_vS[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],
					                  0,0,NULL,NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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

					for (dim = 0; dim < dE; dim++) {
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],
						                  ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->D_vGs_vCs[1][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGs_vCs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],
						                  ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->D_vGs_vIs[1][Pb][0][0],
						                  NIn,NOut,OP,dE,dim,Eclass);
						D_vGs_vIs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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
							get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
							                  ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ds_Weak_VV[P][Pb][0][0],
							                  NIn,NOut,OP,dE,dim,Eclass);
							Ds_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
							                  ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Dc_Weak_VV[P][Pb][0][0],
							                  NIn,NOut,OP,dE,dim,Eclass);
							Dc_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							if (Collocated) {
								convert_to_CSR_d(NvnS[P],NvnIs[Pb],Ds_Weak_VV[P][Pb][0][dim],&Ds_Weak_VV_sp[P][Pb][0][dim]); // keep
								convert_to_CSR_d(NvnS[P],NvnIc[Pb],Dc_Weak_VV[P][Pb][0][dim],&Dc_Weak_VV_sp[P][Pb][0][dim]); // keep
							}
						} else {
							printf("Error: Unsupported EFE.\n"), EXIT_MSG;
						}
					}
				}

				for (f = 0; f < Nf; f++) {
				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],   ELEMENTclass[0]->ChiS_vS[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnS[Pb][0],ELEMENTclass[0]->ChiS_fS[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->ChiS_vIs[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->ChiS_fIs[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->ChiS_vIc[P][Pb],
					                   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->ChiS_fIc[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (Collocated || VFPartUnity[Eclass]) {
						convert_to_CSR_d(NfnIs[Pb][0],NvnS[P],ChiS_fIs[P][Pb][Vf],&ChiS_fIs_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NfnIc[Pb][0],NvnS[P],ChiS_fIc[P][Pb][Vf],&ChiS_fIc_sp[P][Pb][Vf]); // keep
					}

					if (P == Pb) {
						// Outside of fh == 0 if condition as operator is used for testing h-adaptation connectivity.
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnS[Pb],   ELEMENTclass[0]->I_vGs_vS[1][Pb],
						                   ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnS[Pb][0],ELEMENTclass[0]->I_vGs_fS[1][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fS[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}

					get_sf_parametersF(ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb],
					                   ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_FF[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					Is_Weak_FF[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parametersF(ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb],
					                   ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_FF[P][Pb],
					                   NIn,NOut,OP,dE,Vf,Eclass);
					Ic_Weak_FF[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (Collocated || VFPartUnity[Eclass]) {
						convert_to_CSR_d(NvnS[P],NfnIs[Pb][0],Is_Weak_FF[P][Pb][Vf],&Is_Weak_FF_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NvnS[P],NfnIc[Pb][0],Ic_Weak_FF[P][Pb][Vf],&Ic_Weak_FF_sp[P][Pb][Vf]); // keep
					}

					if (fh == 0) {
						if (P == Pb) {
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->I_vGs_vIs[1][Pb],
							                   ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->I_vGs_fIs[1][Pb],
							                   NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIs[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->I_vGs_vIc[1][Pb],
							                   ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->I_vGs_fIc[1][Pb],
							                   NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIc[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[Pb],   ELEMENTclass[0]->I_vGc_vS[P][Pb],
						                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnS[Pb][0],ELEMENTclass[0]->I_vGc_fS[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->I_vGc_vIs[P][Pb],
						                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->I_vGc_fIs[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->I_vGc_vIc[P][Pb],
						                   ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->I_vGc_fIc[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnS[Pb],   ELEMENTclass[0]->I_vCs_vS[P][Pb],
						                   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnS[Pb][0],ELEMENTclass[0]->I_vCs_fS[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIs[Pb],   ELEMENTclass[0]->I_vCs_vIs[P][Pb],
						                   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnIs[Pb][0],ELEMENTclass[0]->I_vCs_fIs[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIs[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NvnIc[Pb],   ELEMENTclass[0]->I_vCs_vIc[P][Pb],
						                   ELEMENTclass[0]->NvnCs[P],ELEMENTclass[0]->NfnIc[Pb][0],ELEMENTclass[0]->I_vCs_fIc[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fIc[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NvnS[Pb],   ELEMENTclass[0]->I_vCc_vS[P][Pb],
						                   ELEMENTclass[0]->NvnCc[P],ELEMENTclass[0]->NfnS[Pb][0],ELEMENTclass[0]->I_vCc_fS[P][Pb],
						                   NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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
		unsigned int NIn0, NIn1, NOut0, NOut1;
		double       *OP0, *OP1, **OPF0, **OPF1;

		NvnGs[1] = (ELEMENTclass[0]->NvnGs[1])*(ELEMENTclass[1]->NvnGs[1]);

		for (P = PSMin; P <= PSMax; P++) {
			NvnGc[P] = (ELEMENTclass[0]->NvnGc[P])*(ELEMENTclass[1]->NvnGc[P]);
			NvnCs[P] = (ELEMENTclass[0]->NvnCs[P])*(ELEMENTclass[1]->NvnCs[P]);
			NvnCc[P] = (ELEMENTclass[0]->NvnCc[P])*(ELEMENTclass[1]->NvnCc[P]);
			NvnS[P]  = (ELEMENTclass[0]->NvnS[P])*(ELEMENTclass[1]->NvnS[P]);

			get_Pb_range(P,&PbMin,&PbMax);
			for (Pb = PbMin; Pb <= PbMax; Pb++) {
				NvnIs[Pb] = (ELEMENTclass[0]->NvnIs[Pb])*(ELEMENTclass[1]->NvnIs[Pb]);
				NvnIc[Pb] = (ELEMENTclass[0]->NvnIc[Pb])*(ELEMENTclass[1]->NvnIc[Pb]);

				NfnS[Pb][0]  = (ELEMENTclass[0]->NfnS[Pb][0])* (ELEMENTclass[1]->NvnS[Pb]);
				NfnS[Pb][1]  = (ELEMENTclass[0]->NvnS[Pb])*    (ELEMENTclass[1]->NfnS[Pb][0]);
				NfnIs[Pb][0] = (ELEMENTclass[0]->NfnIs[Pb][0])*(ELEMENTclass[1]->NvnIs[Pb]);
				NfnIs[Pb][1] = (ELEMENTclass[0]->NvnIs[Pb])*   (ELEMENTclass[1]->NfnIs[Pb][0]);
				NfnIc[Pb][0] = (ELEMENTclass[0]->NfnIc[Pb][0])*(ELEMENTclass[1]->NvnIc[Pb]);
				NfnIc[Pb][1] = (ELEMENTclass[0]->NvnIc[Pb])*   (ELEMENTclass[1]->NfnIc[Pb][0]);

				get_sf_parameters(1,ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->w_vIs[Pb],
				                  1,ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->w_vIs[Pb],NIn,NOut,OP,dE,3,Eclass);
				if (w_vIs[Pb] == NULL)
					w_vIs[Pb] = sf_assemble_d(NIn,NOut,dE,OP);
				get_sf_parameters(1,ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->w_vIc[Pb],
				                  1,ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->w_vIc[Pb],NIn,NOut,OP,dE,3,Eclass);
				if (w_vIc[Pb] == NULL)
					w_vIc[Pb] = sf_assemble_d(NIn,NOut,dE,OP);

				for (vh = 0; vh < Nvref; vh++) {
					if (vh == 0 || P == Pb) {
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->Ihat_vS_vS[P][Pb],
						                   ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[Pb],ELEMENTclass[1]->Ihat_vS_vS[P][Pb],
						                   NIn,NOut,OP,dE,vh,Eclass);
						Ihat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parametersV(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[Pb],ELEMENTclass[0]->L2hat_vS_vS[P][Pb],
						                   ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[Pb],ELEMENTclass[1]->L2hat_vS_vS[P][Pb],
						                   NIn,NOut,OP,dE,vh,Eclass);
						L2hat_vS_vS[P][Pb][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (P == PGlobal && Pb == PGlobal) {
							get_sf_parametersV(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->I_vGs_vGs[1][1],
							                   ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->I_vGs_vGs[1][1],
							                   NIn,NOut,OP,dE,vh,Eclass);
							I_vGs_vGs[1][1][vh] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
					}
				}

				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->ChiS_vP[P][Pb][0],
				                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnP[Pb],ELEMENTclass[1]->ChiS_vP[P][Pb][0],
				                  NIn,NOut,OP,dE,3,Eclass);
				ChiS_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->ChiS_vIc[P][Pb][0],
				                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->ChiS_vIc[P][Pb][0],
				                  NIn,NOut,OP,dE,3,Eclass);
				ChiS_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
				get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnP[Pb],ELEMENTclass[0]->I_vGc_vP[P][Pb][0],
				                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnP[Pb],ELEMENTclass[1]->I_vGc_vP[P][Pb][0],
				                  NIn,NOut,OP,dE,3,Eclass);
				I_vGc_vP[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep

				// Note: Most VOLUME operators need not interpolate between different orders
				if (P == Pb) {
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnP[P],ELEMENTclass[0]->I_vGs_vP[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnP[P],ELEMENTclass[1]->I_vGs_vP[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vP[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->ChiInvS_vS[P][Pb][0],
					                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->ChiInvS_vS[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					ChiInvS_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->ChiS_vIs[P][Pb][0],
					                  ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->ChiS_vIs[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					ChiS_vIs[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnGc[Pb],ELEMENTclass[0]->I_vGs_vGc[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnGc[Pb],ELEMENTclass[1]->I_vGs_vGc[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vGc[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnCs[Pb],ELEMENTclass[1]->I_vGs_vCs[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vCs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGs_vS[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->I_vGs_vS[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vS[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],
					                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->I_vGs_vIs[1][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGs_vIs[1][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnCc[Pb],ELEMENTclass[0]->I_vGc_vCc[P][Pb][0],
					                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnCc[Pb],ELEMENTclass[1]->I_vGc_vCc[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vCc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->I_vGc_vS[P][Pb][0],
					                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->I_vGc_vS[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vS[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					get_sf_parameters(ELEMENTclass[0]->NvnGc[P],ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->I_vGc_vIc[P][Pb][0],
					                  ELEMENTclass[1]->NvnGc[P],ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->I_vGc_vIc[P][Pb][0],
					                  NIn,NOut,OP,dE,3,Eclass);
					I_vGc_vIc[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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

					if (EFE) {
						get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Is_Weak_VV[P][Pb][0],
						                  ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->Is_Weak_VV[P][Pb][0],
						                  NIn,NOut,OP,dE,3,Eclass);
						Is_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],
						                  ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->NvnS[P],ELEMENTclass[1]->Ic_Weak_VV[P][Pb][0],
						                  NIn,NOut,OP,dE,3,Eclass);
						Ic_Weak_VV[P][Pb][0] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					} else {
						printf("Error: Unsupported EFE.\n"), EXIT_MSG;
					}

					for (dim = 0; dim < dE; dim++) {
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGs_vCs[1][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGs_vCs[1][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGs_vCs[1][Pb][0],      OP1 = ELEMENTclass[1]->D_vGs_vCs[1][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnCs[Pb],OP0,
						                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnCs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGs_vCs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						if (dim < 2) OP0 = ELEMENTclass[0]->D_vGs_vIs[1][Pb][0][dim], OP1 = ELEMENTclass[1]->I_vGs_vIs[1][Pb][0];
						else         OP0 = ELEMENTclass[0]->I_vGs_vIs[1][Pb][0],      OP1 = ELEMENTclass[1]->D_vGs_vIs[1][Pb][0][0];
						get_sf_parameters(ELEMENTclass[0]->NvnGs[1],ELEMENTclass[0]->NvnIs[Pb],OP0,
						                  ELEMENTclass[1]->NvnGs[1],ELEMENTclass[1]->NvnIs[Pb],OP1,NIn,NOut,OP,dE,3,Eclass);
						D_vGs_vIs[1][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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
							get_sf_parameters(ELEMENTclass[0]->NvnIs[Pb],ELEMENTclass[0]->NvnS[P],OP0,
							                  ELEMENTclass[1]->NvnIs[Pb],ELEMENTclass[1]->NvnS[P],OP1,NIn,NOut,OP,dE,3,Eclass);
							Ds_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							if (dim < 2) OP0 = ELEMENTclass[0]->Dc_Weak_VV[P][Pb][0][dim], OP1 = ELEMENTclass[1]->Ic_Weak_VV[P][Pb][0];
							else         OP0 = ELEMENTclass[0]->Ic_Weak_VV[P][Pb][0],      OP1 = ELEMENTclass[1]->Dc_Weak_VV[P][Pb][0][0];
							get_sf_parameters(ELEMENTclass[0]->NvnIc[Pb],ELEMENTclass[0]->NvnS[P],OP0,
							                  ELEMENTclass[1]->NvnIc[Pb],ELEMENTclass[1]->NvnS[P],OP1,NIn,NOut,OP,dE,3,Eclass);
							Dc_Weak_VV[P][Pb][0][dim] = sf_assemble_d(NIn,NOut,dE,OP); // keep

							if (Collocated) {
								convert_to_CSR_d(NvnS[P],NvnIs[Pb],Ds_Weak_VV[P][Pb][0][dim],&Ds_Weak_VV_sp[P][Pb][0][dim]); // keep
								convert_to_CSR_d(NvnS[P],NvnIc[Pb],Dc_Weak_VV[P][Pb][0][dim],&Dc_Weak_VV_sp[P][Pb][0][dim]); // keep
							}
						} else {
							;
						}
					}
				}

				for (f = 0; f < Nf; f++) {
				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;

					if (f < 3) { OPF0  = ELEMENTclass[0]->ChiS_fS[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_vS[P][Pb];
					             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],    NOut1 = ELEMENTclass[1]->NvnS[Pb];
					} else {     OPF0  = ELEMENTclass[0]->ChiS_vS[P][Pb], OPF1  = ELEMENTclass[1]->ChiS_fS[P][Pb];
					             NOut0 = ELEMENTclass[0]->NvnS[Pb],       NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
					get_sf_parametersF(ELEMENTclass[0]->NvnS[P],NOut0,OPF0,
					                   ELEMENTclass[1]->NvnS[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					ChiS_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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

					if (Collocated || VFPartUnity[Eclass]) {
						if (f < 3) IndClass = 0;
						else       IndClass = 1;

						convert_to_CSR_d(NfnIs[Pb][IndClass],NvnS[P],ChiS_fIs[P][Pb][Vf],&ChiS_fIs_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NfnIc[Pb][IndClass],NvnS[P],ChiS_fIc[P][Pb][Vf],&ChiS_fIc_sp[P][Pb][Vf]); // keep
					}

					if (P == Pb) {
						// Outside of fh == 0 if condition as operator is used for testing h-adaptation connectivity.
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fS[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vS[1][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGs_vS[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fS[1][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],NOut0,OPF0,
						                   ELEMENTclass[1]->NvnGs[1],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGs_fS[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					}

					if (f < 3) { OPF0 = ELEMENTclass[0]->Is_Weak_FF[P][Pb], OPF1 = ELEMENTclass[1]->Is_Weak_VV[P][Pb];
					             NIn0 = ELEMENTclass[0]->NfnIs[Pb][0],      NIn1 = ELEMENTclass[1]->NvnIs[Pb];
					} else {     OPF0 = ELEMENTclass[0]->Is_Weak_VV[P][Pb], OPF1 = ELEMENTclass[1]->Is_Weak_FF[P][Pb];
					             NIn0 = ELEMENTclass[0]->NvnIs[Pb],         NIn1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
					get_sf_parametersF(NIn0,ELEMENTclass[0]->NvnS[P],OPF0,
					                   NIn1,ELEMENTclass[1]->NvnS[P],OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					Is_Weak_FF[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
					if (f < 3) { OPF0 = ELEMENTclass[0]->Ic_Weak_FF[P][Pb], OPF1 = ELEMENTclass[1]->Ic_Weak_VV[P][Pb];
					             NIn0 = ELEMENTclass[0]->NfnIc[Pb][0],      NIn1 = ELEMENTclass[1]->NvnIc[Pb];
					} else {     OPF0 = ELEMENTclass[0]->Ic_Weak_VV[P][Pb], OPF1 = ELEMENTclass[1]->Ic_Weak_FF[P][Pb];
					             NIn0 = ELEMENTclass[0]->NvnIc[Pb],         NIn1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
					get_sf_parametersF(NIn0,ELEMENTclass[0]->NvnS[P],OPF0,
					                   NIn1,ELEMENTclass[1]->NvnS[P],OPF1,NIn,NOut,OP,dE,Vf,Eclass);
					Ic_Weak_FF[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep

					if (Collocated || VFPartUnity[Eclass]) {
						if (f < 3) IndClass = 0;
						else       IndClass = 1;

						convert_to_CSR_d(NvnS[P],NfnIs[Pb][IndClass],Is_Weak_FF[P][Pb][Vf],&Is_Weak_FF_sp[P][Pb][Vf]); // keep
						convert_to_CSR_d(NvnS[P],NfnIc[Pb][IndClass],Ic_Weak_FF[P][Pb][Vf],&Ic_Weak_FF_sp[P][Pb][Vf]); // keep
					}

					if (fh == 0) {
						if (P == Pb) {
							if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fIs[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vIs[1][Pb];
							             NOut0 = ELEMENTclass[0]->NfnIs[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIs[Pb];
							} else {     OPF0  = ELEMENTclass[0]->I_vGs_vIs[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fIs[1][Pb];
							             NOut0 = ELEMENTclass[0]->NvnIs[Pb],        NOut1 = ELEMENTclass[1]->NfnIs[Pb][0]; }
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],NOut0,OPF0,
							                   ELEMENTclass[1]->NvnGs[1],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIs[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
							if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGs_fIc[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_vIc[1][Pb];
							             NOut0 = ELEMENTclass[0]->NfnIc[Pb][0],     NOut1 = ELEMENTclass[1]->NvnIc[Pb];
							} else {     OPF0  = ELEMENTclass[0]->I_vGs_vIc[1][Pb], OPF1  = ELEMENTclass[1]->I_vGs_fIc[1][Pb];
							             NOut0 = ELEMENTclass[0]->NvnIc[Pb],        NOut1 = ELEMENTclass[1]->NfnIc[Pb][0]; }
							get_sf_parametersF(ELEMENTclass[0]->NvnGs[1],NOut0,OPF0,
							                   ELEMENTclass[1]->NvnGs[1],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
							I_vGs_fIc[1][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
						}
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vGc_fS[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_vS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vGc_vS[P][Pb], OPF1  = ELEMENTclass[1]->I_vGc_fS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnGc[P],NOut0,OPF0,
						                   ELEMENTclass[1]->NvnGc[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vGc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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

						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCs_fS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_vS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCs_vS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCs_fS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCs[P],NOut0,OPF0,
						                   ELEMENTclass[1]->NvnCs[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCs_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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
						if (f < 3) { OPF0  = ELEMENTclass[0]->I_vCc_fS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_vS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NfnS[Pb][0],     NOut1 = ELEMENTclass[1]->NvnS[Pb];
						} else {     OPF0  = ELEMENTclass[0]->I_vCc_vS[P][Pb], OPF1  = ELEMENTclass[1]->I_vCc_fS[P][Pb];
						             NOut0 = ELEMENTclass[0]->NvnS[Pb],        NOut1 = ELEMENTclass[1]->NfnS[Pb][0]; }
						get_sf_parametersF(ELEMENTclass[0]->NvnCc[P],NOut0,OPF0,
						                   ELEMENTclass[1]->NvnCc[P],NOut1,OPF1,NIn,NOut,OP,dE,Vf,Eclass);
						I_vCc_fS[P][Pb][Vf] = sf_assemble_d(NIn,NOut,dE,OP); // keep
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
		printf("Error: Unsupported Eclass.\n"), EXIT_MSG;
	}
	free(ones_Nf);
}

static void setup_L2_projection_preoperators(const unsigned int EType)
{
	// Returned operators
	unsigned int *NvnIc, *NvnS;
	double       **w_vIc, ****ChiS_vIc;

	// Initialize DB Parameters
	unsigned int PGs          = DB.PGs,
	             NP           = DB.NP;

	char         *BasisType     = DB.BasisType,
	             ***NodeTypeIvc = DB.NodeTypeIvc,
	             ***NodeTypeS   = DB.NodeTypeS;

	// Standard datatypes
	unsigned int i, iMax, dE, P, vh, PSMin, PSMax, Pb, PbMin, PbMax, PIvc[NEC],
	             Nve, Nbf, Eclass, Nvref, NEhref, Indh, *Nvve, *EType_h, dummy_ui, *dummyPtr_ui[2];
	double       *E_rst_vC, *rst_vC, **VeV, **rst_vIc, *rst_vS,
	             *IGs, *IS, *TS, *ChiRefS_vS, *ChiRefInvGs_vGs, *ChiRefInvS_vS,
	             *ChiRefGs_vGs, *ChiRefGs_vIc, *ChiRefS_vIc, *ChiS_vS, *dummyPtr_d;

	struct BCoords {
		double **Ic;
	} **BCoords_V;
	struct S_ELEMENT *ELEMENT, *ELEMENT_h;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;
	grad_basis_tdef grad_basis;

	// silence
	E_rst_vC = ChiS_vS = NULL;

	PIvc[C_SI]  = PIvcMaxTET;
	PIvc[C_PYR] = PIvcMaxPYR;
	if (!(EType == TET || EType == PYR))
		printf("Error: Unsupported EType.\n"), EXIT_MSG;

	ELEMENT   = get_ELEMENT_type(EType);
	ELEMENT_h = ELEMENT;

	Eclass = get_Eclass(EType);
	setup_ELEMENT_VeV(EType);

	dE      = ELEMENT->d;
	Nve     = ELEMENT->Nve;
	VeV     = ELEMENT->VeV;
	Nvve    = ELEMENT->Nvve;
	Nvref   = ELEMENT->Nvref;
	NEhref  = ELEMENT->NEhref;
	EType_h = ELEMENT->type_h;

	select_functions(&basis,&grad_basis,&cubature,EType);


	// Stored operators
	NvnIc    = ELEMENT->NvnIc;
	NvnS     = ELEMENT->NvnS;
	w_vIc    = ELEMENT->w_vIc;
	ChiS_vIc = ELEMENT->ChiS_vIc;

	// Allocate memory for arrays with multiple levels of dereferencing
	rst_vIc = malloc(NVREFMAX * sizeof *rst_vIc); // free

	BCoords_V = malloc(NEhref * sizeof *BCoords_V); // free
	for (i = 0, iMax = NEhref; i < iMax; i++) {
		BCoords_V[i]     = malloc(     sizeof *BCoords_V[i]);       // free
		BCoords_V[i]->Ic = calloc(NP , sizeof *(BCoords_V[i]->Ic)); // free
	}

	rst_vC = malloc(Nve*dE * sizeof *rst_vC); // free

	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {

		cubature(&rst_vS,&dummyPtr_d,&dummyPtr_ui[0],&NvnS[P],&dummy_ui,0,P,dE,NodeTypeS[P][Eclass]); free(dummyPtr_ui[0]); // free

		IS         = identity_d(NvnS[P]);             // free
		ChiRefS_vS = basis(P,rst_vS,NvnS[P],&Nbf,dE); // free
		free(rst_vS);
		if (strstr(BasisType,"Modal")) {
			ChiS_vS = ChiRefS_vS;
		} else if (strstr(BasisType,"Nodal")) {
			ChiS_vS = IS;
		}

		ChiRefInvS_vS = inverse_d(NvnS[P],NvnS[P],ChiRefS_vS,IS);                                     // free
		TS            = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[P],NvnS[P],NvnS[P],1.0,ChiRefInvS_vS,ChiS_vS); // free

		free(IS);
		free(ChiRefS_vS);
		free(ChiRefInvS_vS);

		get_Pb_range(P,&PbMin,&PbMax);
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			if (w_vIc[Pb])
				free(w_vIc[Pb]);
			cubature(&rst_vIc[0],&w_vIc[Pb],&dummyPtr_ui[0],&NvnIc[Pb],&dummy_ui,1,PIvc[Eclass],dE,NodeTypeIvc[Pb][Eclass]); free(dummyPtr_ui[0]); // free
			free(rst_vIc[0]);

			for (i = iMax = NEhref; i--; ) {
				ELEMENT = get_ELEMENT_type(EType_h[i]);

				Eclass = get_Eclass(EType_h[i]);
				select_functions(&basis,&grad_basis,&cubature,EType_h[i]);

				Nve = ELEMENT->Nve;

				cubature(&rst_vIc[0],&dummyPtr_d,&dummyPtr_ui[0],&NvnIc[Pb],&dummy_ui,0,PIvc[Eclass],dE,NodeTypeIvc[Pb][Eclass]); free(dummyPtr_ui[0]); // free

				E_rst_vC        = get_rst_vC(ELEMENT);                     // free
				IGs             = identity_d(Nve);                         // free
				ChiRefGs_vGs    = basis(PGs,E_rst_vC,Nve,&Nbf,dE);         // free
				ChiRefInvGs_vGs = inverse_d(Nve,Nve,ChiRefGs_vGs,IGs);     // free
				ChiRefGs_vIc    = basis(PGs,rst_vIc[0],NvnIc[Pb],&Nbf,dE); // free

				BCoords_V[i]->Ic[Pb] = mm_Alloc_d(CBCM,CBT,CBT,NvnIc[Pb],Nve,Nve,1.0,ChiRefGs_vIc,ChiRefInvGs_vGs); // free

				free(IGs);
				free(ChiRefGs_vGs);
				free(ChiRefInvGs_vGs);
				free(ChiRefGs_vIc);

				if (i) {
					free(E_rst_vC);
					free(rst_vIc[0]);
				}
			}

			for (vh = 0; vh < Nvref; vh++) {
				mm_CTN_d(Nvve[vh],dE,Nve,VeV[vh],E_rst_vC,rst_vC);
				if (vh) {
					Indh = get_IndEhref(EType,vh);
					if (EType == PYR)
						ELEMENT_h = get_ELEMENT_type(EType_h[Indh]);
					else
						ELEMENT_h = ELEMENT;

					rst_vIc[vh] = mm_Alloc_d(CBCM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],dE,Nvve[vh],1.0,BCoords_V[Indh]->Ic[Pb],rst_vC); // free
				}

				ChiRefS_vIc         = basis(P,rst_vIc[vh],ELEMENT_h->NvnIc[Pb],&Nbf,dE);                                  // free
				ChiS_vIc[P][Pb][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,ELEMENT_h->NvnIc[Pb],NvnS[P],NvnS[P],1.0,ChiRefS_vIc,TS); // keep

				free(ChiRefS_vIc);
			}

			free(E_rst_vC);
			for (vh = 0; vh < Nvref; vh++)
				free(rst_vIc[vh]);

			for (i = 0; i < NEhref; i++)
				free(BCoords_V[i]->Ic[Pb]);
		}
		free(TS);
	}
	free(rst_vC);
	free(rst_vIc);
	for (i = 0; i < NEhref; i++) {
		free(BCoords_V[i]->Ic);
		free(BCoords_V[i]);
	}
	free(BCoords_V);
}

static void setup_L2_projection_operators(const unsigned int EType)
{
	/*
	 *	Comments:
	 *		This function is separated from the setup_ELEMENT_operators function as the order P operators defined here
	 *		rely on order Pb operators for their assembly.
	 *		L2 projection operators are only needed for coarsening operations. During refinement, these operators reduce
	 *		to standard interpolations.
	 *		For PYR ELEMENTs, the order used for the cubature nodes has a signifiant impact on the accuracy of the L2
	 *		projections (several orders of magnitude), which remain quite poor even when using the highest order
	 *		cubature nodes supported (TET: P10, PYR: P6). As the cubature order is not seen in the L2 projection
	 *		operator, it is advantageous to compute it using the highest possible cubature order. However, this order is
	 *		not required for other operators in the code and is thus redefined appropriately in setup_operators.
	 */

	// Returned Operators
	double ****L2hat_vS_vS;

	// Initialize DB Parameters
	unsigned int Adapt = DB.Adapt;

	// Standard datatypes
	unsigned int P, Pb, vh, PSMin, PSMax, PbMin, PbMax,
	             Indh, Nvref, *NvnS, *NvnIc, *EType_h;
	double       L2_scale, *IS, *ChiTW, *M, *MInv, *S, *diag_w_vIc,
	             ****ChiS_vIc, **w_vIc;

	struct S_ELEMENT *ELEMENT, *ELEMENT_h;

	if (Adapt == ADAPT_0)
		return;

	ELEMENT   = get_ELEMENT_type(EType);

	if (EType == TET || EType == PYR)
		setup_L2_projection_preoperators(EType);

	ELEMENT_h = ELEMENT;

	// Stored operators
	L2hat_vS_vS = ELEMENT->L2hat_vS_vS;

	EType_h = ELEMENT->type_h;

	Nvref = ELEMENT->Nvref;
	NvnS  = ELEMENT->NvnS;
	NvnIc = ELEMENT->NvnIc;

	ChiS_vIc = ELEMENT->ChiS_vIc;
	w_vIc    = ELEMENT->w_vIc;

	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {
		get_Pb_range(P,&PbMin,&PbMax);
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			diag_w_vIc = diag_d(w_vIc[Pb],NvnIc[Pb]); // free

			IS    = identity_d(NvnS[Pb]);                                                                       // free
			ChiTW = mm_Alloc_d(CBRM,CBT, CBNT,NvnS[Pb],NvnIc[Pb],NvnIc[Pb],1.0,ChiS_vIc[Pb][Pb][0],diag_w_vIc); // free
			M     = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[Pb],NvnS[Pb], NvnIc[Pb],1.0,ChiTW,ChiS_vIc[Pb][Pb][0]);      // free
			MInv  = inverse_d(NvnS[Pb],NvnS[Pb],M,IS);                                                          // free

			free(diag_w_vIc);
			free(IS);
			free(ChiTW);
			free(M);
			for (vh = 0; vh < Nvref; vh++) {
				if (vh == 0 || P == Pb) {
					Indh = get_IndEhref(EType,vh);
					if (EType == PYR)
						ELEMENT_h = get_ELEMENT_type(EType_h[Indh]);
					else
						ELEMENT_h = ELEMENT;

					diag_w_vIc = diag_d(ELEMENT_h->w_vIc[Pb],ELEMENT_h->NvnIc[Pb]); // free

					ChiTW = mm_Alloc_d(CBRM,CBT, CBNT,NvnS[Pb],ELEMENT_h->NvnIc[Pb],ELEMENT_h->NvnIc[Pb],1.0,ChiS_vIc[Pb][Pb][vh],diag_w_vIc);     // free
					S     = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[Pb],ELEMENT_h->NvnS[P],  ELEMENT_h->NvnIc[Pb],1.0,ChiTW,ELEMENT_h->ChiS_vIc[P][Pb][0]); // free

					L2_scale = get_L2_scaling(EType,vh);

					// Returned Operators
					L2hat_vS_vS[P][Pb][vh] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS[Pb],ELEMENT_h->NvnS[P],NvnS[Pb],L2_scale,MInv,S); // keep

					free(diag_w_vIc);
					free(ChiTW);
					free(S);
				}
			}
			free(MInv);
		}
	}


/*
	// Returned Operators
	double ****L2hat_vS_vS,
	       ****GfS_fIs, ****GfS_fIc; // Change G to L2 if these operators are kept (ToBeDeleted)

	// Initialize DB Parameters
	unsigned int Adapt = DB.Adapt;

	// Standard datatypes
	unsigned int P, Pb, vref, f, fh, PSMin, PSMax, PbMin, PbMax, Nf, fhMax,
	             Vf, PM, *Nfref, Eclass, IndFType, Nvref,
	             NvnS, NvnSN, NvnI, NfnM, NfnF, NfnIs, NfnIc,
	             *ones_Nf;
	double       L2_scale, *ChiSN_vI, *wSN_vI, *diag_wSN_vI, *ISN, *ChiTW_SN, *MSN, *MSNInv, *ChiS_vI, *SSNS,
	             *one_d, *ChiS_vF, *ChiM_vI, *ChiF_vIs, *ChiF_vIc,
	             *ChiInvF_vF, *ChiInvM_vM, *wM_vI, *diag_wM_vI, *IF, *IM, *ChiTW_M, *ChiTW_F, *ChiTW_SNr,
	             *SMF, *SFM, *MM, *MF, *MMInv, *MFInv, *MInvSMF, *MInvSFM,
	             *IhatS_fS, *L2hatvS_fS, *L2hatfS_fI;

	struct S_ELEMENT *ELEMENT, *ELEMENT_F;

	if (Adapt == ADAPT_0)
		return;

	// silence
	IF = IM = NULL;
	diag_wM_vI = ChiTW_M = NULL;
	MM = MMInv = NULL;

	one_d = malloc(1 * sizeof *one_d); // free
	one_d[0] = 1.0;

	ELEMENT = get_ELEMENT_type(EType);
	Eclass = get_Eclass(EType);

	Nf    = ELEMENT->Nf;
	Nfref = ELEMENT->Nfref;
	Nvref = ELEMENT->Nvref;

	// Stored operators
	L2hat_vS_vS = ELEMENT->L2hat_vS_vS;
	GfS_fIs     = ELEMENT->GfS_fIs;
	GfS_fIc     = ELEMENT->GfS_fIc;


	ones_Nf = malloc(Nf * sizeof *ones_Nf); // free
	for (f = 0; f < Nf; f++)
		ones_Nf[f] = 1;

	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {
		NvnS = ELEMENT->NvnS[P];

		get_Pb_range(P,&PbMin,&PbMax);
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			PM = max(P,Pb);

			NvnSN = ELEMENT->NvnS[Pb]; // N: New
			NvnI  = ELEMENT->NvnIc[Pb];

			ChiSN_vI = ELEMENT->ChiS_vIc[Pb][Pb][0];

			wSN_vI      = ELEMENT->w_vIc[Pb];
			diag_wSN_vI = diag_d(wSN_vI,NvnI); // free

			ISN      = identity_d(NvnSN); // free
			ChiTW_SN = mm_Alloc_d(CBRM,CBT, CBNT,NvnSN,NvnI, NvnI,1.0,ChiSN_vI,diag_wSN_vI); // free
			MSN      = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSN,NvnSN,NvnI,1.0,ChiTW_SN,ChiSN_vI);    // free
			MSNInv   = inverse_d(NvnSN,NvnSN,MSN,ISN);                                       // free

			ChiS_vI  = ELEMENT->ChiS_vIc[P][Pb][0];
			for (vref = 0; vref < Nvref; vref++) {
				ChiSN_vI = ELEMENT->ChiS_vIc[Pb][Pb][vref];

				ChiTW_SNr = mm_Alloc_d(CBRM,CBT, CBNT,NvnSN,NvnI,NvnI,1.0,ChiSN_vI, diag_wSN_vI); // free
				SSNS      = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSN,NvnS,NvnI,1.0,ChiTW_SNr,ChiS_vI);     // free

				L2_scale = get_L2_scaling(EType,vref);

				// Returned Operators
				L2hat_vS_vS[P][Pb][vref] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSN,NvnS,NvnSN,L2_scale,MSNInv,SSNS); // keep

				free(ChiTW_SNr);
				free(SSNS);
			}

			free(diag_wSN_vI);
			free(ISN);
			free(ChiTW_SN);
			free(MSN);
			free(MSNInv);
		}

		for (Pb = P; Pb <= PbMax; Pb++) {
			PM = max(P,Pb);
			for (f = 0; f < Nf; f++) {
				IndFType = get_IndFType(Eclass,f);
				ELEMENT_F = get_ELEMENT_FACET(EType,IndFType);

				if (EType == LINE) {
					NfnM = 1; // (M)ortar
					NfnF = 1; // (F)acet
					NfnIs = 1;
					NfnIc = 1;

//					ChiM_vS   = one_d;
				} else {
					NfnM = ELEMENT_F->NvnS[PM];
					NfnF = ELEMENT_F->NvnS[Pb];
					NfnIs = ELEMENT_F->NvnIs[PM];
					NfnIc = ELEMENT_F->NvnIc[PM];

//					ChiM_vS  = ELEMENT_F->ChiS_vS[PM][PM][0];
					ChiM_vI  = ELEMENT_F->ChiS_vIc[PM][PM][0];
					wM_vI    = ELEMENT_F->w_vIc[PM];

					ChiInvF_vF = ELEMENT_F->ChiInvS_vS[Pb][Pb][0];
					ChiInvM_vM = ELEMENT_F->ChiInvS_vS[PM][PM][0];

					IF   = identity_d(NfnF); // free
					IM   = identity_d(NfnM); // free

					diag_wM_vI = diag_d(wM_vI,NfnIc); // free
					ChiTW_M    = mm_Alloc_d(CBRM,CBT,CBNT,NfnM,NfnIc,NfnIc,1.0,ChiM_vI,diag_wM_vI); // free
					MM         = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NfnM,NfnIc,1.0,ChiTW_M,ChiM_vI);    // free
					MMInv      = inverse_d(NfnM,NfnM,MM,IM);                                        // free
				}

				for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
					Vf = f*NFREFMAX+fh;
//printf("P, Pb, f, fh, Vf: %d %d %d %d %d\n",P,Pb,f,fh,Vf);

					ChiS_vF = ELEMENT->ChiS_fS[P][Pb][f*NFREFMAX];

					if (EType == LINE) {
						ChiF_vIs  = one_d;
						ChiF_vIc  = one_d;
						L2hatvS_fS = mm_Alloc_d(CBRM,CBNT,CBNT,1,NvnS,1,1.0,one_d,ChiS_vF); // free
						L2hatfS_fI = mm_Alloc_d(CBRM,CBNT,CBNT,1,1,1,1.0,one_d,one_d);      // free
					} else {
						ChiF_vIs = ELEMENT_F->ChiS_vIs[Pb][PM][fh];
						ChiF_vIc = ELEMENT_F->ChiS_vIc[Pb][PM][fh];

						ChiTW_F  = mm_Alloc_d(CBRM,CBT,CBNT,NfnF,NfnIc,NfnIc,1.0,ChiF_vIc,diag_wM_vI); // free

						SMF      = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NfnF,NfnIc,1.0,ChiTW_M,ChiF_vIc);  // free
						SFM      = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnM,NfnIc,1.0,ChiTW_M,ChiF_vIc);  // free
						MF       = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnF,NfnIc,1.0,ChiTW_F,ChiF_vIc);  // free
						MFInv    = inverse_d(NfnF,NfnF,MF,IF);                                       // free
						MInvSMF  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NfnF,NfnM,1.0,MMInv,SMF);          // free
						MInvSFM  = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnM,NfnF,1.0,MFInv,SFM);          // free
						IhatS_fS = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NvnS,NfnF,1.0,ChiInvF_vF,ChiS_vF); // free

						L2hatvS_fS = mm_Alloc_d(CBRM,CBNT,CBNT,NfnM,NvnS,NfnF,1.0,MInvSMF,IhatS_fS);   // free
						L2hatfS_fI = mm_Alloc_d(CBRM,CBNT,CBNT,NfnF,NfnM,NfnM,1.0,MInvSFM,ChiInvM_vM); // free

						free(ChiTW_F);
						free(SMF), free(SFM), free(MF);
						free(MFInv), free(MInvSMF), free(MInvSFM);
						free(IhatS_fS);
					}

					// Returned Operators
					GfS_fIs[P][Pb][Vf]   = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIs,NfnM,NfnF,1.0,ChiF_vIs,L2hatfS_fI); // keep
					GfS_fIc[P][Pb][Vf]   = mm_Alloc_d(CBRM,CBNT,CBNT,NfnIc,NfnM,NfnF,1.0,ChiF_vIc,L2hatfS_fI); // keep

					free(L2hatvS_fS);
					free(L2hatfS_fI);
				}
				free(IF), free(IM);
				free(diag_wM_vI), free(ChiTW_M);
				free(MM);
				free(MMInv);
			}
		}
	}
	free(one_d);
	free(ones_Nf);
*/
}

static void setup_ELEMENT_FACET_ordering(const unsigned int FType)
{
	// Returned operators
	unsigned int ***nOrd_fS, ***nOrd_fIs, ***nOrd_fIc;

	// Initialize DB Parameters
	const unsigned int d    = DB.d,
	                   PMax = DB.PMax;

	// Standard datatypes
	unsigned int P, NOrd, IndOrd;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(FType);

	if (FType != POINT) {
		// Initialize DB Parameters
		unsigned int **PIfs         = DB.PIfs,
		             **PIfc         = DB.PIfc;
		char         ***NodeTypeS   = DB.NodeTypeS,
		             ***NodeTypeIfs = DB.NodeTypeIfs,
		             ***NodeTypeIfc = DB.NodeTypeIfc;

		// Standard datatypes
		unsigned int dE, NfnS, NfnIs, NfnIc, NsS, NsIs, NsIc, Eclass,
		             *symms_fS, *symms_fIs, *symms_fIc;
		double       *rst_fS, *rst_fIs, *rst_fIc, *w;

		// Function pointers
		cubature_tdef   cubature;
		basis_tdef      basis;
		grad_basis_tdef grad_basis;

		select_functions(&basis,&grad_basis,&cubature,FType);

		// silence
		NfnS = NfnIs = NfnIc = NsS = NsIs = NsIc = 0;

		Eclass  = get_Eclass(FType);

		dE = d-1;

		if      (FType == LINE)  NOrd = 2;
		else if (FType == QUAD)  NOrd = 8;
		else if (FType == TRI)   NOrd = 6;
		else    printf("Error: Unsupported FType.\n"), EXIT_MSG;

		nOrd_fS  = ELEMENT->nOrd_fS;
		nOrd_fIs = ELEMENT->nOrd_fIs;
		nOrd_fIc = ELEMENT->nOrd_fIc;

		for (P = 0; P <= PMax; P++) {
			cubature(&rst_fS, &w,&symms_fS, &NfnS, &NsS, 0,P,              dE,NodeTypeS[P][Eclass]);   // free
			cubature(&rst_fIs,&w,&symms_fIs,&NfnIs,&NsIs,0,PIfs[P][Eclass],dE,NodeTypeIfs[P][Eclass]); // free
			cubature(&rst_fIc,&w,&symms_fIc,&NfnIc,&NsIc,0,PIfc[P][Eclass],dE,NodeTypeIfc[P][Eclass]); // free

			for (IndOrd = 0; IndOrd < NOrd; IndOrd++) {
				nOrd_fS[P][IndOrd]  = malloc(NfnS  * sizeof ***nOrd_fS);  //  keep
				nOrd_fIs[P][IndOrd] = malloc(NfnIs * sizeof ***nOrd_fIs); //  keep
				nOrd_fIc[P][IndOrd] = malloc(NfnIc * sizeof ***nOrd_fIc); //  keep

				get_facet_ordering(d,IndOrd,FType,NfnS, NsS, symms_fS, rst_fS, nOrd_fS[P][IndOrd]);
				get_facet_ordering(d,IndOrd,FType,NfnIs,NsIs,symms_fIs,rst_fIs,nOrd_fIs[P][IndOrd]);
				get_facet_ordering(d,IndOrd,FType,NfnIc,NsIc,symms_fIc,rst_fIc,nOrd_fIc[P][IndOrd]);
			}
			free(rst_fS);
			free(rst_fIs);
			free(rst_fIc);
			free(symms_fS);
			free(symms_fIs);
			free(symms_fIc);
		}
	} else {
		NOrd = 1;

		nOrd_fS  = ELEMENT->nOrd_fS;
		nOrd_fIs = ELEMENT->nOrd_fIs;
		nOrd_fIc = ELEMENT->nOrd_fIc;

		for (P = 0; P <= PMax; P++) {
		for (IndOrd = 0; IndOrd < NOrd; IndOrd++) {
			nOrd_fS[P][IndOrd]  = malloc(1 * sizeof ***nOrd_fS);  //  keep
			nOrd_fIs[P][IndOrd] = malloc(1 * sizeof ***nOrd_fIs); //  keep
			nOrd_fIc[P][IndOrd] = malloc(1 * sizeof ***nOrd_fIc); //  keep

			get_facet_ordering(d,0,0,0,0,NULL,NULL,nOrd_fS[P][IndOrd]);
			get_facet_ordering(d,0,0,0,0,NULL,NULL,nOrd_fIs[P][IndOrd]);
			get_facet_ordering(d,0,0,0,0,NULL,NULL,nOrd_fIc[P][IndOrd]);
		}}
	}
}

void setup_operators(void)
{
	/*
	 *	Comments:
	 *		TETs are set up if PYRs are present as PYRs h-refine into TETs.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int EType;

	// POINT
	EType = POINT;
	if (d == 1)
		setup_ELEMENT_FACET_ordering(EType);

	// LINE (Includes TP Class)
	if (!DB.MPIrank && !DB.Testing)
		printf("    LINE\n");
	EType = LINE;
	setup_ELEMENT_VeF(EType);
	setup_ELEMENT_plotting(EType);
	setup_ELEMENT_normals(EType);
	setup_ELEMENT_operators(EType);
	setup_L2_projection_operators(EType);
	if (d == 2)
		setup_ELEMENT_FACET_ordering(EType);

	// QUAD
	EType = QUAD;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    QUAD\n");
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
		if (!DB.MPIrank && !DB.Testing)
			printf("    HEX\n");
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_TP_operators(EType);
	}

	// TRI
	EType = TRI;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    TRI\n");
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);
		setup_L2_projection_operators(EType);
		if (d == 3)
			setup_ELEMENT_FACET_ordering(EType);
	}

	// TET
	EType = TET;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    TET\n");

		if (is_ELEMENT_present(PYR)) {
			setup_L2_projection_operators(EType);
			setup_L2_projection_operators(PYR);

			memory_destructor_L2_projection(EType);
			memory_destructor_L2_projection(PYR);
		} else {
			setup_L2_projection_operators(EType);
			memory_destructor_L2_projection(EType);
		}
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);
	}

	// PYR
	EType = PYR;
	if (is_ELEMENT_present(EType)) {
		if (!is_ELEMENT_present(TET))
			printf("Error: TETs must be set up if PYRs are used.\n"), EXIT_MSG;

		if (!DB.MPIrank && !DB.Testing)
			printf("    PYR\n");
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_ELEMENT_operators(EType);
	}

	// WEDGE
	EType = WEDGE;
	if (is_ELEMENT_present(EType)) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    WEDGE\n");
		setup_ELEMENT_VeF(EType);
		setup_ELEMENT_plotting(EType);
		setup_ELEMENT_normals(EType);
		setup_TP_operators(EType);
	}

	// Free unused operators (including unused lower dimensional operators (if applicable) and sum factorized operators) (ToBeDeleted)
}
