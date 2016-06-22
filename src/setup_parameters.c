// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up parameters based on inputs obtained in initialization.c.
 *
 *	Comments:
 *		FACET integration nodes must be consistent between FACETs of different element types.
 *
 *		Guidelines: (ToBeModified)
 *			PF           >= P
 *			PFr          >= P (PFr = PF+PC is the order needed to represent Fr exactly)
 *			PI(v/f)(s/c) >= 2*P
 *
 *			TP Elements:
 *				GL  : Best cubature nodes
 *				GLL : Best interpolating nodes
 *
 *			SI Elements:
 *				AO        : Best interpolating nodes
 *
 *				2D (TRI):
 *					WS : Best collocated nodes? => Compare with Taylor(2007) (ToBeModified)
 *					WV : Best cubature nodes?   => Compare with Xiao(2010) (ToBeModified)
 *
 *				3D (TET):
 *					SH : Best collocated nodes? (ToBeModified)
 *					WV : Best cubature nodes?   (ToBeModified)
 *
 *				Note: For simplicity of treating all elements, WS (2D) and SH (3D) are both stored as WSH (2D-3D).
 *
 *			WEDGE Elements:
 *				Combination of TP and SI Element nodes
 *
 *			PYR Elements:
 *				ToBeModified
 *
 *	Notation:
 *		NP       : (N)umber of (P)olynomial orders available
 *		NEC      : (N)umber of (E)lement (C)lasses
 *		PR       : Order of solution to read in from the (R)estart file.
 *		PP       : Order used for (P)lotting.
 *
 *		PG()     : Order used for representation of element (G)eometry.
 *		PC()[]   : Order used for representation of (C)ofactor matrix terms.
 *		PJ()[]   : Order used for representation of (J)acobian determinant terms.
 *		PFr()[]  : Order used for representation of (F)lux in (r)eference space.
 *		           () : (s)traight, (c)urved
 *		           [] : TP [0], SI [1], PYR [2]
 *		PF       : Order used for representation of the (F)lux.
 *
 *		PI(1)(2) : Order used for integration (cubature order).
 *		           (1) : (v)olume, (f)acet
 *		           (2) : (s)traight, (c)urved
 *
 *		SF_BE[0][1][2] : Flag for when (S)um (F)actorization (B)reaks (E)ven.
 *		                 [0] - Polynomial order
 *		                 [1] - 0: VOLUME, 1: FACET
 *		                 [2] - 0: TP (QUAD/HEX), 1: WEDGE
 *		                 Note: It may be beneficial to also add the capability for the fast interpolation/differentiation
 *		                       using the Fourier transformed operators, certainly if running the code in a very high-order
 *		                       setting. As the test space in DPG requires higher order than the trial space, this may be
 *		                       even more relevant if using DPG. Note that the same "trick" can be played for TP elements
 *		                       using a 2D Fourier Transform.
 *		                       See Hesthaven(2000)-Stable_Spectral_Methods_on_Tetrahedral_Elements for details. Also see
 *		                       Fladrich-Stiller(2008)-Improved_Performance_for_Nodal_Spectral_Element_Operators for
 *		                       computational considerations.
 *
 *		AC              : Specifies whether (a)ll elements are (c)urved or not.
 *		ExactGeom       : Move boundary nodes to exact geometry if enabled.
 *		Parametrization : Type of parametrization used in curved elements.
 *
 *		NodeType()[] : Node type used for each type of node () and each type of element [].
 *		               () : (S)olution, (F)lux, (F)lux in (r)eference space, (I)ntegration
 *		                    (f)acet/(v)olume (s)traight/(c)urved
 *		               [] : TP [0], SI [1], PYR [2]
 *		InviscidFluxType : Type of inviscid numerical flux used.
 *		                   Options: LF, ROE
 *		ExplicitSolverType : Type of explicit timestepping scheme used.
 *		                     Options: RK3_(S)trong(S)tability(Preserving), RK4_(L)ow(S)torage
 *
 *	References:
 *
 */

void setup_parameters()
{
	// Initialize DB Parameters
	unsigned int d          = DB.d,
	             PMax       = DB.PMax,
	             ML         = DB.ML,
	             EFE        = DB.EFE,
	             Collocated = DB.Collocated;

	char         *Parametrization,
	             **NodeTypeG,
	             ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
	             ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;
	unsigned int i, u1, u2,
	             P, NP, NEC, IntOrderfs, IntOrderfc, IntOrdervs, IntOrdervc,
	             ***SF_BE,
	             PGs, *PGc, **PCs, **PCc, **PJs, **PJc,
	             *PF, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;

	if (DB.PGlobal > PMax)
		printf("Error: P must be less than or equal PMax.\n"), exit(1);
	if (PMax == 0)
		printf("Error: Please choose PMax > 0.\n"), exit(1);

	u1 = 1;
	u2 = 2;

	NP  = PMax+1;
	NEC = 3;

	SF_BE = malloc(NP * sizeof *SF_BE); // keep
	PGc   = malloc(NP * sizeof *PGc);   // keep
	PCs   = malloc(NP * sizeof *PCs);   // keep
	PCc   = malloc(NP * sizeof *PCc);   // keep
	PJs   = malloc(NP * sizeof *PJs);   // keep
	PJc   = malloc(NP * sizeof *PJc);   // keep
	PF    = malloc(NP * sizeof *PF);    // keep
	PFrs  = malloc(NP * sizeof *PFrs);  // keep
	PFrc  = malloc(NP * sizeof *PFrc);  // keep
	PIfs  = malloc(NP * sizeof *PIfs);  // keep
	PIfc  = malloc(NP * sizeof *PIfc);  // keep
	PIvs  = malloc(NP * sizeof *PIvs);  // keep
	PIvc  = malloc(NP * sizeof *PIvc);  // keep

	Parametrization = malloc(STRLEN_MAX * sizeof *NodeTypeS); // keep
	NodeTypeG       = malloc(NEC * sizeof *NodeTypeG);        // keep
	NodeTypeS       = malloc(NP  * sizeof *NodeTypeS);        // keep
	NodeTypeF       = malloc(NP  * sizeof *NodeTypeF);        // keep
	NodeTypeFrs     = malloc(NP  * sizeof *NodeTypeFrs);      // keep
	NodeTypeFrc     = malloc(NP  * sizeof *NodeTypeFrc);      // keep
	NodeTypeIfs     = malloc(NP  * sizeof *NodeTypeIfs);      // keep
	NodeTypeIfc     = malloc(NP  * sizeof *NodeTypeIfc);      // keep
	NodeTypeIvs     = malloc(NP  * sizeof *NodeTypeIvs);      // keep
	NodeTypeIvc     = malloc(NP  * sizeof *NodeTypeIvc);      // keep


	// Restart and Plotting
	if (DB.Restart == -1) {
		DB.PR = 0;
	} else if (DB.Restart ==  0) {
		if (DB.PGlobal == 0)
			printf("Invalid entry for Restart.\n"), exit(1);
		else
			DB.PR = DB.PGlobal-1;
	} else {
		DB.PR = DB.PGlobal;
	}

	DB.PP = max(DB.PGlobal,u1);
/*
	if (DB.PGlobal == 0) {
		DB.PP = DB.PGlobal+1;
	} else {
		if (PMax >= ML)
			DB.PP = max(PMax-ML,DB.PGlobal);
		else
			DB.PP = DB.PGlobal;
	}
*/

	// Geometry
	PGs = 1;

	if (strstr(DB.MeshType,"ToBeCurved") != NULL)
		DB.AC = 1, DB.ExactGeom = 0;
	else
		DB.AC = 0, DB.ExactGeom = 1;

	// ToBeModified (likely included in .ctrl file)
	strcpy(Parametrization,"ArcLength");
	//strcpy(Parametrization,"RadialProjection");

	for (i = 0; i < NEC; i++)
		NodeTypeG[i] = malloc(STRLEN_MIN * sizeof **NodeTypeG); // keep

	strcpy(NodeTypeG[0],"GLL");
	strcpy(NodeTypeG[1],"AO");
	strcpy(NodeTypeG[2],"GLL");

	// Order dependent parameters
	for (P = 0; P <= PMax; P++) {
		// Sum Factorization
		SF_BE[P] = malloc(2 * sizeof **SF_BE); 
		for (i = 0; i < 2; i++)
			SF_BE[P][i] = calloc(2 , sizeof ***SF_BE);

		// TP
		if (d == 1 || (d == 2 && P > 9) || (d == 3 && P > 9)) SF_BE[P][0][0] = 1; // ToBeModified
		if (d == 1 || (d == 2 && P > 9) || (d == 3 && P > 9)) SF_BE[P][0][1] = 1; // ToBeModified

		// WEDGE
		if (d == 3 && P > 99) SF_BE[P][1][0] = 1; // ToBeModified
		if (d == 3 && P > 99) SF_BE[P][1][1] = 1; // ToBeModified

		if (SF_BE[P][0][0] || SF_BE[P][0][1] || SF_BE[P][1][0] || SF_BE[P][1][1])
			printf("    Using Sum Factorized Operators for P%d.\n",P);

		// Geometry
		PCs[P] = malloc(NEC * sizeof **PCs); // keep
		PCc[P] = malloc(NEC * sizeof **PCc); // keep
		PJs[P] = malloc(NEC * sizeof **PJs); // keep
		PJc[P] = malloc(NEC * sizeof **PJc); // keep

		// ToBeDeleted: These orders may not be sufficient for 3D. To be investigated.
		PGc[P]    = max(P+1,u2);
		PCs[P][0] = PGs;
		PCs[P][1] = max(PGs-1,u1);
		PCs[P][2] = PGs;             // ToBeModified
		PCc[P][0] = PGc[P];
		PCc[P][1] = max(PGc[P]-1,u1);
PCc[P][2] = PGc[P]-1;          // ToBeModified
		PJs[P][0] = PGs;
		PJs[P][1] = max(PGs-1,u1);
		PJs[P][2] = PGs;             // ToBeModified
		PJc[P][0] = PGc[P];
		PJc[P][1] = max(PGc[P]-1,u1);
		PJc[P][2] = PGc[P];          // ToBeModified

		// Integration & Interpolation
		PFrs[P] = malloc(NEC * sizeof **PFrs); // keep
		PFrc[P] = malloc(NEC * sizeof **PFrc); // keep
		PIfs[P] = malloc(NEC * sizeof **PIfs); // keep
		PIfc[P] = malloc(NEC * sizeof **PIfc); // keep
		PIvs[P] = malloc(NEC * sizeof **PIvs); // keep
		PIvc[P] = malloc(NEC * sizeof **PIvc); // keep

		NodeTypeS[P]   = malloc(NEC * sizeof **NodeTypeS); // keep
		NodeTypeF[P]   = malloc(NEC * sizeof **NodeTypeF); // keep
		NodeTypeFrs[P] = malloc(NEC * sizeof **NodeTypeFrs); // keep
		NodeTypeFrc[P] = malloc(NEC * sizeof **NodeTypeFrc); // keep
		NodeTypeIfs[P] = malloc(NEC * sizeof **NodeTypeIfs); // keep
		NodeTypeIfc[P] = malloc(NEC * sizeof **NodeTypeIfc); // keep
		NodeTypeIvs[P] = malloc(NEC * sizeof **NodeTypeIvs); // keep
		NodeTypeIvc[P] = malloc(NEC * sizeof **NodeTypeIvc); // keep

		for (i = 0; i < NEC; i++) {
			NodeTypeS[P][i]   = malloc(STRLEN_MIN * sizeof ***NodeTypeS);   // keep
			NodeTypeF[P][i]   = malloc(STRLEN_MIN * sizeof ***NodeTypeF);   // keep
			NodeTypeFrs[P][i] = malloc(STRLEN_MIN * sizeof ***NodeTypeFrs); // keep
			NodeTypeFrc[P][i] = malloc(STRLEN_MIN * sizeof ***NodeTypeFrc); // keep
			NodeTypeIfs[P][i] = malloc(STRLEN_MIN * sizeof ***NodeTypeIfs); // keep
			NodeTypeIfc[P][i] = malloc(STRLEN_MIN * sizeof ***NodeTypeIfc); // keep
			NodeTypeIvs[P][i] = malloc(STRLEN_MIN * sizeof ***NodeTypeIvs); // keep
			NodeTypeIvc[P][i] = malloc(STRLEN_MIN * sizeof ***NodeTypeIvc); // keep
		}

		if (!Collocated) {
			if (!EFE)
				PF[P] = P+1;
			else
				PF[P] = P;

			// Interpolation (TP)
			if (!EFE) {
				PFrs[P][0] = PCs[P][0]+PF[P];
				PFrc[P][0] = PCc[P][0]+PF[P];
			} else {
				PFrs[P][0] = P;
				PFrc[P][0] = P;
			}

			if (strstr(DB.NodeType,"GLL") != NULL) {
				if (P == 0)
					strcpy(NodeTypeS[P][0],"GL");
				else
					strcpy(NodeTypeS[P][0],"GLL");

				if (PF[P] == 0)
					strcpy(NodeTypeF[P][0],"GL");
				else
					strcpy(NodeTypeF[P][0],"GLL");

				strcpy(NodeTypeFrs[P][0],"GLL");
				strcpy(NodeTypeFrc[P][0],"GLL");
			} else {
				strcpy(NodeTypeS  [P][0],"GL");
				strcpy(NodeTypeF  [P][0],"GL");
				strcpy(NodeTypeFrs[P][0],"GL");
				strcpy(NodeTypeFrc[P][0],"GL");
			}

			// Interpolation (SI)
			if (!EFE) {
				PFrs[P][1] = PCs[P][1]+PF[P];
				PFrc[P][1] = PCc[P][1]+PF[P];
			} else {
				PFrs[P][1] = P;
				PFrc[P][1] = P;
			}

			if (strstr(DB.NodeType,"AO") != NULL) {
				if (P == 0) {
					strcpy(NodeTypeS[P][1],"WSH");
				} else {
					strcpy(NodeTypeS[P][1],"AO");
				}

				if (PF[P] == 0) {
					strcpy(NodeTypeF[P][1],"WSH");
				} else {
					strcpy(NodeTypeF[P][1],"AO");
				}

				strcpy(NodeTypeFrs[P][1],"AO");
				strcpy(NodeTypeFrc[P][1],"AO");
			} else {
				strcpy(NodeTypeS  [P][1],"WSH");
				strcpy(NodeTypeF  [P][1],"WSH");
				strcpy(NodeTypeFrs[P][1],"WSH");
				strcpy(NodeTypeFrc[P][1],"WSH");
			}

			// Interpolation (PYR)
			if (strstr(DB.NodeType,"GLL") != NULL) {
				if (P == 0)
					strcpy(NodeTypeS[P][2],"GL");
				else
					strcpy(NodeTypeS[P][2],"GLL");

				if (PF[P] == 0)
					strcpy(NodeTypeF[P][2],"GL");
				else
					strcpy(NodeTypeF[P][2],"GLL");

				strcpy(NodeTypeFrs[P][2],"GLL");
				strcpy(NodeTypeFrc[P][2],"GLL");
			} else {
				strcpy(NodeTypeS  [P][2],"GL");
				strcpy(NodeTypeF  [P][2],"GL");
				strcpy(NodeTypeFrs[P][2],"GL");
				strcpy(NodeTypeFrc[P][2],"GL");
			}

			// Integration
			IntOrderfs = max(2*P,u1);
			IntOrderfc = max(2*P,u1);
			IntOrdervs = max(2*P,u1);
			IntOrdervc = max(2*P,u1);

			// TP
			strcpy(NodeTypeIfs[P][0],"GL");
			strcpy(NodeTypeIfc[P][0],"GL");
			strcpy(NodeTypeIvs[P][0],"GL");
			strcpy(NodeTypeIvc[P][0],"GL");

			PIfs[P][0] = floor(1.0*IntOrderfs/2.0);
			PIfc[P][0] = floor(1.0*IntOrderfc/2.0);
			PIvs[P][0] = floor(1.0*IntOrdervs/2.0);
			PIvc[P][0] = floor(1.0*IntOrdervc/2.0);

			// SI
			if (d == 2) {
				strcpy(NodeTypeIfs[P][1],"GL");
				strcpy(NodeTypeIfc[P][1],"GL");
				strcpy(NodeTypeIvs[P][1],"WV");
				strcpy(NodeTypeIvc[P][1],"WV");

				PIfs[P][1] = floor(1.0*IntOrderfs/2.0);
				PIfc[P][1] = floor(1.0*IntOrderfc/2.0);
				PIvs[P][1] = IntOrdervs;
				PIvc[P][1] = IntOrdervc;
			} else if (d == 3) {
				strcpy(NodeTypeIfs[P][1],"WV");
				strcpy(NodeTypeIfc[P][1],"WV");
				strcpy(NodeTypeIvs[P][1],"WV");
				strcpy(NodeTypeIvc[P][1],"WV");

				PIfs[P][1] = IntOrderfs;
				PIfc[P][1] = IntOrderfc;
				PIvs[P][1] = IntOrdervs;
				PIvc[P][1] = IntOrdervc;
			}

			// PYR
			strcpy(NodeTypeIfs[P][2],"WV");
			strcpy(NodeTypeIfc[P][2],"WV");
			strcpy(NodeTypeIvs[P][2],"GLW");
			strcpy(NodeTypeIvc[P][2],"GLW");

			PIfs[P][2] = IntOrderfs;
PIfc[P][2] = IntOrderfc+1;
			PIvs[P][2] = floor(1.0*IntOrdervs/2.0);
PIvc[P][2] = floor(1.0*IntOrdervc/2.0+1.0);
		} else {
			// THESE PARAMETERS CANNOT BE MODIFIED.
			if (strstr(DB.BasisType,"Nodal") == NULL) {
				printf("Selected BasisType: %s.\n",DB.BasisType);
				printf("Error: Nodal BasisType must be selected to use a collocated scheme\n"), exit(1);
			}

			PF[P] = P;

			// TP
			PFrs[P][0] = P;
			PFrc[P][0] = P;

			if (strstr(DB.NodeType,"GLL") != NULL && P > 0) {
				/*
				 * Brian's original code parameters
				 * This is the only version of the scheme where all nodes are collocated (i.e. FACET nodes are
				 * collocated with VOLUME nodes as well). Consequently, there is negligible additional cost for the
				 * strong form as compared to the weak form as the discontinuous boundary flux term is already evaluated
				 * in the VOLUME term.
				 */

				strcpy(NodeTypeS  [P][0],"GLL");
				strcpy(NodeTypeF  [P][0],"GLL");
				strcpy(NodeTypeFrs[P][0],"GLL");
				strcpy(NodeTypeFrc[P][0],"GLL");

				strcpy(NodeTypeIfs[P][0],"GLL");
				strcpy(NodeTypeIfc[P][0],"GLL");
				strcpy(NodeTypeIvs[P][0],"GLL");
				strcpy(NodeTypeIvc[P][0],"GLL");
			} else {
				strcpy(NodeTypeS  [P][0],"GL");
				strcpy(NodeTypeF  [P][0],"GL");
				strcpy(NodeTypeFrs[P][0],"GL");
				strcpy(NodeTypeFrc[P][0],"GL");

				strcpy(NodeTypeIfs[P][0],"GL");
				strcpy(NodeTypeIfc[P][0],"GL");
				strcpy(NodeTypeIvs[P][0],"GL");
				strcpy(NodeTypeIvc[P][0],"GL");
			}

			// SI
			PFrs[P][1] = P;
			PFrc[P][1] = P;

			strcpy(NodeTypeS  [P][1],"WSH");
			strcpy(NodeTypeF  [P][1],"WSH");
			strcpy(NodeTypeFrs[P][1],"WSH");
			strcpy(NodeTypeFrc[P][1],"WSH");

			strcpy(NodeTypeIvs[P][1],"WSH");
			strcpy(NodeTypeIvc[P][1],"WSH");
			if (d == 2) {
				if (strstr(DB.NodeType,"GLL") != NULL && P > 0) {
					strcpy(NodeTypeIfs[P][1],"GLL");
					strcpy(NodeTypeIfc[P][1],"GLL");
				} else {
					strcpy(NodeTypeIfs[P][1],"GL");
					strcpy(NodeTypeIfc[P][1],"GL");
				}
			} else if (d == 3) {
				strcpy(NodeTypeIfs[P][1],"WV");
				strcpy(NodeTypeIfc[P][1],"WV");
			}

			// PYR
			PFrs[P][2] = P;
			PFrc[P][2] = P;

			// Note: Collocated interpolation/integration nodes for PYR elements are currently unavailable
			strcpy(NodeTypeS  [P][2],"NOT_SUPPORTED");
			strcpy(NodeTypeF  [P][2],"NOT_SUPPORTED");
			strcpy(NodeTypeFrs[P][2],"NOT_SUPPORTED");
			strcpy(NodeTypeFrc[P][2],"NOT_SUPPORTED");

			strcpy(NodeTypeIfs[P][2],"NOT_SUPPORTED");
			strcpy(NodeTypeIfc[P][2],"NOT_SUPPORTED");
			strcpy(NodeTypeIvs[P][2],"NOT_SUPPORTED");
			strcpy(NodeTypeIvc[P][2],"NOT_SUPPORTED");

			// For collocated interpolation and integration nodes, a desired VOLUME integration order cannot be
			// specified. Further, if using GLL NodeType, FACET integration order also cannot be specified.
			for (i = 0; i < NEC; i++) {
				PIfs[P][i] = P;
				PIfc[P][i] = P;
				PIvs[P][i] = P;
				PIvc[P][i] = P;
			}
		}
	}

	// Solver
//	DB.InviscidFluxType = FLUX_LF;
	DB.InviscidFluxType = FLUX_ROE;

 	DB.ExplicitSolverType = RK3_SSP;
// 	DB.ExplicitSolverType = RK4_LS;

	// hp adaptation
	DB.DOFcap_frac = 10.0;
	DB.refine_frac = 0.3;
	DB.coarse_frac = 0.1;

	// Assign DB Parameters
	DB.NP    = NP;
	DB.NEC   = NEC;

	DB.SF_BE = SF_BE;
	DB.PGs   = PGs;
	DB.PGc   = PGc;
	DB.PCs   = PCs;
	DB.PCc   = PCc;
	DB.PJs   = PJs;
	DB.PJc   = PJc;
	DB.PF    = PF;
	DB.PFrs  = PFrs;
	DB.PFrc  = PFrc;
	DB.PIfs  = PIfs;
	DB.PIfc  = PIfc;
	DB.PIvs  = PIvs;
	DB.PIvc  = PIvc;

	DB.Parametrization = Parametrization;
	DB.NodeTypeG       = NodeTypeG;
	DB.NodeTypeS       = NodeTypeS;
	DB.NodeTypeF       = NodeTypeF;
	DB.NodeTypeFrs     = NodeTypeFrs;
	DB.NodeTypeFrc     = NodeTypeFrc;
	DB.NodeTypeIfs     = NodeTypeIfs;
	DB.NodeTypeIfc     = NodeTypeIfc;
	DB.NodeTypeIvs     = NodeTypeIvs;
	DB.NodeTypeIvc     = NodeTypeIvc;
}
