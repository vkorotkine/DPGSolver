// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "initialize_test_case.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"
#include "exact_solutions.h"
#include "adaptation.h"
#include "compute_errors.h"
#include "array_print.h"
#include "output_to_paraview.h"
#include "update_VOLUMEs.h"
#include "setup_geom_factors.h"
#include "array_free.h"

/*
 *	Purpose:
 *		Initialize solution on all VOLUMEs for each test case.
 *
 *	Comments:
 *		For the Poisson case, the solution is initialized to 0 such that the update in the linear solve returns the
 *		exact solution (i.e. du = u-u0 = u).
 *
 *		For unsteady cases, it would be advantageous to write an adaptive initialization function. This could most
 *		easily be done by updating the mesh to be in a much finer space than the original, initializing the solution
 *		there, and only coarsening in regions where the initial solution error is low (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NvnI, NvnG, NvnC, NvnSPm1;
	double       *w_vI, *I_vG_vI, *I_vG_vS, *ChiS_vI, *ChiInvS_vS, *I_vC_vS, **D_vG_vC, **GradChiS_vS,
	             **Ihat_vS_vS, **L2hat_vS_vS;
};

struct S_VInfo {
	int          P, level;
	unsigned int type, refine_p, refine_h, adapt_class, adapt_type,
	             fh_range[NFMAX*2], neigh[NFMAX*NFREFMAX];
	double       L2s;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME)
{
	// Standard datatypes
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT_OPS;

	P = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT_OPS = get_ELEMENT_type(type);

	OPS->NvnS        = ELEMENT_OPS->NvnS[P];
//	OPS->NvnSPm1     = ELEMENT_OPS->NvnS[P-1];
	OPS->ChiInvS_vS  = ELEMENT_OPS->ChiInvS_vS[P][P][0];
	OPS->GradChiS_vS = ELEMENT_OPS->GradChiS_vS[P][P][0];
//	OPS->Ihat_vS_vS  = ELEMENT_OPS->Ihat_vS_vS[P-1][P];
//	OPS->L2hat_vS_vS = ELEMENT_OPS->L2hat_vS_vS[P][P-1];
	if (!curved) {
		OPS->NvnI = ELEMENT_OPS->NvnIs[P];
		OPS->NvnG = ELEMENT_OPS->NvnGs[1];
		OPS->NvnC = ELEMENT_OPS->NvnCs[P];

		OPS->w_vI = ELEMENT_OPS->w_vIs[P];

		OPS->I_vG_vI = ELEMENT_OPS->I_vGs_vIs[1][P][0];
		OPS->I_vG_vS = ELEMENT_OPS->I_vGs_vS[1][P][0];
		OPS->I_vC_vS = ELEMENT_OPS->I_vCs_vS[P][P][0];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIs[P][P][0];

		OPS->D_vG_vC = ELEMENT_OPS->D_vGs_vCs[1][P][0];
	} else {
		OPS->NvnI = ELEMENT_OPS->NvnIc[P];
		OPS->NvnG = ELEMENT_OPS->NvnGc[P];
		OPS->NvnC = ELEMENT_OPS->NvnCc[P];

		OPS->w_vI = ELEMENT_OPS->w_vIc[P];

		OPS->I_vG_vI = ELEMENT_OPS->I_vGc_vIc[P][P][0];
		OPS->I_vG_vS = ELEMENT_OPS->I_vGc_vS[P][P][0];
		OPS->I_vC_vS = ELEMENT_OPS->I_vCc_vS[P][P][0];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIc[P][P][0];

		OPS->D_vG_vC = ELEMENT_OPS->D_vGc_vCc[P][P][0];
	}
}

void compute_solution(const unsigned int Nn, double *XYZ, double *UEx, const unsigned int solved);
static void adapt_initial(unsigned int *adapt_update);
static void check_levels_refine(const unsigned int indexg, struct S_VInfo **VInfo_list, const unsigned int adapt_class);

void initialize_test_case_parameters(void)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase,
	             *Geometry = DB.Geometry;
	unsigned int d         = DB.d;

	DB.Nvar = d+2; // Euler and NS Equations
	DB.Neq  = d+2;

	// Standard datatypes
	char         *SolverType;
	unsigned int SourcePresent;

	if (strstr(TestCase,"Poisson")) {
		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
		strcpy(SolverType,"Implicit");
		SourcePresent = 1;

		if (strstr(Geometry,"dm1-Spherical_Section")) {
			DB.rIn  = 0.5;
			DB.rOut = 1.0;
		} else if (strstr(Geometry,"Ellipsoidal_Section")) {
			// These parameters must be consistent with the mesh for "ToBeCurved" meshes
			DB.rIn  = 0.5;
			DB.rOut = 1.0;

			// These parameters must be consistent with the mesh for "Curved" meshes
			DB.aIn  = 0.50; DB.bIn  = 0.50; DB.aOut = 1.00;

//			DB.bOut = 1.00; // POISSON_SCALE = 0.5
//			DB.bOut = 2.00; // POISSON_SCALE = 0.25
			DB.bOut = 3.00; // POISSON_SCALE = 0.125

			DB.cIn  = 1.50*DB.rIn;
			DB.cOut = 1.50*DB.rOut;
		} else if (strstr(Geometry,"Ringleb")) {
			DB.Q0   = 0.5;
			DB.KMin = 0.7;
			DB.KMax = 1.5;
		} else if (strstr(Geometry,"HoldenRamp")) {
			DB.lHR = 1.0;
//			DB.rIn = 2.54161;   // L/C ~= 1
//			DB.rIn = 1.52613;   // L/C ~= 2
			DB.rIn = 0.363683;  // L/C ~= 10
//			DB.rIn = 0.0380061; // L/C ~= 100
//			DB.rIn = 0.0038178; // L/C ~= 1000
		} else if (strstr(Geometry,"GaussianBump")) {
			DB.GBa = 0.0625;
			DB.GBb = 0.0;
			DB.GBc = 0.2/pow(2.0,1.0);
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}

		DB.Nvar = 1;
		DB.Neq  = 1;
	} else if (strstr(TestCase,"InviscidChannel")) {
		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
		strcpy(SolverType,"Implicit");
		SourcePresent = 0;

		// Equivalent to choosing total pressure/temperature and back pressure
		DB.MInf   = 0.3;
		DB.rhoInf = 1.0;
		DB.pInf   = 1.0;
		DB.cInf   = sqrt(GAMMA*DB.pInf/DB.rhoInf);

		if (strstr(Geometry,"GaussianBump")) {
			unsigned int BumpFactor = 3;

			DB.GBb = 0.0;
//			DB.GBa = 0.0625;
//			DB.GBc = 0.2/pow(2.0,BumpFactor);
			DB.GBa = 0.0625/pow(2.0,BumpFactor);
			DB.GBc = sqrt(1.5*DB.GBa);
		} else if (strstr(Geometry,"NacaSymmetric")) {
			double r = 0.25;

			DB.NSc = 1.0;
			DB.NSt = DB.NSc*sqrt(r)/1.1019;
			DB.NS0 =  0.2969;
			DB.NS1 = -0.1260;
			DB.NS2 = -0.3516;
			DB.NS3 =  0.2843;
			DB.NS4 = -0.1036;
		} else if (strstr(Geometry,"EllipsoidalBump")) {
			DB.aIn = 0.5;
			DB.bIn = DB.aIn/1.0;
		} else if (strstr(Geometry,"JoukowskiSymmetric")) {
			double a, l, t;
			a = 1.0;
			l = a/2.25;
			t = PI;

			double l2, l3, cost2;
			l2    = pow(l,2.0);
			l3    = pow(l,3.0);
			cost2 = pow(cos(t),2.0);

			DB.JSxL = -2.0*a*(l3+(l3+2.0*l2+l)*cost2+l2-(2.0*l3+3.0*l2+2.0*l+1.0)*cos(t)+l)/
			          (2.0*l2-2.0*(l2+l)*cos(t)+2.0*l+1.0);
			DB.JSa = a;
			DB.JSl = l;
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
	} else if (strstr(TestCase,"PeriodicVortex")) {
		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
		strcpy(SolverType,"Explicit");
		SourcePresent = 0;

//		DB.Xc = -DB.PeriodL*0.05;
		DB.Xc =  0.0;
		DB.Yc =  0.0;
		DB.Rc =  0.2;
//		DB.PeriodFraction = 0.5;
		DB.PeriodFraction = 1.0;

		DB.MInf = 0.5;
//		DB.MInf = 0.0;
		DB.pInf = 1.0;
		DB.TInf = 1.0;
		DB.Rg   = 1.0;

		DB.Cscale = 0.1;

		DB.uInf   = DB.MInf*sqrt(GAMMA*DB.Rg*DB.TInf);
		DB.vInf   = 0.1*EPS;
		DB.wInf   = 0.1*EPS;
		DB.VInf   = sqrt(DB.uInf*DB.uInf+DB.vInf*DB.vInf+DB.wInf*DB.wInf);

		if (fabs(DB.VInf) < 10*EPS)
			DB.FinalTime = 0.3;
		else
			DB.FinalTime = DB.PeriodFraction*DB.PeriodL/DB.VInf;
	} else if (strstr(TestCase,"SupersonicVortex")) {
		// Standard datatypes
		double pIn, cIn;

		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
//		strcpy(SolverType,"Explicit");
		strcpy(SolverType,"Implicit");
		SourcePresent = 0;

		DB.rIn  = 1.0;
		DB.rOut = 1.384;

		DB.MIn = 2.25;

		DB.rhoIn = 1.0;
		pIn   = pow(DB.rhoIn,GAMMA)/GAMMA;

		cIn = sqrt(GAMMA*pIn/DB.rhoIn);
		DB.VIn = cIn*DB.MIn/DB.rIn;
	} else {
		printf("Error: Unsupported TestCase: %s.\n",TestCase), EXIT_MSG;
	}

	if (strstr(SolverType,"Implicit"))
		DB.FinalTime  = 1e10;

	DB.SolverType    = SolverType;
	DB.SourcePresent = SourcePresent;
}

static void compute_gradient_polynomial(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int i, j, row, col, dim, dim2, IndD, IndC, NvnS, NvnG, NvnC;
	double       *XYZ_vS, *XYZ, *J_vC, *J_vS, *C_vS, *detJV_vS, *q, **D, *Dxyz;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	init_ops(OPS,VOLUME);

	// Initialize Geometric Transformation Operators
	NvnS   = OPS->NvnS;
	XYZ_vS = malloc(NvnS*d * sizeof *XYZ_vS); // free

	XYZ = VOLUME->XYZ;
	mm_CTN_d(NvnS,d,VOLUME->NvnG,OPS->I_vG_vS,XYZ,XYZ_vS);

	NvnG = OPS->NvnG;
	NvnC = OPS->NvnC;
	J_vC = malloc(NvnC*d*d * sizeof *J_vC); // free

	for (row = 0; row < d; row++) {
	for (col = 0; col < d; col++) {
		mm_CTN_d(NvnC,1,NvnG,OPS->D_vG_vC[col],&XYZ[NvnG*row],&J_vC[NvnC*(d*row+col)]);
	}}

	J_vS = malloc(NvnS*d*d * sizeof *J_vS); // free
	C_vS = malloc(NvnS*d*d * sizeof *C_vS); // free
	mm_CTN_d(NvnS,d*d,NvnC,OPS->I_vC_vS,J_vC,J_vS);
	mm_CTN_d(NvnS,d*d,NvnC,OPS->I_vC_vS,VOLUME->C_vC,C_vS);
	free(J_vC);

	detJV_vS = malloc(NvnS * sizeof *detJV_vS);
	compute_detJV(NvnS,J_vS,detJV_vS);

	// Compute solution gradients
	D = malloc(d * sizeof *D); // free
	for (dim = 0; dim < d; dim++) {
		D[dim] = malloc(NvnS*NvnS * sizeof *D[dim]); // free
		for (i = 0; i < NvnS*NvnS; i++)
			D[dim][i] = OPS->GradChiS_vS[dim][i];
	}

	for (dim = 0; dim < d; dim++) {
		Dxyz = calloc(NvnS*NvnS , sizeof *Dxyz); // free
		for (dim2 = 0; dim2 < d; dim2++) {
			IndD = 0;
			IndC = (dim+dim2*d)*NvnS;
			for (i = 0; i < NvnS; i++) {
				for (j = 0; j < NvnS; j++) {
					Dxyz[IndD+j] += D[dim2][IndD+j]*C_vS[IndC+j];
				}
				IndD += NvnS;
			}
		}
		IndD = 0;
		for (i = 0; i < NvnS; i++) {
			for (j = 0; j < NvnS; j++)
				Dxyz[IndD+j] /= detJV_vS[i];
			IndD += NvnS;
		}

		q = malloc(NvnS * sizeof *q); // free
		mm_CTN_d(NvnS,1,NvnS,Dxyz,VOLUME->uhat,q);
		mm_CTN_d(NvnS,1,NvnS,OPS->ChiInvS_vS,q,VOLUME->qhat[dim]);
		free(q);
		free(Dxyz);
	}
	array_free2_d(d,D);

	free(J_vS);
	free(C_vS);

	free(OPS);
}

static void compute_gradient_L2proj(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	// Standard datatypes
	unsigned int n, dim, NvnS, NvnI, P;
	double       *w_vI, *XYZ_vI, *q, *detJV_vI, *q_inter, *qhat_inter, *FilterP;

	struct S_OPERATORS *OPS;
	struct S_ELEMENT   *ELEMENT;

	OPS = malloc(sizeof *OPS); // free

	if (DB.Adapt == ADAPT_0)
		printf("Error: Unsupported.\n"), EXIT_MSG;
	else if (VOLUME->P == 0)
		printf("Error: Increased order (i.e. use P > 0) required.\n"), EXIT_MSG;

	init_ops(OPS,VOLUME);

	NvnS = OPS->NvnS;
	NvnI = OPS->NvnI;
	w_vI = OPS->w_vI;

	XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
	mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

	q = malloc(NvnI*Nvar*d * sizeof *q); // free

	compute_exact_gradient(NvnI,XYZ_vI,q);
	free(XYZ_vI);

	detJV_vI = VOLUME->detJV_vI;
	for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NvnI; n++)
			q[NvnI*dim+n] *= w_vI[n]*detJV_vI[n];
	}

	q_inter = malloc(NvnS*Nvar*d * sizeof *q_inter); // free
	mm_d(CBCM,CBNT,CBNT,NvnS,Nvar*d,NvnI,1.0,0.0,OPS->ChiS_vI,q,q_inter);
	free(q);

	compute_inverse_mass(VOLUME);

	// Get L2 projection of order P
	qhat_inter = malloc(NvnS*Nvar*d * sizeof *qhat_inter); // free
	for (dim = 0; dim < d; dim++)
		mm_d(CBCM,CBNT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,VOLUME->MInv,&q_inter[NvnS*Nvar*dim],&qhat_inter[NvnS*Nvar*dim]);
	free(q_inter);

	if (VOLUME->MInv) {
		free(VOLUME->MInv);
		VOLUME->MInv = NULL;
	}

	// Filter order P terms to get q accurate to order P-1
	ELEMENT = get_ELEMENT_type(VOLUME->type);
	P = VOLUME->P;

	OPS->NvnSPm1     = ELEMENT->NvnS[P-1];
	OPS->Ihat_vS_vS  = ELEMENT->Ihat_vS_vS[P-1][P];
	OPS->L2hat_vS_vS = ELEMENT->L2hat_vS_vS[P][P-1];

	FilterP = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,OPS->NvnSPm1,1.0,OPS->Ihat_vS_vS[0],OPS->L2hat_vS_vS[0]);

	for (dim = 0; dim < d; dim++)
		mm_CTN_d(NvnS,Nvar,NvnS,FilterP,&qhat_inter[NvnS*Nvar*dim],VOLUME->qhat[dim]);
	free(qhat_inter);
	free(FilterP);

	free(OPS);
}

void initialize_test_case(const unsigned int adapt_update_MAX)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d,
	             Nvar      = DB.Nvar,
	             Adapt     = DB.Adapt;

	DB.OutputInterval = 1e3;

	// Standard datatypes
	unsigned int DOF0 = 0, PolyGradient = 0;
	unsigned int n, dim, NvnS, NvnI, adapt_update, adapt_count;
	double       *XYZ_vS, *XYZ_vI, *U, *W, *What, *Q, *u_inter, *w_vI, *detJV_vI;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	adapt_count = 0;
	adapt_update = 1;
	while (adapt_update) {
		adapt_update = 0;
		if (strstr(TestCase,"PeriodicVortex") ||
		    strstr(TestCase,"SupersonicVortex") ||
			strstr(TestCase,"InviscidChannel")) {

			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				init_ops(OPS,VOLUME);

				NvnS         = OPS->NvnS;
				VOLUME->NvnS = NvnS;

				free(VOLUME->What);
				free(VOLUME->RES);

				VOLUME->What = malloc(NvnS*Nvar * sizeof *(VOLUME->What)); // keep
				VOLUME->RES  = calloc(NvnS*Nvar , sizeof *(VOLUME->RES));  // keep
				What = VOLUME->What;

				XYZ_vS = malloc(NvnS*d    * sizeof *XYZ_vS); // free

				mm_CTN_d(NvnS,d,VOLUME->NvnG,OPS->I_vG_vS,VOLUME->XYZ,XYZ_vS);

				U = malloc(NvnS*NVAR3D * sizeof *U); // free
				W = malloc(NvnS*Nvar   * sizeof *W); // free

				compute_solution(NvnS,XYZ_vS,U,0);

				convert_variables(U,W,3,d,NvnS,1,'p','c');
				mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,W,What);

				free(XYZ_vS);
				free(U);
				free(W);
			}
		} else if (strstr(TestCase,"Poisson")) {
			// Initializing with the L2 projection (Update other TestCases above), ToBeModified
			adapt_count = adapt_update_MAX; // No need for updating

			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				init_ops(OPS,VOLUME);

				NvnS = OPS->NvnS;
				NvnI = OPS->NvnI;
				w_vI = OPS->w_vI;

				VOLUME->NvnS = NvnS;

				free(VOLUME->uhat);
				for (dim = 0; dim < d; dim++)
					free(VOLUME->qhat[dim]);

				VOLUME->uhat = calloc(NvnS*Nvar , sizeof *(VOLUME->uhat)); // keep
				for (dim = 0; dim < d; dim++)
					VOLUME->qhat[dim] = calloc(NvnS*Nvar , sizeof *(VOLUME->qhat[dim])); // keep

				XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
				mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

				U = malloc(NvnI*Nvar * sizeof *U); // free

				compute_exact_solution(NvnI,XYZ_vI,U,0);
				free(XYZ_vI);

				detJV_vI = VOLUME->detJV_vI;
				for (n = 0; n < NvnI; n++)
					U[n] *= w_vI[n]*detJV_vI[n];

				u_inter = malloc(NvnS*Nvar * sizeof *u_inter); // free
				mm_d(CBCM,CBNT,CBNT,NvnS,Nvar,NvnI,1.0,0.0,OPS->ChiS_vI,U,u_inter);
				free(U);

				compute_inverse_mass(VOLUME);
				mm_d(CBCM,CBNT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,VOLUME->MInv,u_inter,VOLUME->uhat);
				free(u_inter);
				if (VOLUME->MInv) {
					free(VOLUME->MInv);
					VOLUME->MInv = NULL;
				}

				if (PolyGradient) {
					compute_gradient_polynomial(VOLUME);
				} else {
					if (Adapt != ADAPT_0) {
						compute_gradient_L2proj(VOLUME);
					} else {
						Q = calloc(NvnI*Nvar*d , sizeof *Q); // free
						for (dim = 0; dim < d; dim++)
							mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,&Q[NvnS*Nvar*dim],VOLUME->qhat[dim]);
						free(Q);
					}
				}
			}
		} else {
			printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
		}

		if (adapt_count < adapt_update_MAX) {
			switch (Adapt) {
			default: // ADAPT_P, ADAPT_H, ADAPT_HP
				adapt_initial(&adapt_update);
				break;
			case ADAPT_0:
				// Exit while loop.
				break;
			}
		}
		adapt_count++;
	}
	free(OPS);

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
		DOF0 += VOLUME->NvnS;
	DB.DOF0 = DOF0;

	// Ensure that update_VOLUME_Ops will be executed
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
		VOLUME->update = 1;

	// Output initial solution to paraview
//	output_to_paraview("ZTest_Sol_Init");
}

static void compute_uniform_solution(const unsigned int Nn, double *U)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	// Standard datatypes
	unsigned int n;

	if (strstr(TestCase,"InviscidChannel")) {
		double rhoInf, pInf, MInf, cInf, uInf;

		rhoInf = DB.rhoInf;
		pInf   = DB.pInf;
		MInf   = DB.MInf;
		cInf   = DB.cInf;

		uInf = MInf*cInf;

		for (n = 0; n < Nn; n++) {
			U[0*Nn+n] = rhoInf;
			U[1*Nn+n] = uInf;
			U[2*Nn+n] = 0.0;
			U[3*Nn+n] = 0.0;
			U[4*Nn+n] = pInf;
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

void compute_solution(const unsigned int Nn, double *XYZ, double *UEx, const unsigned int solved)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	if (strstr(TestCase,"PeriodicVortex") || strstr(TestCase,"SupersonicVortex")) {
		compute_exact_solution(Nn,XYZ,UEx,solved);
	} else if (strstr(TestCase,"InviscidChannel")) {
		compute_uniform_solution(Nn,UEx);
	} else {
		printf("Error: Unsupported TestCase: %s.\n",TestCase), EXIT_MSG;
	}
}

static void adapt_initial(unsigned int *adapt_update)
{
	/*
	 *	Purpose:
	 *		Perform h/p-adaptation of the initial mesh such that the discretization is either at the maximum allowed
	 *		resolution or the L2 error of the initial solution is below REFINE_TOL.
	 *
	 *	Comments:
	 *		Simultaneous h and p refinement is allowed for VOLUMEs which obtain refinement indication based on
	 *		refinement propagation.
	 *
	 *		This function is similar to adapt_hp, but only employs refinement and uses the solution L2 error as the
	 *		indicator.
	 *		The function is structured in preparation for MPI where structs of limited parts of VOLUMEs will be passed
	 *		between processors.
	 */

	// Initialize DB Parameters
	unsigned int NVglobal  = DB.NVglobal,
	             Adapt     = DB.Adapt,
	             PMax      = DB.PMax,
	             LevelsMax = DB.LevelsMax;

	// Standard datatypes
	unsigned int i, f, Nf, fh, fhMin, fhMax, Vf, indexg,
	             adapt_type, p_allow, h_allow, dummy_ui;
	double       *L2Error2, dummy_d;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;
	struct S_VInfo   *VInfo, **VInfo_list;

	VInfo_list = malloc(NVglobal * sizeof *VInfo_list); // free

	L2Error2 = malloc((NVAR3D+1) * sizeof *L2Error2); // free

	// Initialize VInfo structs
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;

		VInfo = malloc(sizeof *VInfo); // free

		VInfo->P          = (int) VOLUME->P;
		VInfo->level      = (int) VOLUME->level;
		VInfo->type       = VOLUME->type;
		VInfo->refine_p   = 0;
		VInfo->refine_h   = 0;
		VInfo->adapt_type = ADAPT_0;

		ELEMENT = get_ELEMENT_type(VInfo->type);
		Nf = ELEMENT->Nf;

		for (f = 0; f < Nf; f++) {
			get_fh_range(VOLUME,f,&fhMin,&fhMax);
			VInfo->fh_range[f*2  ] = fhMin;
			VInfo->fh_range[f*2+1] = fhMax;
			for (fh = fhMin; fh <= fhMax; fh++) {
				Vf = f*NFREFMAX+fh;
				VInfo->neigh[Vf] = VOLUME->neigh[Vf];
			}
		}
		VInfo_list[indexg] = VInfo;
	}

	// Compute L2 Errors
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		compute_errors(VOLUME,L2Error2,&dummy_d,&dummy_ui,0);

		VInfo = VInfo_list[VOLUME->indexg];
		VInfo->L2s = sqrt(L2Error2[NVAR3D]);

		switch (Adapt) {
		default: // ADAPT_HP
			// h vs p indicator goes here (ToBeDeleted).
			break;
		case ADAPT_P: VInfo->adapt_class = ADAPT_P; break;
		case ADAPT_H: VInfo->adapt_class = ADAPT_H; break;
		}
	}

	// MPI communication goes here. (ToBeDeleted)
	for (i = 0; i < NVglobal; i++) {
		VInfo = VInfo_list[i];
		if (!VInfo)
			printf("Error: MPI not supported.\n"), EXIT_MSG;
	}

	// Mark VOLUMEs for refinement (No limits on refinement)
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;

		VInfo = VInfo_list[indexg];
		if (VInfo->L2s > REFINE_TOL)
			check_levels_refine(indexg,VInfo_list,VInfo->adapt_class);
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		VInfo = VInfo_list[indexg];
		if (strstr(DB.TestCase,"Test")) {
			*adapt_update = 1;
			VOLUME->Vadapt = 1;
			if (Adapt == ADAPT_P)
				VOLUME->adapt_type = PREFINE;
			else if (Adapt == ADAPT_H)
				VOLUME->adapt_type = HREFINE;
			else
				printf("Error: Unsupported Adapt = %d.\n",Adapt), EXIT_MSG;
		} else {
			if (VInfo->refine_p || VInfo->refine_h) {
				adapt_type = VInfo->adapt_type;
				switch (Adapt) {
				default: // ADAPT_HP
					p_allow = (unsigned int) VInfo->P < PMax;
					h_allow = (unsigned int) VInfo->level < LevelsMax;
					break;
				case ADAPT_P:
					p_allow = (unsigned int) VInfo->P < PMax;
					h_allow = 0;
					break;
				case ADAPT_H:
					p_allow = 0;
					h_allow = (unsigned int) VInfo->level < LevelsMax;
					break;
				}
				if (p_allow || h_allow) {
					*adapt_update = 1;
					VOLUME->Vadapt = 1;
					if (adapt_type == HPREFINE) {
						if (p_allow && h_allow)
							VOLUME->adapt_type = VInfo->adapt_type;
						else if (p_allow)
							VOLUME->adapt_type = PREFINE;
						else
							VOLUME->adapt_type = HREFINE;
					} else {
						VOLUME->adapt_type = VInfo->adapt_type;
					}
				}
			}
		}
		free(VInfo);
	}

	free(VInfo_list);
	free(L2Error2);

	mesh_update();
}

static void check_levels_refine(const unsigned int indexg, struct S_VInfo **VInfo_list, const unsigned int adapt_class)
{
	unsigned int f, Nf, fh, fhMin, fhMax, indexg_n;

	struct S_ELEMENT *ELEMENT;
	struct S_VInfo   *VInfo, *VInfo_n;

	VInfo = VInfo_list[indexg];

	ELEMENT = get_ELEMENT_type(VInfo->type);
	Nf = ELEMENT->Nf;

	switch (adapt_class) {
	default: // ADAPT_P
		VInfo->refine_p   = 1;
		VInfo->adapt_type = PREFINE;
		if (VInfo->refine_h)
			VInfo->adapt_type = HPREFINE;

		for (f = 0; f < Nf; f++) {
			indexg_n = VInfo->neigh[f*NFREFMAX];
			VInfo_n  = VInfo_list[indexg_n];
			if (!(VInfo_n->refine_p) && (VInfo->P - VInfo_n->P) > 0)
				check_levels_refine(indexg_n,VInfo_list,adapt_class);
		}
		break;
	case ADAPT_H:
		VInfo->refine_h   = 1;
		VInfo->adapt_type = HREFINE;
		if (VInfo->refine_p)
			VInfo->adapt_type = HPREFINE;

		for (f = 0; f < Nf; f++) {
			fhMin = VInfo->fh_range[f*2  ];
			fhMax = VInfo->fh_range[f*2+1];
			for (fh = fhMin; fh <= fhMax; fh++) {
				indexg_n = VInfo->neigh[f*NFREFMAX+fh];
				VInfo_n  = VInfo_list[indexg_n];
				if (!(VInfo_n->refine_h) && (VInfo->level - VInfo_n->level) > 0)
					check_levels_refine(indexg_n,VInfo_list,adapt_class);
			}
		}
		break;
	}
}
