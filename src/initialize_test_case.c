// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "initialize_test_case.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>

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
#include "update_FACEs.h"
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
	char         *TestCase      = DB.TestCase,
	             *PDE           = DB.PDE,
	             *PDESpecifier  = DB.PDESpecifier,
	             *Geometry      = DB.Geometry,
	             *GeomSpecifier = DB.GeomSpecifier;
	unsigned int d              = DB.d;

	// Initialize Geometry related parameters which are used for all TestCase types
	bool InitializedGeometry = 1;
	if (strstr(Geometry,"GaussianBump")) {
		DB.GBa = 0.0625;
		DB.GBb = 0.0;

		double BExp = 0.0;
		if      (strstr(GeomSpecifier,"/0/"))   BExp = 0.0;
		else if (strstr(GeomSpecifier,"/0-5/")) BExp = 0.5;
		else if (strstr(GeomSpecifier,"/1/"))   BExp = 1.0;
		else
			EXIT_UNSUPPORTED;

		DB.GBc = 0.2/pow(2.0,BExp);
	} else if (strstr(Geometry,"Ringleb")) {
		DB.Q0   = 0.5;
		DB.KMin = 0.7;
		DB.KMax = 1.5;
	} else if (strstr(Geometry,"NacaSymmetric")) {
		double r = 0.25;

		DB.NSc = 1.0;
		DB.NSt = DB.NSc*sqrt(r)/1.1019;
		DB.NS0 =  0.2969;
		DB.NS1 = -0.1260;
		DB.NS2 = -0.3516;
		DB.NS3 =  0.2843;
		DB.NS4 = -0.1036;
	} else if (strstr(Geometry,"JoukowskiSymmetric")) {
		double a, l, t;
		a = 1.0;

		double ratio = 0.0;
		if      (strstr(GeomSpecifier,"/1/"))    ratio = 1.0;
		else if (strstr(GeomSpecifier,"/2-25/")) ratio = 2.25;
		else
			EXIT_UNSUPPORTED;

//			l = a/1.00; // MInf = 0.20
//			l = a/2.25; // MInf = 0.25
		l = a/ratio;
		t = PI;

		double l2, l3, cost2;
		l2    = pow(l,2.0);
		l3    = pow(l,3.0);
		cost2 = pow(cos(t),2.0);

		DB.JSxL = -2.0*a*(l3+(l3+2.0*l2+l)*cost2+l2-(2.0*l3+3.0*l2+2.0*l+1.0)*cos(t)+l)/
		          (2.0*l2-2.0*(l2+l)*cos(t)+2.0*l+1.0);
		DB.JSa = a;
		DB.JSl = l;
	} else if (strstr(Geometry,"HoldenRamp")) {
		// Current not used.
		printf("Warning: Ensure correct implementation before using.\n"), EXIT_MSG;

		DB.lHR = 1.0;

		printf("Error: Update this to depend on GeomSpecifier.\n"), EXIT_MSG;
//		DB.rIn = 2.54161;   // L/C ~= 1
//		DB.rIn = 1.52613;   // L/C ~= 2
		DB.rIn = 0.363683;  // L/C ~= 10
//		DB.rIn = 0.0380061; // L/C ~= 100
//		DB.rIn = 0.0038178; // L/C ~= 1000
	} else {
		if (strstr(DB.MeshType,"Curved"))
			InitializedGeometry = 0;
	}

	DB.Nvar = 0;
	DB.Neq  = 0;

	DB.SolverType = malloc(STRLEN_MIN * sizeof *(DB.SolverType)); // keep
	if (strstr(PDE,"Advection")) {
		// Currently requires div (dot) b = 0
		DB.PDE_index        = PDE_ADVECTION;
		DB.InviscidFluxType = FLUX_UPWIND;

		DB.Nvar = 1;
		DB.Neq  = 1;

		DB.Viscous = 0;
		DB.SourcePresent = 0;

		if (strstr(PDESpecifier,"Steady"))
			strcpy(DB.SolverType,"Implicit");
		else if (strstr(PDESpecifier,"Unsteady"))
			strcpy(DB.SolverType,"Explicit");
		else
			EXIT_UNSUPPORTED;

		if (strstr(PDESpecifier,"Default") || strstr(PDESpecifier,"Peterson")) {
			DB.ADV_b[0] = 0.0;
			DB.ADV_b[1] = 1.0;
			DB.ADV_b[2] = 0.0;
		} else {
			EXIT_UNSUPPORTED;
		}

		if (strstr(Geometry,"n-Cube")) {
			if (strstr(GeomSpecifier,"YL")) {
				DB.ADV_XYZB[0] = -1.0;
				DB.ADV_XYZB[1] = -1.0;
				DB.ADV_XYZB[2] = -1.0;
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(PDE,"Poisson")) {
		DB.PDE_index = PDE_POISSON;

		DB.Nvar = 1;
		DB.Neq  = 1;

		DB.Viscous = 1;

		DB.SourcePresent = 1;
		strcpy(DB.SolverType,"Implicit");

		// Default
		DB.Poisson_scale = 0.5;

		if (!InitializedGeometry) {
			if (strstr(Geometry,"n-Ball_HollowSection")) {
				DB.rIn  = 0.5;
				DB.rOut = 1.0;
			} else if (strstr(Geometry,"n-Ellipsoid")) {
				// These parameters must be consistent with the mesh for "ToBeCurved" meshes
				DB.rIn  = 0.5;
				DB.rOut = 1.0;

				DB.aIn  = 0.50;
				DB.bIn  = 0.50;
				DB.aOut = 1.00;

				if      (strstr(GeomSpecifier,"AR_1")) {
					DB.bOut = 1.0;
					DB.Poisson_scale = 0.5;
				} else if (strstr(GeomSpecifier,"AR_2")) {
					DB.bOut = 2.0;
					DB.Poisson_scale = 0.25;
				} else if (strstr(GeomSpecifier,"AR_3")) {
					DB.bOut = 3.0;
					DB.Poisson_scale = 0.125;
				} else {
					EXIT_UNSUPPORTED;
				}

				DB.cIn  = 1.50*DB.rIn;
				DB.cOut = 1.50*DB.rOut;
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	} else if (strstr(PDE,"Euler")) {
		DB.PDE_index = PDE_EULER;

		DB.Nvar = d+2;
		DB.Neq  = d+2;

		DB.Viscous = 0;
		DB.SourcePresent = 0;

		if (strstr(PDESpecifier,"Periodic")) {
			if (strstr(TestCase,"Vortex")) {
				strcpy(DB.SolverType,"Explicit");

				DB.Xc =  0.0;
				DB.Yc =  0.0;
				DB.Rc =  0.1;

				DB.pInf = 1.0;
				DB.TInf = 1.0;
				DB.Rg   = 1.0;

				DB.MInf    = 0.5;
				DB.uInf    = DB.MInf*sqrt(GAMMA*DB.Rg*DB.TInf);
				DB.vInf    = 1e-1*EPS;
				DB.wInf    = 1e-1*EPS;
				DB.VInf    = sqrt(DB.uInf*DB.uInf+DB.vInf*DB.vInf+DB.wInf*DB.wInf);
				DB.PeriodL = 2.0;

				DB.Cscale = 0.1*DB.VInf;

//				DB.PeriodFraction = 1.0;
				DB.PeriodFraction = 0.1;
				DB.FinalTime      = DB.PeriodFraction*DB.PeriodL/DB.VInf;

				// Update values for the stationary case
				if (strstr(TestCase,"Stationary")) {
					DB.MInf = 0.0;
					DB.uInf = DB.MInf*sqrt(GAMMA*DB.Rg*DB.TInf);
					DB.VInf = sqrt(DB.uInf*DB.uInf+DB.vInf*DB.vInf+DB.wInf*DB.wInf);

					DB.PeriodFraction = 0.0;
					// Uses FinalTime from moving vortex case.
				}
			} else {
				printf("%s\n",TestCase);
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(PDESpecifier,"Internal")) {
			if (!InitializedGeometry) {
				if (strstr(Geometry,"n-Cylinder")) {
					DB.rIn  = 1.0;
					DB.rOut = 1.384;
				} else if (strstr(Geometry,"EllipsoidalSection")) {
					if (strstr(GeomSpecifier,"Annular")) {
						// These parameters must be consistent with the mesh for "ToBeCurved" meshes
						DB.rIn  = 0.50;
						DB.rOut = 1.00;

						// These parameters must be consistent with the mesh for "Curved" meshes
						DB.aIn  = 0.50; DB.bIn  = 0.50; DB.aOut = 1.00;

						if      (strstr(GeomSpecifier,"/1/")) DB.bOut = 1.0;
						else if (strstr(GeomSpecifier,"/2/")) DB.bOut = 2.0;
						else if (strstr(GeomSpecifier,"/3/")) DB.bOut = 3.0;
						else
							EXIT_UNSUPPORTED;
					} else {
						DB.aIn = 0.5;
						double ratio = 0.0;
						if (strstr(GeomSpecifier,"/3/")) ratio = 3.0;
						else
							EXIT_UNSUPPORTED;

						DB.bIn = DB.aIn/ratio;
						DB.cIn = DB.aIn;
					}
				} else if (strstr(Geometry,"n-Cube")) {
					if (!(strstr(TestCase,"EllipticPipe") ||
					      strstr(TestCase,"ParabolicPipe") ||
					      strstr(TestCase,"SinusoidalPipe"))){

						double a_1 = 0.5, a_2 = 2.5, b = 0.5;
						DB.geo_store[0] = a_1;
						DB.geo_store[1] = a_2;
						DB.geo_store[2] = b;
						}

				} else {
					printf("%s\n",Geometry);
					EXIT_UNSUPPORTED;
				}
			}

			strcpy(DB.SolverType,"Implicit");
			if (strstr(TestCase,"SupersonicVortex")) {
				DB.MIn   = 2.25;
				DB.rhoIn = 1.0;

				double pIn, cIn;
				pIn = pow(DB.rhoIn,GAMMA)/GAMMA;
				cIn = sqrt(GAMMA*pIn/DB.rhoIn);

				DB.VIn = cIn*DB.MIn/DB.rIn;

//				printf("Add documentation.\n"), EXIT_BASIC;
			} else if (strstr(TestCase,"InviscidChannel")) {
				if (strstr(PDESpecifier,"Supersonic")) {
					DB.rhoInf = 1.0;
					DB.pInf   = 1.0;
					DB.MInf   = 1.01;
					DB.cInf   = sqrt(GAMMA*DB.pInf/DB.rhoInf);
				} else if (strstr(PDESpecifier,"Subsonic")) {
					DB.MInf    = 0.0;
					DB.p_Total = 1.0;
					DB.T_Total = 1.0;
					DB.Rg      = 1.0;
					DB.pBack   = 0.99*DB.p_Total;

					DB.rhoInf = DB.p_Total/(DB.Rg*DB.T_Total);
					DB.pInf   = DB.p_Total;

					DB.MInf   = 0.0*sqrt(2.0/GM1*(pow((DB.pBack/DB.p_Total),-GM1/GAMMA)-1.0));
					DB.cInf   = sqrt(GAMMA*DB.pInf/DB.rhoInf);
				} else {
					EXIT_UNSUPPORTED;
				}
			} else if (strstr(TestCase,"EllipticPipe")) {
				DB.SourcePresent = 1;
				        int i;
                                        double r_par[5] = {5, 2, 1, 1, 1}, p_par[5] = {1000, 250, 1, 250, 1};
					for (i = 0; i < 5; i++) {
      						DB.rho_store[i] = r_par[i];
						DB.p_store[i] = p_par[i];}

			} else if (strstr(TestCase,"ParabolicPipe")) {
				DB.SourcePresent = 1;

		      /*double r_par[7] = {5, -1, PI/2, 1, PI/5, 0, 0.25},
				       p_par[7] = {5, 1, PI/2, 1, PI/5, 0, 0.25},
				       u_par[7] = {3, -1, PI/2, -1, PI/5, 0, 0.25},
				       v_par[7] = {3, -1, PI/2, -1, PI/5, 0, 0.25};*/

			    double r_par[7] = {2, 0, 0.5, 0, PI/2, 0, PI/2},
				       p_par[7] = {2, 0, 0.5, 0, PI/2, 0, PI/2},
				       u_par[7] = {2, 0, 0.5, 0, PI/2, 0, PI/2},
				       v_par[7] = {2, 0, 0.5, 0, PI/2, 0, PI/2};

				for (int i = 0; i < 7; i++) {
					DB.rho_store[i] = r_par[i];
					DB.p_store[i]   = p_par[i];
					DB.u_store[i]   = u_par[i];
					DB.v_store[i]   = v_par[i];
				}

			} else if (strstr(TestCase,"SinusoidalPipe")) {
				DB.SourcePresent = 1;
                                        int i;
                                        double r_par[5] = {3, 1, 1, -1, 1}, p_par[5] = {6000, 2000, 1, -2000, 1}, w_par[5] = {30, 10, 0.5, 10, 0.5};
                                        for (i = 0; i < 5; i++) {
                                                DB.rho_store[i] = r_par[i];
                                                DB.p_store[i] = p_par[i];
                                                DB.w_store[i] = w_par[i];}
                        } else {
				printf("%s\n",TestCase);
				EXIT_UNSUPPORTED;
			}
		} else if (strstr(PDESpecifier,"External")) {
			if (strstr(TestCase,"PrandtlMeyer")) {
				strcpy(DB.SolverType,"Implicit");

				double l = 1.0, cIn;

				DB.aIn = 2.0*l;
				DB.bIn = l;

				// Compute at a reasonable angle (15 degrees) and then a 90 degree angle in preparation for ellipse.
				// ToBeDeleted
				printf("Error: Add dependence on GeomSpecifier.\n"), EXIT_MSG;
				DB.MIn   = 1.41421356237;
				DB.rhoIn = 1.0;
				DB.pIn   = 1.0;

				cIn = sqrt(GAMMA*DB.pIn/DB.rhoIn);
				DB.VIn = cIn*DB.MIn;

				// See: http://www.potto.org/fluidMech/2Dgd2.php
//				printf("Add documentation.\n"), EXIT_BASIC;
			} else {
				DB.MInf   = 0.10;
				DB.rhoInf = 1.0;
				DB.pInf   = 1.0;
				DB.cInf   = sqrt(GAMMA*DB.pInf/DB.rhoInf);
// ToBeDeleted
//// Initialized for Testing
//DB.p_Total = 1.0;
//DB.T_Total = 1.0;
//DB.Rg      = 1.0;
//DB.pBack   = 0.99*DB.p_Total;
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(PDE,"NavierStokes")) {
		DB.PDE_index = PDE_NAVIERSTOKES;

		DB.Nvar = d+2;
		DB.Neq  = d+2;

		DB.Viscous = 1;
		DB.Pr      = 0.72;
		DB.Rg      = 1.0;

		DB.Const_mu = 0; // Default

		DB.SourcePresent = 0;
		if (strstr(PDESpecifier,"Internal")) {
			if (!InitializedGeometry) {
				if (strstr(Geometry,"n-Cylinder")) {
					DB.rIn  = 0.5;
					DB.rOut = 1.0;
				} else {
					printf("%s\n",Geometry); EXIT_UNSUPPORTED;
				}
			}

			strcpy(DB.SolverType,"Explicit"); DB.FinalTime = 1e10; DB.ExplicitSolverType = EULER;
			strcpy(DB.SolverType,"Implicit");
			if (strstr(TestCase,"TaylorCouette")) {
				DB.omega = 1.0;
				DB.TIn   = 1.0;
				DB.pIn   = 1.0;
				DB.rhoIn = DB.pIn/(DB.Rg*DB.TIn);
				DB.mu    = 1e-3;

				DB.Const_mu = 1;
			} else if (strstr(TestCase,"PlaneCouette")) {
				DB.uIn   = 1.0;
				DB.TIn   = 1.0;
				DB.pIn   = 1.0;
				DB.rhoIn = DB.pIn/(DB.Rg*DB.TIn);
				DB.mu    = 1e-0;

				DB.Const_mu = 1;
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
		DB.Cp    = GAMMA/GM1*DB.Rg;
		DB.kappa = DB.mu*DB.Cp/DB.Pr;
	} else {
		printf("PDE: %s\n",PDE);
		EXIT_UNSUPPORTED;
	}

	if (strstr(DB.SolverType,"Implicit"))
		DB.FinalTime = 1e10;
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
		mm_CTN_d(NvnS,1,NvnS,Dxyz,VOLUME->What,q);
		mm_CTN_d(NvnS,1,NvnS,OPS->ChiInvS_vS,q,VOLUME->Qhat[dim]);
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
	/*
	 *	Comments:
	 *		Requires use of P (or HP) adaptation as this function performs a projection between different orders.
	 */

	if (!(DB.Adapt == ADAPT_P || DB.Adapt == ADAPT_HP))
		EXIT_UNSUPPORTED;

	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	// Standard datatypes
	unsigned int n, dim, NvnS, NvnI, P;
	double       *w_vI, *XYZ_vI, *q, *detJV_vI, *q_inter, *Qhat_inter, *FilterP;

	struct S_OPERATORS *OPS;
	struct S_ELEMENT   *ELEMENT;

	OPS = malloc(sizeof *OPS); // free
	init_ops(OPS,VOLUME);

	NvnS = OPS->NvnS;

	if (DB.Adapt == ADAPT_0)
		printf("Error: Unsupported.\n"), EXIT_MSG;
	else if (VOLUME->P == 0) {
		for (dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnS*Nvar; i++)
				VOLUME->Qhat[dim][i] = 0.0;
		}
		free(OPS);
		return;
	}

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
	Qhat_inter = malloc(NvnS*Nvar*d * sizeof *Qhat_inter); // free
	for (dim = 0; dim < d; dim++)
		mm_d(CBCM,CBNT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,VOLUME->MInv,&q_inter[NvnS*Nvar*dim],&Qhat_inter[NvnS*Nvar*dim]);
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
		mm_CTN_d(NvnS,Nvar,NvnS,FilterP,&Qhat_inter[NvnS*Nvar*dim],VOLUME->Qhat[dim]);
	free(Qhat_inter);
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

	DB.OutputInterval = 2e4;

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
		if (strstr(TestCase,"Advection") || strstr(TestCase,"Euler") || strstr(TestCase,"NavierStokes")) {
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				init_ops(OPS,VOLUME);

				NvnS         = OPS->NvnS;
				VOLUME->NvnS = NvnS;

				free(VOLUME->What);
				free(VOLUME->RES);

				VOLUME->What = malloc(NvnS*Nvar * sizeof *(VOLUME->What)); // keep
				VOLUME->RES  = calloc(NvnS*Nvar , sizeof *(VOLUME->RES));  // keep
				What = VOLUME->What;

				XYZ_vS = malloc(NvnS*d * sizeof *XYZ_vS); // free

				mm_CTN_d(NvnS,d,VOLUME->NvnG,OPS->I_vG_vS,VOLUME->XYZ,XYZ_vS);

				U = calloc(NvnS*NVAR3D , sizeof *U); // free
				W = malloc(NvnS*Nvar   * sizeof *W); // free

				if (strstr(TestCase,"Euler") || strstr(TestCase,"NavierStokes")) {
					compute_solution(NvnS,XYZ_vS,U,0);
					convert_variables(U,W,3,d,NvnS,1,'p','c');
				} else if (strstr(TestCase,"Advection")) {
					compute_solution(NvnS,XYZ_vS,W,0);
				}
				mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,W,What);

				free(XYZ_vS);
				free(U);
				free(W);
			}
		} else if (strstr(TestCase,"Poisson")) {
			// Initializing with the L2 projection.
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				init_ops(OPS,VOLUME);

				NvnS = OPS->NvnS;
				NvnI = OPS->NvnI;
				w_vI = OPS->w_vI;

				VOLUME->NvnS = NvnS;

				free(VOLUME->What);
				for (dim = 0; dim < d; dim++)
					free(VOLUME->Qhat[dim]);

				VOLUME->What = calloc(NvnS*Nvar , sizeof *(VOLUME->What)); // keep
				for (dim = 0; dim < d; dim++)
					VOLUME->Qhat[dim] = calloc(NvnS*Nvar , sizeof *(VOLUME->Qhat[dim])); // keep

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
				mm_d(CBCM,CBNT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,VOLUME->MInv,u_inter,VOLUME->What);
				free(u_inter);
				if (VOLUME->MInv) {
					free(VOLUME->MInv);
					VOLUME->MInv = NULL;
				}

				if (PolyGradient) {
					compute_gradient_polynomial(VOLUME);
				} else {
					if (Adapt == ADAPT_P || Adapt == ADAPT_HP) {
						compute_gradient_L2proj(VOLUME);
					} else {
						Q = calloc(NvnI*Nvar*d , sizeof *Q); // free
						for (dim = 0; dim < d; dim++)
							mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,&Q[NvnS*Nvar*dim],VOLUME->Qhat[dim]);
						free(Q);
					}
				}
			}
		} else {
			printf("Error: Unsupported TestCase (%s).\n",TestCase), EXIT_MSG;
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

	update_memory_VOLUMEs();
	update_memory_FACEs();
	// Output initial solution to paraview
//	output_to_paraview("ZTest_Sol_Init");
}

static void compute_uniform_solution(const unsigned int Nn, const double *XYZ, double *U)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int n;
	double       rhoInf, pInf, MInf, cInf, VInf;

	rhoInf = DB.rhoInf;
	pInf   = DB.pInf;
	MInf   = DB.MInf;
	cInf   = DB.cInf;

	VInf = MInf*cInf;
	if (strstr(TestCase,"InviscidChannel")) {
		for (n = 0; n < Nn; n++) {
			U[0*Nn+n] = rhoInf;
			U[1*Nn+n] = VInf;
			U[2*Nn+n] = 0.0;
			U[3*Nn+n] = 0.0;
			U[4*Nn+n] = pInf;
		}
	} else if (strstr(TestCase,"SubsonicNozzle")) {

		if (d == 3)
			printf("Add support.\n"), EXIT_MSG;

		if (strstr(DB.PDESpecifier,"Subsonic")) {
			for (n = 0; n < Nn; n++) {
				U[0*Nn+n] = DB.rhoInf;
				U[1*Nn+n] = 0.0;
				U[2*Nn+n] = 0.0;
				U[3*Nn+n] = 0.0;
				U[4*Nn+n] = DB.pInf;
			}
		} else if (strstr(DB.PDESpecifier,"Supersonic")) {
			// Define the initial solution such that the velocity vector points in approximately the correct direction.
			const double *X, *Y;

			X = &XYZ[0*Nn];
			Y = &XYZ[1*Nn];

			double aIn  = DB.aIn,
				   aOut = DB.aOut,
				   bIn  = DB.bIn,
				   bOut = DB.bOut;
//printf("\n");
			for (n = 0; n < Nn; n++) {
				// Find the equation of the ellipse on which the point lies
				// Note: Points on polynomial curved edges may lie outside of the analytical domain and will not be
				//       found by the algorithm below even for very high tolerances.
				double da, db, a, b, a_sign;

				db = bOut-bIn;

				da = 0.5*(aOut-aIn);
				a_sign = 1.0;
				a = aIn;

				unsigned int count, countMax = 60;
//printf("\n\n");
				for (count = 0; count < countMax; count++) {
					double y;
					a += a_sign*da;
					b = bIn+db*(a-aIn)/(aOut-aIn);

					da *= 0.5;

					if (X[n] > a) {
//printf("%2d % .3e % .3e % .3e % .3e\n",count,X[n],Y[n],a,X[n]-a);
						a_sign = 1.0;
						if (da > EPS)
							continue;

						if (fabs(X[n]-a) < 1e1*EPS)
							break;
//printf("daX\n");
					}

					y = b*sqrt(1.0-pow(X[n]/a,2.0));
//printf("cxy: %2d % .3e % .3e % .3e % .3e % .3e % .3e % .3e\n",count,X[n],Y[n],da,a,b,y,y-Y[n]);

					if (fabs(y-Y[n]) < 1e-5)
						break;

					if (da < EPS) {
						if (!(fabs(a-aIn) < 1e1*EPS || fabs(a-aOut) < 1e1*EPS))
							count = countMax;
						break;
					}

					if (Y[n] < y)
						a_sign = -1.0;
					else
						a_sign = 1.0;
				}

				if (count == countMax)
					printf("Error: Did not find ellipse (% .3e % .3e % .3e % .3e % .3e).\n",X[n],Y[n],da,a,b), EXIT_MSG;

				// Find the tangent to the ellipse in the direction of the flow
				double dydx, t1, t2, tNorm;

//			if (X[n] > a)
				if (Y[n] < EPS)
					dydx = -1e16;
				else
					dydx = 0.5*b/sqrt(1.0-pow(X[n]/a,2.0))*(-2.0*X[n]/(a*a));

				t1 = -1.0;
				t2 = -dydx;

				tNorm = sqrt(t1*t1+t2*t2);
				t1 /= tNorm;
				t2 /= tNorm;

				// Find the local Mach Number based on the area ratio relation
				double t, A, AIn, Achoke;

				t = atan2(Y[n],X[n]);

				A   = sqrt(pow((aOut-aIn)*cos(t),2.0)+pow((bOut-bIn)*sin(t),2.0));
				AIn = aOut-aIn;

				double MIn = DB.MInf;

				Achoke = AIn/(1.0/MIn*pow(2.0/(GAMMA+1)*(1+0.5*GM1*MIn*MIn),(GAMMA+1)/(2*GM1)));

				// Iterate to find M (Newton's method)
				double f, dfdM, M;

				M = MIn;
				for (count = 0; count < countMax; count++) {
					double update;
					f    = 1.0/M*pow(2.0/(GAMMA+1)*(1+0.5*GM1*M*M),(GAMMA+1)/(2*GM1)) - A/Achoke;
					dfdM = 2.0*(M*M-1.0)/((GAMMA-1)*pow(M,4.0)+2.0*M*M)*
						   pow(((GAMMA-1)*M*M+2.0)/(GAMMA+1),0.5*GAMMA/(GAMMA-1)+0.5/(GAMMA-1));

					update = -f/dfdM;
					if (fabs(update) < 1e2*EPS)
						break;

					M += update;
				}

				if (count == countMax)
					printf("Error: Newton's method not converging (% .3e % .3e % .3e % .3e).\n",f,dfdM,M,A/Achoke), EXIT_MSG;

//if (A/AIn > 1.0+EPS) {
//	printf("% .3e % .3e\n",A/AIn,M);
//	EXIT_MSG;
//}

				// Initialize the flow using Isentropic relations based on local Mach number
				double rhoIn = DB.rhoInf,
					   pIn   = DB.pInf;

				double rho, p, c, V;

				rho = rhoIn*pow((1+0.5*GM1*M*M)/(1+0.5*GM1*MIn*MIn),-1.0/GM1);
				p   = pIn*pow(rho/rhoIn,GAMMA);
				c   = sqrt(GAMMA*p/rho);
				V   = M*c;

				U[0*Nn+n] = rho;
				U[1*Nn+n] = V*t1;
				U[2*Nn+n] = V*t2;
				U[3*Nn+n] = 0.0;
				U[4*Nn+n] = p;
//printf("%2d % .3e % .3e % .3e % .3e % .3e\n",n,rho,V,t1,t2,p);
			}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
	} else if (strstr(TestCase,"PrandtlMeyer")) {
		// Standard datatypes
		const double *X, *Y;

		X = &XYZ[0*Nn];
		Y = &XYZ[1*Nn];
		if (d == 3)
			printf("Add support.\n"), EXIT_MSG;

		// Define the initial solution such that the velocity vector points in approximately the correct direction.
		for (n = 0; n < Nn; n++) {
			double t;

			t = atan2(Y[n],X[n]);
			U[0*Nn+n] = DB.rhoIn;
			U[4*Nn+n] = DB.pIn;
			if (X[n] < EPS) {
				U[1*Nn+n] = DB.VIn;
				U[2*Nn+n] = 0.0;
			} else if (Y[n] < EPS) {
				U[1*Nn+n] = 0.0;
				U[2*Nn+n] = -DB.VIn;
			} else {
				U[1*Nn+n] =  sin(t)*DB.VIn;
				U[2*Nn+n] = -cos(t)*DB.VIn;
			}
//			U[1*Nn+n] = 0.0;
//			U[2*Nn+n] = 0.0;
			U[3*Nn+n] = 0.0*t;
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

void compute_solution(const unsigned int Nn, double *XYZ, double *UEx, const unsigned int solved)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	if (strstr(TestCase,"Advection") ||
	    strstr(TestCase,"PeriodicVortex") ||
	    strstr(TestCase,"SupersonicVortex") ||
	    strstr(TestCase,"PlaneCouette") ||
	    strstr(TestCase,"TaylorCouette") ||
	    strstr(TestCase,"EllipticPipe") ||
	    strstr(TestCase,"ParabolicPipe") ||
	    strstr(TestCase,"SinusoidalPipe")) {
		compute_exact_solution(Nn,XYZ,UEx,solved);
	} else if (strstr(TestCase,"InviscidChannel") ||
	           strstr(TestCase,"PrandtlMeyer") ||
	           strstr(TestCase,"SubsonicNozzle")) {
		compute_uniform_solution(Nn,XYZ,UEx);
	} else {
		EXIT_UNSUPPORTED;
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
	L2Error2 = malloc((NVAR3D+1) * sizeof *L2Error2); // free
	if (strstr(DB.TestCase,"Euler")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			compute_errors(VOLUME,L2Error2,&dummy_d,&dummy_ui,0);

			VInfo = VInfo_list[VOLUME->indexg];
			VInfo->L2s = sqrt(L2Error2[NVAR3D]);

			switch (Adapt) {
			default: // ADAPT_HP
				VInfo->adapt_class = ADAPT_P; // ToBeModified
				// h vs p indicator goes here (ToBeDeleted).
				break;
			case ADAPT_P: VInfo->adapt_class = ADAPT_P; break;
			case ADAPT_H: VInfo->adapt_class = ADAPT_H; break;
			}
		}
	}

	// MPI communication goes here. (ToBeDeleted)
	for (i = 0; i < NVglobal; i++) {
		VInfo = VInfo_list[i];
		if (!VInfo)
			printf("Error: MPI not supported.\n"), EXIT_MSG;
	}

	// Mark VOLUMEs for refinement (No limits on refinement)
	if (strstr(DB.TestCase,"Euler")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			indexg = VOLUME->indexg;

			VInfo = VInfo_list[indexg];
			if (VInfo->L2s > REFINE_TOL)
				check_levels_refine(indexg,VInfo_list,VInfo->adapt_class);
		}
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
			else if (Adapt == ADAPT_HP) // Default to h-refinement for the time being. (ToBeModified)
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
