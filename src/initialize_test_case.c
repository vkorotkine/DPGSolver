// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

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
	unsigned int NvnS;
	double       *I_vG_vS, *ChiInvS_vS;
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

	OPS->NvnS       = ELEMENT_OPS->NvnS[P];
	OPS->ChiInvS_vS = ELEMENT_OPS->ChiInvS_vS[P][P][0];
	if (!curved) {
		OPS->I_vG_vS = ELEMENT_OPS->I_vGs_vS[1][P][0];
	} else {
		OPS->I_vG_vS = ELEMENT_OPS->I_vGc_vS[P][P][0];
	}
}

static void compute_initial_solution(const unsigned int Nn, double *XYZ, double *UEx);
static void adapt_initial(unsigned int *adapt_update);
static void check_levels_refine(const unsigned int indexg, struct S_VInfo **VInfo_list, const unsigned int adapt_class);

void initialize_test_case_parameters(char *TestCase)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	DB.Nvar = d+2; // Euler and NS Equations
	DB.Neq  = d+2;

	// Standard datatypes
	char         *SolverType;
	unsigned int SourcePresent;

	if (strstr(TestCase,"Poisson")) {
		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
		strcpy(SolverType,"Implicit");
		SourcePresent = 1;

		DB.Nvar = 1;
		DB.Neq  = 1;
	} else if (strstr(TestCase,"dSphericalBump")) {
		EXIT_MSG;
	} else if (strstr(TestCase,"GaussianBump")) {
		EXIT_MSG;
	} else if (strstr(TestCase,"PeriodicVortex") ||
	           strstr(TestCase,"Test_L2_proj")   ||
	           strstr(TestCase,"Test_update_h")) {
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
	} else if (strstr(TestCase,"PolynomialBump")) {
		EXIT_MSG;
	} else if (strstr(TestCase,"SupersonicVortex") ||
	           strstr(TestCase,"Test_linearization")) {
		// Standard datatypes
		double pIn, cIn;

		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
//		strcpy(SolverType,"Explicit");
		strcpy(SolverType,"Implicit");
		SourcePresent = 0;

		DB.rIn  = 1.0;
//		rOut = 1.384;

		DB.MIn = 2.25;

		DB.rhoIn = 1.0;
		pIn   = pow(DB.rhoIn,GAMMA)/GAMMA;

		cIn = sqrt(GAMMA*pIn/DB.rhoIn);
		DB.VIn = cIn*DB.MIn/DB.rIn;

		DB.FinalTime  = 1e10;
	} else {
		printf("Error: Unsupported TestCase: %s.\n",TestCase), EXIT_MSG;
	}
	DB.SolverType    = SolverType;
	DB.SourcePresent = SourcePresent;
}

void initialize_test_case(const unsigned int adapt_update_MAX)
{
	initialize_test_case_parameters(DB.TestCase);

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d,
	             Nvar      = DB.Nvar,
	             Adapt     = DB.Adapt;

	DB.OutputInterval = 1e3;

	// Standard datatypes
	unsigned int DOF0 = 0;
	unsigned int NvnS, adapt_update, adapt_count;
	double       *XYZ_vS, *U, *W, *What;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	adapt_count = 0;
	adapt_update = 1;
	while (adapt_update) {
		adapt_update = 0;
		if (strstr(TestCase,"PeriodicVortex")   ||
		    strstr(TestCase,"SupersonicVortex") ||
		    strstr(TestCase,"Test")) {

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

				U   = malloc(NvnS*NVAR3D * sizeof *U); // free
				W   = malloc(NvnS*Nvar   * sizeof *W); // free

				compute_initial_solution(NvnS,XYZ_vS,U);

				convert_variables(U,W,3,d,NvnS,1,'p','c');
				mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,W,What);

				free(XYZ_vS);
				free(U);
				free(W);
			}
		} else if (strstr(TestCase,"Poisson")) {
			adapt_count = adapt_update_MAX; // No need for updating

			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				init_ops(OPS,VOLUME);

				NvnS         = OPS->NvnS;
				VOLUME->NvnS = NvnS;

				VOLUME->uhat = calloc(NvnS*Nvar , sizeof *(VOLUME->uhat)); // keep

 // ToBeDeleted
XYZ_vS = malloc(NvnS*d * sizeof *XYZ_vS); // free
mm_CTN_d(NvnS,d,VOLUME->NvnG,OPS->I_vG_vS,VOLUME->XYZ,XYZ_vS);

double *u = malloc(NvnS * sizeof *u); // free
compute_initial_solution(NvnS,XYZ_vS,u);
mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,u,VOLUME->uhat);
free(XYZ_vS);
free(u);
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

static void compute_initial_solution(const unsigned int Nn, double *XYZ, double *UEx)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	if (strstr(TestCase,"PeriodicVortex") || strstr(TestCase,"SupersonicVortex") || strstr(TestCase,"Test")) {
		compute_exact_solution(Nn,XYZ,UEx,0);
} else if (strstr(TestCase,"Poisson")) { // ToBeDeleted
	compute_exact_solution(Nn,XYZ,UEx,0);
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
