// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_update_h.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "mkl.h"
#include <mpi.h>
#include <petscksp.h>

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "test_support.h"
#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "adaptation.h"
#include "array_norm.h"
#include "array_print.h"
#include "memory_free.h"
#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of update_VOLUME_hp and update_FACET_hp for h-adaptation. Also ensures that
 *		Jacobians of refined VOLUMEs are positive.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_Limits {
	unsigned int index; // 0: < (All dims), 1: > (All dims)
	char         type;  // (a)nd, (d)iagonal, (o)r
	double       XYZ[3];
};

static void code_startup(int nargc, char **argv, const unsigned int Nref);
static void code_cleanup(const unsigned int final);
static void mark_VOLUMEs(const unsigned int adapt_type, const struct S_Limits *XYZ_lim);
static void check_correspondence(unsigned int *pass);
static void check_Jacobians(unsigned int *pass);
static void run_test(unsigned int *pass, const char *test_type);

void test_integration_update_h(int nargc, char **argv)
{
	unsigned int pass = 0;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *
	 *		Mesh file to be h-refined/coarsened.
	 *
	 *	Expected Output:
	 *
	 *		Interpolation of VOLUME geometry nodes to maximally 1-irregular FACETs from both sides should match.
	 *
	 */

	struct S_Limits *XYZ_lim;

	XYZ_lim = malloc(sizeof *XYZ_lim); // free

	// **************************************************************************************************** //
	// LINEs


	// **************************************************************************************************** //
	// TRIs
	strcpy(argvNew[1],"test/Test_update_h_TRI");

	code_startup(nargc,argvNew,2);
	if (DB.PGlobal <= 1)
		printf("Please increase PGlobal above 1 in the ctrl file (%s.ctrl)\n",argvNew[1]), TestDB.Nwarnings++;

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (TRI,   FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	XYZ_lim->XYZ[0] = -0.75; XYZ_lim->XYZ[1] = -0.75; XYZ_lim->type = 'd'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->type = 'd'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	mesh_update();

	XYZ_lim->XYZ[0] = -0.50; XYZ_lim->XYZ[1] = -0.50; XYZ_lim->type = 'd'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->type = 'o'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	pass = 1;
	check_Jacobians(&pass);
	//     0         10        20        30        40        50
	printf("update_h (       Jacobians):                     ");
	test_print(pass);

	code_cleanup(0);


	// **************************************************************************************************** //
	// QUADs
	strcpy(argvNew[1],"test/Test_update_h_QUAD");

	code_startup(nargc,argvNew,3);
	if (DB.PGlobal <= 1)
		printf("Please increase PGlobal above 1 in the ctrl file (%s.ctrl)\n",argvNew[1]), TestDB.Nwarnings++;

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (QUAD,  FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	XYZ_lim->XYZ[0] = -0.50; XYZ_lim->XYZ[1] = -0.50; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->type = 'a'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	mesh_update();

	XYZ_lim->XYZ[0] = -0.25; XYZ_lim->XYZ[1] = -0.25; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->type = 'o'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	pass = 1;
	check_Jacobians(&pass);
	//     0         10        20        30        40        50
	printf("update_h (       Jacobians):                     ");
	test_print(pass);

	code_cleanup(0);


	// **************************************************************************************************** //
	// TETs
	strcpy(argvNew[1],"test/Test_update_h_TET");

	code_startup(nargc,argvNew,2);
	if (DB.PGlobal <= 1)
		printf("Please increase PGlobal above 1 in the ctrl file (%s.ctrl)\n",argvNew[1]), TestDB.Nwarnings++;

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (TET,   FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	XYZ_lim->XYZ[0] =  1.00; XYZ_lim->XYZ[1] =  1.00; XYZ_lim->XYZ[2] =  0.00; XYZ_lim->type = 'd'; XYZ_lim->index = 1;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] = -1.00; XYZ_lim->XYZ[1] = -1.00; XYZ_lim->XYZ[2] =  1.00; XYZ_lim->type = 'd'; XYZ_lim->index = 0;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	mesh_update();

	XYZ_lim->XYZ[0] =  1.00; XYZ_lim->XYZ[1] =  1.00; XYZ_lim->XYZ[2] = -0.50; XYZ_lim->type = 'd'; XYZ_lim->index = 1;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  1.00; XYZ_lim->XYZ[1] =  1.00; XYZ_lim->XYZ[2] = -1.00; XYZ_lim->type = 'd'; XYZ_lim->index = 0;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	pass = 1;
	check_Jacobians(&pass);
	//     0         10        20        30        40        50
	printf("update_h (       Jacobians):                     ");
	test_print(pass);

	code_cleanup(0);


	// **************************************************************************************************** //
	// HEXs
	strcpy(argvNew[1],"test/Test_update_h_HEX");

	code_startup(nargc,argvNew,3);
	if (DB.PGlobal <= 1)
		printf("Please increase PGlobal above 1 in the ctrl file (%s.ctrl)\n",argvNew[1]), TestDB.Nwarnings++;

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (HEX,   FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	XYZ_lim->XYZ[0] = -0.50; XYZ_lim->XYZ[1] = -0.50; XYZ_lim->XYZ[2] = -0.50; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->XYZ[2] =  0.00; XYZ_lim->type = 'a'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	mesh_update();

	XYZ_lim->XYZ[0] = -0.25; XYZ_lim->XYZ[1] = -0.25; XYZ_lim->XYZ[2] = -0.25; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->XYZ[2] =  0.00; XYZ_lim->type = 'o'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	pass = 1;
	check_Jacobians(&pass);
	//     0         10        20        30        40        50
	printf("update_h (       Jacobians):                     ");
	test_print(pass);

	code_cleanup(0);


	// **************************************************************************************************** //
	// WEDGEs
	strcpy(argvNew[1],"test/Test_update_h_WEDGE");

	code_startup(nargc,argvNew,3);
	if (DB.PGlobal <= 1)
		printf("Please increase PGlobal above 1 in the ctrl file (%s.ctrl)\n",argvNew[1]), TestDB.Nwarnings++;

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (WEDGE, FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	XYZ_lim->XYZ[0] = -0.50; XYZ_lim->XYZ[1] = -0.50; XYZ_lim->XYZ[2] = -0.50; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->XYZ[2] =  0.00; XYZ_lim->type = 'a'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	mesh_update();

	XYZ_lim->XYZ[0] = -0.25; XYZ_lim->XYZ[1] = -0.25; XYZ_lim->XYZ[2] = -0.25; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->XYZ[2] =  0.00; XYZ_lim->type = 'o'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	pass = 1;
	check_Jacobians(&pass);
	//     0         10        20        30        40        50
	printf("update_h (       Jacobians):                     ");
	test_print(pass);

	code_cleanup(0);


	// **************************************************************************************************** //
	// PYRs
	strcpy(argvNew[1],"test/Test_update_h_PYR");

	code_startup(nargc,argvNew,2);
	if (DB.PGlobal <= 1)
		printf("Please increase PGlobal above 1 in the ctrl file (%s.ctrl)\n",argvNew[1]), TestDB.Nwarnings++;

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (PYR,   FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	XYZ_lim->XYZ[0] = -0.50; XYZ_lim->XYZ[1] = -0.50; XYZ_lim->XYZ[2] = -0.50; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->XYZ[2] =  0.00; XYZ_lim->type = 'a'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	mesh_update();

	XYZ_lim->XYZ[0] = -0.25; XYZ_lim->XYZ[1] = -0.25; XYZ_lim->XYZ[2] = -0.25; XYZ_lim->type = 'a'; XYZ_lim->index = 0;
	mark_VOLUMEs(HREFINE,XYZ_lim);
	XYZ_lim->XYZ[0] =  0.00; XYZ_lim->XYZ[1] =  0.00; XYZ_lim->XYZ[2] =  0.00; XYZ_lim->type = 'o'; XYZ_lim->index = 1;
	mark_VOLUMEs(HCOARSE,XYZ_lim);

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	pass = 1;
	check_Jacobians(&pass);
	//     0         10        20        30        40        50
	printf("update_h (       Jacobians):                     ");
	test_print(pass);

	code_cleanup(0);

//output_to_paraview("ZTest_Geomadapt");
//output_to_paraview("ZTest_Normals");
//EXIT_MSG;

	free(argvNew[0]); free(argvNew[1]); free(argvNew);
	free(XYZ_lim);
}

static void code_startup(int nargc, char **argv, const unsigned int Nref)
{
	int  MPIrank, MPIsize;

	// Start MPI and PETSC
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	DB.MPIsize = MPIsize;
	DB.MPIrank = MPIrank;

	// Initialization
	initialization(nargc,argv);
	setup_parameters();
	setup_mesh();
	setup_operators();
	setup_structures();
	setup_geometry();

	initialize_test_case(Nref);
}

static void code_cleanup(const unsigned int final)
{
	mesh_to_level(0);
	memory_free();
	if (final)
		PetscFinalize();
}

static void mark_VOLUMEs(const unsigned int adapt_type, const struct S_Limits *XYZ_lim)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int i, dim, Nve, IndXYZ, update;
	double       XYZ_cent[3], *XYZ_vC, XYZ_lim_sum, XYZ_cent_sum;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		// Calculate centroid
		ELEMENT = get_ELEMENT_type(VOLUME->type);
		Nve = ELEMENT->Nve;

		XYZ_vC = VOLUME->XYZ_vC;

		for (i = 0; i < 3; i++)
			XYZ_cent[i] = 0.0;

		for (dim = 0; dim < d; dim++) {
			IndXYZ = dim*Nve;
			for (i = 0; i < Nve; i++)
				XYZ_cent[dim] += XYZ_vC[IndXYZ+i];
			XYZ_cent[dim] /= Nve;
		}

		// Mark VOLUME if centroid is within limits
		switch (XYZ_lim->type) {
		case 'a':
			update = 1;
			for (dim = 0; dim < d; dim++) {
				switch (XYZ_lim->index) {
				case 0: if (!(XYZ_cent[dim] < XYZ_lim->XYZ[dim])) update = 0; break;
				case 1: if (!(XYZ_cent[dim] > XYZ_lim->XYZ[dim])) update = 0; break;
				default:
					printf("Error: Unsupported index in mark_VOLUMEs for type (a).\n"), EXIT_MSG;
					break;
				}
			}
			break;
		case 'o':
			update = 0;
			for (dim = 0; dim < d; dim++) {
				switch (XYZ_lim->index) {
				case 0: if (XYZ_cent[dim] < XYZ_lim->XYZ[dim]) update = 1; break;
				case 1: if (XYZ_cent[dim] > XYZ_lim->XYZ[dim]) update = 1; break;
				default:
					printf("Error: Unsupported index in mark_VOLUMEs for type (o).\n"), EXIT_MSG;
					break;
				}
			}
			break;
		case 'd':
			update = 1;
			XYZ_lim_sum  = 0.0;
			XYZ_cent_sum = 0.0;
			for (dim = 0; dim < d; dim++) {
				XYZ_lim_sum  += XYZ_lim->XYZ[dim];
				XYZ_cent_sum += XYZ_cent[dim];
			}
			switch (XYZ_lim->index) {
			case 0: if (!(XYZ_cent_sum < XYZ_lim_sum)) update = 0; break;
			case 1: if (!(XYZ_cent_sum > XYZ_lim_sum)) update = 0; break;
			default:
				printf("Error: Unsupported index in mark_VOLUMEs for type (d).\n"), EXIT_MSG;
				break;
			}
			break;
		default:
			printf("Error: Unsupported type in mark_VOLUMEs.\n"), EXIT_MSG;
			break;
		}

		if (update) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = adapt_type;
		}
	}
}


struct S_OPERATORS {
	unsigned int NfnS, NvnGs, *nOrdInOut, *nOrdOutIn;
	double       **I_vGs_fS;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndFType)
{
	unsigned int PF, VType, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACET;

//	PV    = VOLUME->P;
	VType = VOLUME->type;

	PF = FACET->P;
	IndOrdInOut = FACET->IndOrdInOut;
	IndOrdOutIn = FACET->IndOrdOutIn;

	ELEMENT = get_ELEMENT_type(VType);
	ELEMENT_FACET = get_ELEMENT_FACET(VType,IndFType);

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NfnS  = ELEMENT_FACET->NvnS[PF];

	OPS->I_vGs_fS = ELEMENT->I_vGs_fS[1][PF];

	OPS->nOrdInOut = ELEMENT_FACET->nOrd_fS[PF][IndOrdInOut];
	OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fS[PF][IndOrdOutIn];
}

static void check_correspondence(unsigned int *pass)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int Vf, IndFType, NfnS, *nOrdInOut, *nOrdOutIn,
	             dim, n, Indd, BC, FACET_is_internal;
	double       *XYZ_fSIn, *XYZ_fSOut, *XYZ_fSInOut, *XYZ_fSOutIn;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;
	struct S_FACET     *FACET;

	OPS = malloc(sizeof *OPS); // free

	*pass = 1;
	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
		VOLUME = FACET->VIn;
		Vf     = FACET->VfIn;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACET,IndFType);

		NfnS = OPS->NfnS;

		nOrdInOut = OPS->nOrdInOut;
		nOrdOutIn = OPS->nOrdOutIn;

		XYZ_fSIn = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vC); // free

		VOLUME = FACET->VOut;
		Vf     = FACET->VfOut;

		IndFType = get_IndFType(VOLUME->Eclass,Vf/NFREFMAX);
		init_ops(OPS,VOLUME,FACET,IndFType);

		XYZ_fSOut = mm_Alloc_d(CBCM,CBT,CBNT,NfnS,d,OPS->NvnGs,1.0,OPS->I_vGs_fS[Vf],VOLUME->XYZ_vC); // free

		XYZ_fSInOut = malloc(NfnS*d * sizeof *XYZ_fSInOut); // free
		XYZ_fSOutIn = malloc(NfnS*d * sizeof *XYZ_fSOutIn); // free

		for (dim = 0; dim < d; dim++) {
			Indd = dim*NfnS;
			for (n = 0; n < NfnS; n++) {
				XYZ_fSInOut[Indd+n] = XYZ_fSIn[Indd+nOrdInOut[n]];
				XYZ_fSOutIn[Indd+n] = XYZ_fSOut[Indd+nOrdOutIn[n]];
			}
		}

		BC = FACET->BC;
		FACET_is_internal = (BC == 0 || (BC % BC_STEP_SC > 50));

		if (FACET_is_internal && (array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSOutIn,"Inf")  > 10*EPS ||
		                          array_norm_diff_d(NfnS*d,XYZ_fSInOut,XYZ_fSOut,"Inf") > 10*EPS)) {
				*pass = 0;
				printf("Problem in check_correspondence\n");
printf("%d %d %d\n",FACET->indexg,FACET->IndOrdInOut,FACET->IndOrdOutIn);
printf("%d %d %d %d\n",FACET->VIn->type,FACET->VIn->indexg,FACET->VfIn,FACET->VIn->level);
printf("%d %d %d %d\n",FACET->VOut->type,FACET->VOut->indexg,FACET->VfOut,FACET->VOut->level);
				printf("Errors: %e %e\n\n",array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSOutIn,"Inf"),
		                                   array_norm_diff_d(NfnS*d,XYZ_fSInOut,XYZ_fSOut,"Inf"));
				array_print_d(NfnS,d,XYZ_fSIn,'C');
				array_print_d(NfnS,d,XYZ_fSOutIn,'C');
				array_print_d(NfnS,d,XYZ_fSOut,'C');
				array_print_d(NfnS,d,XYZ_fSInOut,'C');
				break;
		}

		free(XYZ_fSIn);
		free(XYZ_fSOut);
		free(XYZ_fSInOut);
		free(XYZ_fSOutIn);
	}
	free(OPS);

	if (*pass)
		TestDB.Npass++;
}

static void check_Jacobians(unsigned int *pass)
{
	unsigned int i, P, NvnI;
	double       *detJV_vI;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;

	for (VOLUME = DB.VOLUME; *pass && VOLUME; VOLUME = VOLUME->next) {
		P = VOLUME->P;

		ELEMENT = get_ELEMENT_type(VOLUME->type);

		if (!VOLUME->curved)
			NvnI = ELEMENT->NvnIs[P];
		else
			NvnI = ELEMENT->NvnIc[P];

		detJV_vI = VOLUME->detJV_vI;
		for (i = 0; i < NvnI; i++) {
			if (detJV_vI[i] <= 0.0) {
				*pass = 0;
				break;
			}
		}
	}

	if (*pass)
		TestDB.Npass++;
}

static void run_test(unsigned int *pass, const char *test_type)
{
	struct S_VOLUME *VOLUME;

	if (strstr(test_type,"FullREFINE") != NULL) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HREFINE;
		}
	} else if (strstr(test_type,"FullCOARSE")) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HCOARSE;
		}
	} else {
		// VOLUME processing done outside of this function.
	}
	mesh_update();
	check_correspondence(pass);
}
