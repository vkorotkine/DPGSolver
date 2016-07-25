// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <petscksp.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"
#include "database.h"

/*
 *	Purpose:
 *		Test correctness of implementation of update_VOLUME_hp and update_FACET_hp for h-adaptation.
 *
 *	Comments:
 *		Ensure to test with BCs other than periodic as well. (ToBeDeleted)
 *		Fix memory leaks. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void code_startup(int nargc, char **argv, const unsigned int Nref);
static void code_cleanup(const unsigned int final);
static void check_correspondence(unsigned int *pass);
static void run_test(unsigned int *pass, const char *test_type);

void test_imp_update_h(int nargc, char **argv)
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


	unsigned int indexg;
	struct S_VOLUME *VOLUME;

	// LINEs


	// TRIs
	strcpy(argvNew[1],"test/Test_update_h_TRI");

	code_startup(nargc,argvNew,2);

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (TRI,   FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		if (indexg == 0) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HREFINE;
		} else if (indexg >= 16) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HCOARSE;
		}
	}
	mesh_update();

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		if (indexg <= 6) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HREFINE;
		} else if ((indexg >= 7 && indexg <= 14) || indexg >= 19) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HCOARSE;
		}
	}

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	mesh_to_level(2);
	//     0         10        20        30        40        50
	run_test(&pass,"ToLevel2");
	printf("update_h (       ToLevel2):                      ");
	test_print(pass);

//	output_to_paraview("ZTest_Geom");
//	output_to_paraview("ZTest_Normals");
//	exit(1);

	code_cleanup(0);


	// QUADs
	strcpy(argvNew[1],"test/Test_update_h_QUAD");

	code_startup(nargc,argvNew,3);

	//     0         10        20        30        40        50
	run_test(&pass,"FullREFINE");
	printf("update_h (QUAD,  FullREFINE):                    ");
	test_print(pass);

	//     0         10        20        30        40        50
	run_test(&pass,"FullCOARSE");
	printf("update_h (       FullCOARSE):                    ");
	test_print(pass);

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		if (indexg <= 3) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HREFINE;
		} else if (indexg >= 48) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HCOARSE;
		}
	}
	mesh_update();

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		if (indexg <= 23) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HREFINE;
		} else if ((indexg >= 32 && indexg <= 43) || indexg >= 48) {
			VOLUME->Vadapt = 1;
			VOLUME->adapt_type = HCOARSE;
		}
	}

	//     0         10        20        30        40        50
	run_test(&pass,"Mixed");
	printf("update_h (       Mixed):                         ");
	test_print(pass);

	mesh_to_level(2);
	//     0         10        20        30        40        50
	run_test(&pass,"ToLevel2");
	printf("update_h (       ToLevel2):                      ");
	test_print(pass);

//	output_to_paraview("ZTest_Geom");
//	output_to_paraview("ZTest_Normals");
//	exit(1);

	code_cleanup(0);

	free(argvNew[0]); free(argvNew[1]); free(argvNew);
}

static void code_startup(int nargc, char **argv, const unsigned int Nref)
{
	int  MPIrank, MPIsize;

	// Start MPI and PETSC
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	// Test memory leaks only from Petsc and MPI using valgrind
	//PetscFinalize(), exit(1);

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

		if (FACET_is_internal && (array_norm_diff_d(NfnS*d,XYZ_fSIn,XYZ_fSOutIn,"Inf")  > EPS ||
		                          array_norm_diff_d(NfnS*d,XYZ_fSInOut,XYZ_fSOut,"Inf") > EPS)) {
				*pass = 0;
				printf("Problem in check_correspondence\n");
printf("%d\n",FACET->indexg);
printf("%d %d\n",FACET->VIn->indexg,FACET->VfIn);
printf("%d %d\n",FACET->VOut->indexg,FACET->VfOut);
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
