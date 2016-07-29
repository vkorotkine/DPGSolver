// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <petscksp.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"
#include "database.h"

/*
 *	Purpose:
 *		Test correctness of implementation of L2 projection operators.
 *
 *	Comments:
 *		Fix memory leaks. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void   code_startup  (int nargc, char **argv, const unsigned int Nref);
static void   code_cleanup  (const unsigned int final);
static double *get_L2err    (void);
static void   mark_VOLUMEs  (const unsigned int adapt_type);

void test_imp_L2_projections(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;
	double       *L2err[2];

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *
	 *		Mesh file to be refined/coarsened.
	 *
	 *	Expected Output:
	 *
	 *		The L2 error of the initial solution after refinement and coarsening should not change.
	 *
	 */

	// **************************************************************************************************** //
	// LINEs


	// **************************************************************************************************** //
	// TRIs
	strcpy(argvNew[1],"test/Test_L2_proj_p_TRI");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 10*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (TRI,   ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_TRI");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	// **************************************************************************************************** //
	// QUADs
	strcpy(argvNew[1],"test/Test_L2_proj_p_QUAD");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 10*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (QUAD,  ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_QUAD");

	code_startup(nargc,argvNew,3);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	// **************************************************************************************************** //
	// TETs
	strcpy(argvNew[1],"test/Test_L2_proj_p_TET");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (TET,   ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_TET");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	// **************************************************************************************************** //
	// HEXs
	strcpy(argvNew[1],"test/Test_L2_proj_p_HEX");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (HEX,   ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_HEX");

	code_startup(nargc,argvNew,3);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	// **************************************************************************************************** //
	// WEDGEs
	strcpy(argvNew[1],"test/Test_L2_proj_p_WEDGE");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (WEDGE, ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_WEDGE");

	code_startup(nargc,argvNew,3);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	// **************************************************************************************************** //
	// PYRs
	strcpy(argvNew[1],"test/Test_L2_proj_p_PYR");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (PYR,   ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup(0);


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_PYR");

	code_startup(nargc,argvNew,2);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("L2_projections (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

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
	memory_free();
	if (final)
		PetscFinalize();
}



struct S_OPERATORS {
	unsigned int NfnS, NvnGs, *nOrdInOut, *nOrdOutIn;
	double       **I_vGs_fS;
};

static double *get_L2err(void)
{
	unsigned int i, iMax, dummy_ui;
	double       dummy_d, *L2Error2, *L2Error;

	struct S_VOLUME *VOLUME;

	iMax = NVAR3D+1;
	L2Error2 = malloc(iMax * sizeof *L2Error2); // free
	L2Error  = calloc(iMax , sizeof *L2Error);  // keep (requires external free)
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		compute_errors(VOLUME,L2Error2,&dummy_d,&dummy_ui,0);

//		printf("%d\n",VOLUME->indexg);
//		array_print_d(1,NVAR3D+1,L2Error2,'R');

		for (i = 0; i < iMax; i++)
			L2Error[i] += sqrt(L2Error2[i]);
	}
//	array_print_d(1,iMax,L2Error,'R');

	free(L2Error2);
	return L2Error;
}

static void mark_VOLUMEs(const unsigned int adapt_type)
{
	struct S_VOLUME *VOLUME;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->Vadapt = 1;
		VOLUME->adapt_type = adapt_type;
	}
}
