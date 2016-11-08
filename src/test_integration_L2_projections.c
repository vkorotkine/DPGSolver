// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_L2_projections.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "compute_errors.h"
#include "adaptation.h"
#include "array_norm.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Test correctness of implementation of L2 projection operators.
 *
 *	Comments:
 *		L2 projection from TET to PYR results in large error. The error is reduced as the TET/PYR VOLUME cubature orders
 *		are	increased, but it seems that an impractically high order is required for machine precision. Tests for TET6
 *		and PYR L2 projections thus show high errors. The evidence for this being the likely cause of the error is that
 *		the difference between L2 errors before and after projection is on the order of machine precision if only a
 *		single refinement level is used for TET6 and is higher as soon as multiple levels are employed.
 *
 *	Notation:
 *
 *	References:
 */

static double *get_L2err    (void);
static void   mark_VOLUMEs  (const unsigned int adapt_type);

void test_integration_L2_projections(int nargc, char **argv)
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
	// TRIs
	strcpy(argvNew[1],"test/Test_L2_proj_p_TRI");

	code_startup(nargc,argvNew,2,0);

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

	code_cleanup();


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_TRI");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	// **************************************************************************************************** //
	// QUADs
	strcpy(argvNew[1],"test/Test_L2_proj_p_QUAD");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 10*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (QUAD,  ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_QUAD");

	code_startup(nargc,argvNew,3,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	// **************************************************************************************************** //
	// TETs

	strcpy(argvNew[1],"test/Test_L2_proj_p_TET");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (TET,   ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_TET");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e4*EPS) {
		pass = 1, TestDB.Npass++;
	} else if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e-3) {
		pass = 1, TestDB.Npass++;
		printf("\nWarning: L2 projection test for P%d TETs passing with norm_diff = % .3e\n\n",
		       DB.PGlobal,array_norm_diff_d(1,L2err[0],L2err[1],"Inf"));
		TestDB.Nwarnings++;
	}

	//     0         10        20        30        40        50
	printf("               (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();

	// **************************************************************************************************** //
	// HEXs
	strcpy(argvNew[1],"test/Test_L2_proj_p_HEX");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (HEX,   ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_HEX");

	code_startup(nargc,argvNew,3,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	// **************************************************************************************************** //
	// WEDGEs
	strcpy(argvNew[1],"test/Test_L2_proj_p_WEDGE");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (WEDGE, ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_WEDGE");

	code_startup(nargc,argvNew,3,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	// **************************************************************************************************** //
	// PYRs
	strcpy(argvNew[1],"test/Test_L2_proj_p_PYR");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(PREFINE); mesh_update();
	mark_VOLUMEs(PCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e3*EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("               (PYR,   ADAPT_P):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();


	strcpy(argvNew[0],argv[0]);
	strcpy(argvNew[1],"test/Test_L2_proj_h_PYR");

	code_startup(nargc,argvNew,2,0);

	L2err[0] = get_L2err();
	mark_VOLUMEs(HREFINE); mesh_update();
	mark_VOLUMEs(HCOARSE); mesh_update();
	L2err[1] = get_L2err();

	pass = 0;
	if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e3*EPS) {
		pass = 1, TestDB.Npass++;
	} else if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e-4) {
		pass = 1, TestDB.Npass++;
		printf("\nWarning: L2 projection test for P%d PYRs passing with norm_diff = % .3e\n\n",
		       DB.PGlobal,array_norm_diff_d(1,L2err[0],L2err[1],"Inf"));
		TestDB.Nwarnings++;
	}

	//     0         10        20        30        40        50
	printf("               (       ADAPT_H):                 ");
	test_print(pass);
	free(L2err[0]), free(L2err[1]);

	code_cleanup();

	free(argvNew[0]); free(argvNew[1]); free(argvNew);
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
