// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
 *		L2 projections between TET and PYR ELEMENTs result in significantly higher error than other projections.
 *		Observations:
 *			- This error is on the order of machine precision when the solution is of order P1.
 *			- This error is reduced when the TET/PYR VOLUME cubature orders are increased. Noting that the cubature
 *			  order is not seen in the L2 projection operator (as it is embedded within it), the L2 projection operators
 *			  are computed using the highest available cubature orders for TET and PYR ELEMENTs, with no effect on
 *			  performance during code execution other than in the preprocessing stage.
 *				Include higher order accurate PYR cubature nodes for the computation of L2 projection operators and
 *				check effect on error in this test. (ToBeDeleted)
 *				Note that using PIvcMaxPYR-1 (==5) gives more accurate results for P3 TET and PYR tests than using
 *				PIvcMaxPYR (==6). Reducing PIvcTET always gives less accurate results. (ToBeDeleted)
 *			- The error is also on the order of machine precision for TET6 when only a single refinement/coarsening step
 *			  is performed (No TET to PYR projections), while it is much higher for the PYR case. This suggests that the
 *			  problem is with the TET to PYR L2 projection operator. INVESTIGATE (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 */

struct S_L2proj {
	char         **argvNew, *EName, *CtrlName[2];
	int          nargc;
	unsigned int Nref, update_argv;
};

static double *get_L2err    (void);
static void   mark_VOLUMEs  (const unsigned int adapt_type);

static void test_L2_projection(struct S_L2proj *data)
{
	char *PrintName = malloc(STRLEN_MAX * sizeof *PrintName); // free

	unsigned int TETrefineType = DB.TETrefineType;

	char         *argvNew[2];
	unsigned int pass = 0, refType, NrefTypes;
	double       *L2err[2];

	argvNew[0] = data->argvNew[0];
	argvNew[1] = data->argvNew[1];

	if (strstr(data->EName,"TET"))
		NrefTypes = 3;
	else
		NrefTypes = 1;

	for (refType = 0; refType < NrefTypes; refType++) {
		// p-refinement
		strcpy(argvNew[1],data->CtrlName[0]);
		strcat(argvNew[1],"_p/Test_L2_proj_p_"); strcat(argvNew[1],data->CtrlName[1]);

		code_startup_mod_prmtrs(data->nargc,(char const *const *const) data->argvNew,0,data->update_argv,1);
		if      (refType == 0) DB.TETrefineType = TET8;
		else if (refType == 1) DB.TETrefineType = TET12;
		else if (refType == 2) DB.TETrefineType = TET6;
		code_startup_mod_prmtrs(data->nargc,(char const *const *const) data->argvNew,data->Nref,data->update_argv,2);

		L2err[0] = get_L2err();
		mark_VOLUMEs(PREFINE); mesh_update();
		mark_VOLUMEs(PCOARSE); mesh_update();
		L2err[1] = get_L2err();

		pass = 0;
		if (array_norm_diff_d(NVAR3D+1,L2err[0],L2err[1],"Inf") < 1e2*EPS)
			pass = 1;

		//     0         10        20        30        40        50
		if (strstr(data->EName,"TRI"))
			sprintf(PrintName,"L2_projections (%s%d,  ADAPT_P):",data->EName,refType);
		else
			sprintf(PrintName,"               (%s%d,  ADAPT_P):",data->EName,refType);
		test_print2(pass,PrintName);
		free(L2err[0]), free(L2err[1]);

		code_cleanup();


		// h-refinement
		strcpy(argvNew[1],data->CtrlName[0]);
		strcat(argvNew[1],"_h/Test_L2_proj_h_"); strcat(argvNew[1],data->CtrlName[1]);

		code_startup_mod_prmtrs(data->nargc,(char const *const *const) data->argvNew,0,data->update_argv,1);
		if      (refType == 0) DB.TETrefineType = TET8;
		else if (refType == 1) DB.TETrefineType = TET12;
		else if (refType == 2) DB.TETrefineType = TET6;
		code_startup_mod_prmtrs(data->nargc,(char const *const *const) data->argvNew,data->Nref,data->update_argv,2);

		L2err[0] = get_L2err();
		mark_VOLUMEs(HREFINE); mesh_update();
		mark_VOLUMEs(HCOARSE); mesh_update();
		L2err[1] = get_L2err();

		pass = 0;
		if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1.4e4*EPS) {
			pass = 1;
		} else if (array_norm_diff_d(1,L2err[0],L2err[1],"Inf") < 1e-5) {
			pass = 1;
			printf("\nWarning: h L2 projection test for P%d %ss passing with norm_diff = % .3e\n\n",
				   DB.PGlobal,data->CtrlName[1],array_norm_diff_d(1,L2err[0],L2err[1],"Inf"));
			TestDB.Nwarnings++;
		}

		test_print2(pass,"               (          ADAPT_H):");
		free(L2err[0]), free(L2err[1]);

		code_cleanup();
	}
	DB.TETrefineType = TETrefineType;

	free(PrintName);
}

void test_integration_L2_projections(int nargc, char **argv)
{
	char         **argvNew;

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

	struct S_L2proj *data;

	data = malloc(sizeof *data); // free
	data->EName       = malloc(STRLEN_MIN * sizeof *(data->EName));       // free
	data->CtrlName[0] = malloc(STRLEN_MAX * sizeof *(data->CtrlName[0])); // free
	data->CtrlName[1] = malloc(STRLEN_MAX * sizeof *(data->CtrlName[1])); // free

	strcpy(data->CtrlName[0],"test/L2_proj");
	data->argvNew     = argvNew;
	data->nargc       = nargc;
	data->Nref        = 2;
	data->update_argv = 0;

	// **************************************************************************************************** //
	// TRIs
	strcpy(data->EName,"TRI   ");
	strcpy(data->CtrlName[1],"TRI");

	test_L2_projection(data);

	// **************************************************************************************************** //
	// QUADs
	strcpy(data->EName,"QUAD  ");
	strcpy(data->CtrlName[1],"QUAD");

	test_L2_projection(data);

	// **************************************************************************************************** //
	// TETs
	printf("\nInclude higher order accurate TET/PYR cubature nodes for the computation of L2 projection operators and\n"
 	         "check effect on error in this test.\n\n");
	TestDB.Nwarnings++;

	strcpy(data->EName,"TET   ");
	strcpy(data->CtrlName[1],"TET");

	test_L2_projection(data);

	// **************************************************************************************************** //
	// HEXs
	strcpy(data->EName,"HEX   ");
	strcpy(data->CtrlName[1],"HEX");

	test_L2_projection(data);

	// **************************************************************************************************** //
	// WEDGEs
	strcpy(data->EName,"WEDGE ");
	strcpy(data->CtrlName[1],"WEDGE");

	test_L2_projection(data);

	// **************************************************************************************************** //
	// PYRs
	strcpy(data->EName,"PYR   ");
	strcpy(data->CtrlName[1],"PYR");

	test_L2_projection(data);


	free(data->EName);
	free(data->CtrlName[0]);
	free(data->CtrlName[1]);
	free(data);

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
