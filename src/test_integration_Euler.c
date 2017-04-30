// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_code_integration_equivalence_real_complex.h"
#include "test_code_integration_linearization.h"
#include "test_code_integration_conv_order.h"
#include "test_support.h"

#include "adaptation.h"
#include "explicit_VOLUME_info.h"
#include "explicit_VOLUME_info_c.h"
#include "explicit_FACE_info.h"
#include "explicit_FACE_info_c.h"

#include "array_norm.h"
#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Test various aspects of the Euler solver implementation:
 *			- Equivalence between real and complex versions of functions;
 *			- Equivalence between running using different algorithms (with different flop counts);
 *			- Linearization;
 *			- Optimal convergence orders.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_equivalence {
	bool         TestTRI;
	char         **argvNew, *PrintName;
	unsigned int P, ML, Adapt, PG_add, IntOrder_add, IntOrder_mult;
};

static void set_test_equivalence_data(struct S_equivalence *data, const char *TestName)
{
	// default values
	data->P     = 3;
	data->ML    = 1;
	data->Adapt = ADAPT_HP;

	data->PG_add        = 1;
	data->IntOrder_add  = 0;
	data->IntOrder_mult = 2;

	if (strstr(TestName,"n-Cylinder_HollowSection")) {
		strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void test_equivalence_alg(int const nargc, char const *const *const argvNew, char const *const TestName,
                                 struct S_equivalence *const data)
{
	/*
	 *	Expected Output:
	 *		Correspondence between RHS terms computed using various algotithm options: Standard, sparse, sum factorized.
	 */

	set_test_equivalence_data(data,TestName);

	TestDB.PGlobal = data->P;
	TestDB.ML      = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	unsigned int NAlgs = 3;

	unsigned int *RHS_size, NV = 0;
	double       **RHS[NAlgs];

	for (size_t alg = 0; alg < NAlgs; alg++) {
		code_startup_mod_ctrl(nargc,argvNew,0,1,1);

		// Set necessary parameters for the sparse algorithm
		strcpy(DB.BasisType,"Nodal");
		strcpy(DB.NodeType,"GLL-AO");
		DB.Collocated = 1;

		code_startup_mod_ctrl(nargc,argvNew,0,1,2);
		unsigned int ***SF_BE     = DB.SF_BE,
		             *VFPartUnity = DB.VFPartUnity;
		if (alg == 0) { // Standard
			for (size_t P = 0; P <= DB.PMax; P++) {
			for (size_t i = 0; i < 2; i++) {
			for (size_t j = 0; j < 2; j++) {
				SF_BE[P][i][j] = 0;
			}}}

			for (size_t i = 0; i < NEC+1; i++)
				VFPartUnity[i] = 0;

			DB.AllowSparseVOL  = 0;
			DB.AllowSparseFACE = 0;
		} else if (alg == 1) { // Sparse
			for (size_t P = 0; P <= DB.PMax; P++) {
			for (size_t i = 0; i < 2; i++) {
			for (size_t j = 0; j < 2; j++) {
				SF_BE[P][i][j] = 0;
			}}}
		} else if (alg == 2) { // Sum factorized
			for (size_t P = 0; P <= DB.PMax; P++) {
			for (size_t i = 0; i < 2; i++) {
			for (size_t j = 0; j < 2; j++) {
				SF_BE[P][i][j] = 1;
			}}}

			for (size_t i = 0; i < NEC+1; i++)
				VFPartUnity[i] = 0;
		} else {
			EXIT_UNSUPPORTED;
		}
		code_startup_mod_ctrl(nargc,argvNew,0,1,3);
		mesh_to_level(TestDB.ML);
		mesh_to_order(TestDB.PGlobal);

		// Compute RHS
		explicit_VOLUME_info();
		explicit_FACE_info();

		// Copy VOLUME->RHS to RHS[alg]
		unsigned int Nvar = DB.Nvar;
		if (alg == 0)
			RHS_size = calloc(DB.NV , sizeof *RHS_size); // free
		RHS[alg] = calloc(DB.NV , sizeof *RHS[alg]); // free

		size_t IndV = 0;
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next, IndV++) {
			unsigned int NvnS = VOLUME->NvnS;

			if (alg == 0)
				RHS_size[IndV] = NvnS*Nvar;
			double *RHS_current = malloc(NvnS*Nvar * sizeof *RHS_current); // free
			for (size_t i = 0; i < NvnS*Nvar; i++)
				RHS_current[i] = VOLUME->RHS[i];

			RHS[alg][IndV] = RHS_current;
		}

		if (IndV != DB.NV)
			printf("Error: Wrong number of VOLUMEs found.\n"), EXIT_UNSUPPORTED;

		if (alg == 0) {
			NV = DB.NV;
			set_PrintName("equiv_alg",data->PrintName,&data->TestTRI);
		}
		code_cleanup();
	}

	// Check for equivalence
	unsigned int pass = 1;

	// Check that different algorithms were called (when applicable) by ensuring that some RHS terms are slightly
	// different.
	for (size_t alg = 1; alg < NAlgs; alg++) {
		double RHS_diff = 0.0;
		for (size_t IndV = 0; IndV < NV; IndV++)
			RHS_diff += array_norm_diff_d(RHS_size[IndV],RHS[0][IndV],RHS[alg][IndV],"Inf");

		if (RHS_diff == 0.0) {
			printf("\nWarning: Potentially comparing the same functions.\n\n");
			TestDB.Nwarnings++;
		}
	}

	for (size_t IndV = 0; IndV < NV; IndV++) {
		double RHS_diff = 0.0;
		for (size_t alg = 1; alg < NAlgs; alg++)
			RHS_diff = array_norm_diff_d(RHS_size[IndV],RHS[0][IndV],RHS[alg][IndV],"Inf");

		if (RHS_diff > 1e1*EPS || array_norm_d(RHS_size[IndV],RHS[0][IndV],"Inf") < EPS) {
			printf("%3zu % .3e % .3e\n",IndV,RHS_diff,array_norm_d(RHS_size[IndV],RHS[0][IndV],"Inf"));
			pass = 0;
			break;
		}
	}
	free(RHS_size);
	for (size_t alg = 0; alg < NAlgs; alg++)
		array_free2_d(NV,RHS[alg]);

	test_print2(pass,data->PrintName);
}

void test_integration_Euler(int nargc, char **argv)
{
//	bool const (ToBeModified)
	bool RunTests_equivalence_real_complex = 1,
	     RunTests_equivalence_algorithms   = 1,
	     RunTests_linearization            = 1,
	     RunTests_conv_order               = 1;

	// ToBeDeleted after Manmeet has finished his initial verification.
	bool const PeriodicVortexOnly = 0;
	if (PeriodicVortexOnly) {
		RunTests_equivalence_real_complex = 0;
		RunTests_equivalence_algorithms   = 0;
		RunTests_linearization            = 0;
	}

	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);   // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew);  // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew);  // free
	PrintName  = malloc(STRLEN_MAX * sizeof *PrintName); // free

	// silence
	strcpy(argvNew[0],argv[0]);

	// **************************************************************************************************** //
	// Real/Complex Equivalence
	// **************************************************************************************************** //
	if (RunTests_equivalence_real_complex) {
		struct S_equivalence_rc *const data_rc = calloc(1 , sizeof *data_rc); // free
		data_rc->nargc     = nargc;
		data_rc->argvNew   = argvNew;
		data_rc->PrintName = PrintName;

		test_equivalence_real_complex(data_rc,"Euler_n-Cylinder_HollowSection_CurvedMIXED2D");

		free(data_rc);
	} else {
		test_print_warning("Euler equivalence real/complex testing currently disabled");
	}

	// **************************************************************************************************** //
	// Algorithm equivalence
	// **************************************************************************************************** //
	if (RunTests_equivalence_algorithms) {
		struct S_equivalence *const data_alg = calloc(1 , sizeof *data_alg); // free

		data_alg->argvNew   = argvNew;
		data_alg->PrintName = PrintName;

		test_equivalence_alg(nargc,(char const *const *const) argvNew,"n-Cylinder_HollowSection_CurvedMIXED2D",data_alg);

		// Using a mixed HEX-WEDGE mesh would check everything here (ToBeModified)
		test_print_warning("Equivalence of WEDGE sum factorized computation is currently not being checked");

		free(data_alg);
	} else {
		test_print_warning("Euler equivalence algorithms testing currently disabled");
	}

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (RunTests_linearization) {
		struct S_linearization *const data_l = calloc(1 , sizeof *data_l); // free

		data_l->nargc     = nargc;
		data_l->argvNew   = argvNew;
		data_l->PrintName = PrintName;

		test_linearization(data_l,"Euler_MIXED2D");
bool const test_3D = 1;
if (test_3D) {
		test_linearization(data_l,"Euler_MIXED_TET_PYR");
		test_linearization(data_l,"Euler_MIXED_HEX_WEDGE");
} else {
		test_print_warning("3D Euler linearization testing currently disabled");
}

		free(data_l);
	} else {
		test_print_warning("Euler linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order
	// **************************************************************************************************** //
	if (RunTests_conv_order) {
		struct S_convorder *const data_c = calloc(1 , sizeof *data_c); // free

		data_c->nargc     = nargc;
		data_c->argvNew   = argvNew;
		data_c->PrintName = PrintName;

if (!PeriodicVortexOnly) {
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_CurvedMIXED2D");
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedMIXED2D");
}

bool const test_3D = 0;
if (test_3D) {
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedTET");
		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedHEX"); // ~Optimal
//		test_conv_order(data_c,"Euler_n-Cylinder_HollowSection_ToBeCurvedWEDGE"); // Need to implement operators
} else {
		test_print_warning("3D SupersonicVortex testing is currently disabled");
}

if (PeriodicVortexOnly)
		test_conv_order(data_c,"Euler_n-Cube_QUAD");

		test_print_warning("Add all integration tests for PeriodicVortex case (Stationary and moving)");

		free(data_c);
	} else {
		test_print_warning("Euler convergence order testing currently disabled");
	}

	array_free2_c(2,argvNew);
	free(PrintName);
}
