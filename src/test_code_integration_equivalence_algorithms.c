// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_integration_equivalence_algorithms.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"

#include "solver.h"
#include "adaptation.h"

#include "array_norm.h"
#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for equivalence between algorithms integration testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void set_test_equivalence_data(struct S_equivalence_algs *data, const char *TestName)
{
	// default values
	data->P  = 3;
	data->ML = 1;

	data->PG_add        = 1;
	data->IntOrder_add  = 0;
	data->IntOrder_mult = 2;

	data->NAlgs = 4; // All algorithms

	if (strstr(TestName,"Euler")) {
		if (strstr(TestName,"n-Cylinder_HollowSection")) {
			strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_ToBeCurvedMIXED3D_HW");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"NavierStokes")) {
		data->NAlgs = 3; // WEDGE sum factorization currently not being checked (ToBeModified)
		if (strstr(TestName,"n-Cylinder_Hollow")) {
			strcpy(data->argvNew[1],"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedMIXED2D");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void test_equivalence_algorithms(struct S_equivalence_algs *const data, char const *const TestName)
{
	/*
	 *	Expected Output:
	 *		Correspondence between RHS terms computed using various algotithm options: Standard, sparse, sum factorized.
	 */

	set_test_equivalence_data(data,TestName);

	int  const               nargc   = data->nargc;
	char const *const *const argvNew = (char const *const *const) data->argvNew;

	TestDB.PGlobal = data->P;
	TestDB.ML      = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	unsigned int const NAlgs = data->NAlgs;

	unsigned int *RHS_size = NULL, NV = 0;
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

		DB.AllowSparseVOL  = 0;
		DB.AllowSparseFACE = 0;
		if (alg == 0) { // Standard
			for (size_t P = 0; P <= DB.PMax; P++) {
			for (size_t i = 0; i < 2; i++) {
			for (size_t j = 0; j < 2; j++) {
				SF_BE[P][i][j] = 0;
			}}}

			for (size_t i = 0; i < NEC+1; i++)
				VFPartUnity[i] = 0;
		} else if (alg == 1) { // Sparse
			for (size_t P = 0; P <= DB.PMax; P++) {
			for (size_t i = 0; i < 2; i++) {
			for (size_t j = 0; j < 2; j++) {
				SF_BE[P][i][j] = 0;
			}}}

			DB.AllowSparseVOL  = 1;
			DB.AllowSparseFACE = 1;
		} else if (alg == 2) { // Sum factorized TP
			for (size_t P = 0; P <= DB.PMax; P++) {
			for (size_t i = 0; i < 2; i++) {
				SF_BE[P][0][i] = 1;
			}}

			for (size_t i = 0; i < NEC+1; i++)
				VFPartUnity[i] = 0;
		} else if (alg == 3) { // Sum factorized WEDGE
			for (size_t P = 0; P <= DB.PMax; P++) {
			for (size_t i = 0; i < 2; i++) {
				SF_BE[P][1][i] = 1;
			}}

			for (size_t i = 0; i < NEC+1; i++)
				VFPartUnity[i] = 0;
		} else {
			EXIT_UNSUPPORTED;
		}
		code_startup_mod_ctrl(nargc,argvNew,0,1,3);
		mesh_to_level(TestDB.ML);
		mesh_to_order(TestDB.PGlobal);

if (DB.Method != METHOD_DG)
	EXIT_UNSUPPORTED; // update this (ToBeDeleted)

		// Compute RHS
		struct S_solver_info solver_info = constructor_solver_info(false,false,false,'E',DB.Method);
		compute_RLHS(&solver_info);

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

		if (RHS_diff == 0.0)
			test_print_warning("Potentially comparing the same functions");
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
