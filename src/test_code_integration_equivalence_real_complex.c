// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_integration_equivalence_real_complex.h"

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

#include "adaptation.h"
#include "explicit_GradW.h"
#include "explicit_GradW_c.h"
#include "explicit_VOLUME_info.h"
#include "explicit_VOLUME_info_c.h"
#include "explicit_FACE_info.h"
#include "explicit_FACE_info_c.h"
#include "finalize_RHS.h"
#include "finalize_RHS_c.h"

#include "array_norm.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for equivalence between real and complex integration testing.
 *
 *	Comments:
 *		All supported numerical fluxes are tested.
 *		The initialized What is randomly perturbed such that boundary conditions are not exactly satisfied, leading to
 *		zero jumps in FACE->Qhat terms.
 *		The contributions from the VOLUME and FACE terms to RHS can be isolated by enabling/disabling the VOLUMEOnly
 *		flag.
 *
 *	Notation:
 *
 *	References:
 */

static void set_test_equivalence_data(struct S_equivalence_rc *data, const char *TestName)
{
	// default values
	data->P  = 3;
	data->ML = 1;

	data->PG_add        = 1;
	data->IntOrder_add  = 0;
	data->IntOrder_mult = 2;

	if (strstr(TestName,"Euler")) {
		data->NNumFluxes = TEST_N_INVISCID_FLUXES;
		for (size_t i = 0; i < TEST_N_INVISCID_FLUXES; i++)
			TestDB.EnteredInviscidFlux[i] = 0;
		data->EnteredNumFluxes = TestDB.EnteredInviscidFlux;

		if (strstr(TestName,"n-Cylinder_HollowSection")) {
			strcpy(data->argvNew[1],"test/Euler/Test_Euler_SupersonicVortex_CurvedMIXED2D");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"NavierStokes")) {
		data->NNumFluxes = TEST_N_VISCOUS_FLUXES;
		for (size_t i = 0; i < TEST_N_VISCOUS_FLUXES; i++)
			TestDB.EnteredViscousFlux[i] = 0;
		data->EnteredNumFluxes = TestDB.EnteredViscousFlux;

		if (strstr(TestName,"n-Cylinder_Hollow")) {
			if (strstr(TestName,"CurvedTRI")) {
				strcpy(data->argvNew[1],"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedTRI");
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void test_equivalence_real_complex(struct S_equivalence_rc *const data, char const *const TestName)
{
	/*
	 *	Expected Output:
	 *		Correspondence between RHS terms computed using real and complex functions.
	 */

	bool const VOLUMEOnly = 0;

	set_test_equivalence_data(data,TestName);

	int  const               nargc   = data->nargc;
	char const *const *const argvNew = (char const *const *const) data->argvNew;

	TestDB.PGlobal       = data->P;
	TestDB.ML            = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	unsigned int const NNumFluxes = data->NNumFluxes;

	unsigned int pass = 1;
	for (size_t nNF = 0; nNF < NNumFluxes; nNF++) {
		code_startup_mod_prmtrs(nargc,(char const *const *const) argvNew,0,2,1);

		if (strstr(TestName,"Euler")) {
			if      (nNF == 0) DB.InviscidFluxType = FLUX_LF;
			else if (nNF == 1) DB.InviscidFluxType = FLUX_ROE;
			else               EXIT_UNSUPPORTED;
		} else if (strstr(TestName,"NavierStokes")) {
			if      (nNF == 0) DB.ViscousFluxType = FLUX_CDG2;
			else if (nNF == 1) DB.ViscousFluxType = FLUX_BR2;
			else               EXIT_UNSUPPORTED;
		} else {
			EXIT_UNSUPPORTED;
		}
		code_startup_mod_prmtrs(nargc,(char const *const *const) argvNew,0,2,2);
		mesh_to_level(TestDB.ML);
		mesh_to_order(TestDB.PGlobal);

		unsigned int const Nvar = DB.Nvar;

		// Copy What to What_c
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int const NvnS = VOLUME->NvnS;

			if (VOLUME->What_c)
				free(VOLUME->What_c);
			VOLUME->What_c = malloc(NvnS*Nvar * sizeof *(VOLUME->What_c));

			for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++) {
				VOLUME->What[i] += 1e3*EPS*((double) rand() / ((double) RAND_MAX+1));
				VOLUME->What_c[i] = VOLUME->What[i];
			}
		}

		// Compute RHS terms using the real and complex functions
		explicit_GradW();
		explicit_VOLUME_info();
		explicit_FACE_info();
		if (!VOLUMEOnly)
			finalize_RHS();

		explicit_GradW_c();
		explicit_VOLUME_info_c();
		explicit_FACE_info_c();
		if (!VOLUMEOnly)
			finalize_RHS_c();

		// Check for equivalence
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int const NvnS = VOLUME->NvnS;

			if (array_norm_diff_dc(NvnS*Nvar,VOLUME->RHS,VOLUME->RHS_c,"Inf") > 1e1*EPS ||
				array_norm_d(NvnS*Nvar,VOLUME->RHS,"Inf") < EPS) {
				array_print_d(NvnS,Nvar,VOLUME->RHS,'C');
				array_print_cmplx(NvnS,Nvar,VOLUME->RHS_c,'C');
			    printf("% .3e % .3e\n",array_norm_diff_dc(NvnS*Nvar,VOLUME->RHS,VOLUME->RHS_c,"Inf"),
				                       array_norm_d(NvnS*Nvar,VOLUME->RHS,"Inf"));
				pass = 0;
				break;
			}
		}

		if (nNF == 0)
			set_PrintName("equiv_rc",data->PrintName,&data->TestTRI);
		code_cleanup();
	}

	// Check that all numerical fluxes were entered
	for (size_t i = 0; i < data->NNumFluxes; i++) {
		if (data->EnteredNumFluxes[i] == 0) {
			printf("First unchecked numerical flux index: %zu.\n",i);
			pass = 0;
			break;
		}
	}

	test_print2(pass,data->PrintName);
}
